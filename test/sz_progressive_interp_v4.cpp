#include <compressor/SZProgressiveInterpolationCompressorV4.hpp>
#include <compressor/SZBlockInterpolationCompressor.hpp>
#include <compressor/SZInterpolationCompressor.hpp>
#include <quantizer/IntegerQuantizer.hpp>
#include <predictor/ComposedPredictor.hpp>
#include <predictor/SimplePredictor.hpp>
#include <lossless/Lossless_zstd.hpp>
#include <meta/meta_compress.hpp>
#include <utils/Iterator.hpp>
#include <utils/Verification.hpp>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <memory>
#include <type_traits>
#include <sstream>

template<uint N, class ... Dims>
void interp_compress_decompress(const char *path, double eb, int interp_op, int direction_op,
                                int level_independent, int block_size, int level_fill, Dims ... args) {
    std::string ori_datafile(path);
    std::stringstream ss;
    ss << ori_datafile.substr(ori_datafile.rfind('/') + 1) << ".sz3.out";
    std::string decompressed_file_name = ss.str();
    ss.str(std::string());
    ss << ori_datafile.substr(ori_datafile.rfind('/') + 1) << ".sz3";
    std::string compress_file_name = ss.str();
    std::cout << "decompressed file = " << decompressed_file_name << std::endl;

    std::vector<size_t> compressed_size;
    double compression_ratio;
    auto dims = std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...};
    size_t num = std::accumulate(dims.begin(), dims.end(), (size_t) 1, std::multiplies<size_t>());
    {

        std::cout << "****************** compression ****************" << std::endl;
        std::cout << "Interp op          = " << interp_op << std::endl
                  << "Direction          = " << direction_op << std::endl
                  << "Level independent  = " << level_independent << std::endl
                  << "Block size         = " << block_size << std::endl
                  << "Level fill         = " << level_fill << std::endl;
        struct timespec start, end;
        clock_gettime(CLOCK_REALTIME, &start);


        auto sz = SZ::SZProgressiveInterpolationCompressorV4<float, N, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer<float>(eb),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd(),
                dims, interp_op, direction_op, 32,
                level_independent, block_size, level_fill
        );
        sz.compress(path, compress_file_name.data(), compressed_size);
        size_t total_compressed_size = std::accumulate(compressed_size.begin(), compressed_size.end(), (size_t) 0);

        clock_gettime(CLOCK_REALTIME, &end);
        double compress_time =
                (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
        compression_ratio = num * sizeof(float) * 1.0 / total_compressed_size;
        std::cout << "Compression time = " << compress_time << "s" << std::endl;
        std::cout << "Compressed size = " << total_compressed_size << std::endl;
        std::cout << "Compression ratio = " << compression_ratio << std::endl << std::endl;
    }
    {
        std::cout << "****************** Decompression ****************" << std::endl;

        struct timespec start, end;
        clock_gettime(CLOCK_REALTIME, &start);

        auto dims = std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...};
        auto sz = SZ::SZProgressiveInterpolationCompressorV4<float, N, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer<float>(eb),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd(),
                dims, interp_op, direction_op, 32,
                level_independent, block_size, level_fill
        );
        sz.decompress(compress_file_name.data(), compressed_size, decompressed_file_name.c_str());

        clock_gettime(CLOCK_REALTIME, &end);
        double decompress_time =
                (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
        std::cout << "Decompression time = " << decompress_time << "s" << std::endl;

    }

    double psnr, nrmse;
    SZ::verify<float>(path, decompressed_file_name.c_str(), num, psnr, nrmse);

    printf("PSNR = %f, NRMSE = %.10G, Compression Ratio = %.2f\n", psnr, nrmse, compression_ratio);
    remove(decompressed_file_name.c_str());
    remove(compress_file_name.c_str());

}

template<uint N>
double interp_compress_test_block(float *data, size_t num, double eb, int interp_level, int interp_op, int direction_op,
                                  int block_size, int interp_block_size, std::array<size_t, N> dims) {

    std::cout << "****************** compression ****************" << std::endl;
    std::cout << "Interp Level       = " << interp_level << std::endl
              << "Interp Op          = " << interp_op << std::endl
              << "Direction          = " << direction_op << std::endl
              << "SZ block size      = " << block_size << std::endl
              << "Interp block size  = " << interp_block_size << std::endl;

    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);

    std::vector<float> data1(data, data + num);
    size_t compressed_size = 0;

    SZ::Config<float, N> conf(eb, dims);
    conf.block_size = block_size;
    conf.stride = conf.block_size;
    auto sz = SZ::make_sz_fast_block_interpolation_compressor(
            conf,
            SZ::SimplePredictor<float, N>(eb),
            SZ::LinearQuantizer<float>(eb),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd(),
            interp_op,
            direction_op,
            interp_level
    );

    sz.compress(data1.data(), compressed_size);
    clock_gettime(CLOCK_REALTIME, &end);
    double compression_time =
            (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
    auto compression_ratio = num * sizeof(float) * 1.0 / compressed_size;
    std::cout << "Compressed size = " << compressed_size << std::endl;
    std::cout << "Compression ratio = " << compression_ratio << std::endl;
    std::cout << "Tuning Interp Compression time = " << compression_time
              << " Ratio = " << compression_ratio
              << " Params = " << interp_level
              << " " << interp_op
              << " " << direction_op
              << " " << block_size
              << " " << interp_block_size
              << std::endl;

    std::cout << "****************** end ****************" << std::endl;
    return compression_ratio;
}

template<uint N>
float cal_sampling_ratio(size_t block, size_t n, size_t dmin, std::array<size_t, N> dims) {
    size_t sample_n = 1;
    for (auto dim: dims) {
        sample_n *= dim / dmin * 2 * block;
    }
    return sample_n * 1.0 / n;
}

template<uint N>
std::pair<int, int> interp_tuning(char *path, double eb, std::array<size_t, N> dims) {
    std::cout << "================================ BEGIN TUNING ================================" << std::endl;


    size_t num = 0;
    auto data = SZ::readfile<float>(path, num);
    std::cout << "Read " << num << " elements\n";

    size_t dmin = *std::min_element(dims.begin(), dims.end());
    int sampling_block = dmin;
    while (cal_sampling_ratio<3>(sampling_block, num, dmin, dims) > 0.035) {
        sampling_block--;
    }
    if (sampling_block * 2 > dmin) {
        sampling_block = dmin / 2;
    }
    size_t b1 = dims[0] / dmin;
    size_t b2 = dims[1] / dmin;
    size_t b3 = dims[2] / dmin;
    size_t r1 = b1 * 2 * sampling_block;
    size_t r2 = b2 * 2 * sampling_block;
    size_t r3 = b3 * 2 * sampling_block;
    size_t di, dj, dk;
    size_t sampling_num = r1 * r2 * r3;
    std::array<size_t, 3> sampling_dims{r1, r2, r3};
    std::vector<float> sampling_data(sampling_num, 0);
    for (size_t bi = 0; bi < b1; bi++) {
        for (size_t bj = 0; bj < b2; bj++) {
            for (size_t bk = 0; bk < b3; bk++) {
                for (size_t i = 0; i < 2 * sampling_block; i++) {
                    for (size_t j = 0; j < 2 * sampling_block; j++) {
                        for (size_t k = 0; k < 2 * sampling_block; k++) {
                            di = i < sampling_block ? i + sampling_block : dmin - 3 * sampling_block + i;
                            dj = j < sampling_block ? j + sampling_block : dmin - 3 * sampling_block + j;
                            dk = k < sampling_block ? k + sampling_block : dmin - 3 * sampling_block + k;
                            auto d = data[(bi * dmin + di) * dims[1] * dims[2] + (bj * dmin + dj) * dims[2] +
                                          bk * dmin + dk];
                            sampling_data[(bi * 2 * sampling_block + i) * r2 * r3
                                          + (bj * 2 * sampling_block + j) * r3
                                          + bk * 2 * sampling_block + k] = d;
                        }
                    }
                }
            }
        }
    }
    printf("Tuning sampling block = %d percent = %.3f%% Time = %.3f \n",
           sampling_block, sampling_num * 100.0 / num, 0.0);

    int interp_level = -1, interp_op, direction_op = 0, block_size = sampling_block, interp_block_size = sampling_block;

    double best_interp_ratio = 0, ratio;
    for (int i = 0; i < 2; i++) {
        ratio = interp_compress_test_block<N>(sampling_data.data(), sampling_num, eb, interp_level, i, direction_op,
                                              block_size, interp_block_size, sampling_dims);
        if (ratio > best_interp_ratio) {
            best_interp_ratio = ratio;
            interp_op = i;
        }
    }
    std::cout << "Tuning interp best interp_op = " << interp_op << " , best ratio = " << best_interp_ratio << std::endl;

    ratio = interp_compress_test_block<N>(sampling_data.data(), sampling_num, eb, interp_level, interp_op, 5,
                                          block_size, interp_block_size, sampling_dims);
    if (ratio > best_interp_ratio * 1.02) {
        direction_op = 5;
    }
    std::cout << "Tuning interp interp_op = " << interp_op << ", direction_op = " << direction_op << std::endl;
    std::cout << "================================ END TUNING ================================" << std::endl;

    return std::pair<int, int>{interp_op, direction_op};

}

int main(int argc, char **argv) {


    int dim = atoi(argv[2] + 1);
    assert(1 <= dim && dim <= 4);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }
    float eb = atof(argv[argp++]);

    int interp_op = 0; // linear
    int direction_op = 0; // dimension high -> low
    if (argp < argc) {
        interp_op = atoi(argv[argp++]);
    }
    if (argp < argc) {
        direction_op = atoi(argv[argp++]);
    }
    if (interp_op == -1 || direction_op == -1) {
        if (dim != 3) {
            std::cout << "Tuning only support 3D array.\n";
            return 0;
        }
        std::array<size_t, 3> dims_array{dims[0], dims[1], dims[2]};
        auto x = interp_tuning<3>(argv[1], eb, dims_array);
        interp_op = x.first;
        direction_op = x.second;
    }

    int level_fill = 0, level_independent = 0;
    int block_size = 128;
    if (argp < argc) {
        level_independent = atoi(argv[argp++]);
    }

    if (argp < argc) {
        block_size = atoi(argv[argp++]);
    }

    if (argp < argc) {
        level_fill = atoi(argv[argp++]);
    }


    if (dim == 1) {
        interp_compress_decompress<1>(argv[1], eb, interp_op, direction_op, level_independent,
                                      block_size, level_fill, dims[0]);
    } else if (dim == 2) {
        interp_compress_decompress<2>(argv[1], eb, interp_op, direction_op, level_independent,
                                      block_size, level_fill, dims[0], dims[1]);
    } else if (dim == 3) {
        interp_compress_decompress<3>(argv[1], eb, interp_op, direction_op, level_independent,
                                      block_size, level_fill, dims[0], dims[1], dims[2]);
    } else if (dim == 4) {
        interp_compress_decompress<4>(argv[1], eb, interp_op, direction_op, level_independent,
                                      block_size, level_fill, dims[0], dims[1], dims[2], dims[3]);
    }


    return 0;
}
