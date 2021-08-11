#include <compressor/SZProgressiveInterpolationCompressorV3.hpp>
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

enum SZ_Option {
    SZ_LR, // SZ with Lorenzo + Regression
    SZ_L,  // SZ with 1D lorenzo (previous datapoint as predicted value)
};

template<uint N, class ... Dims>
META::meta_compress_info
interp_compress_decompress(char *path, float *data, size_t num, double eb, int interp_level, int interp_op,
                           int direction_op,
                           int block_size, int interp_block_size, SZ_Option sz_op, Dims ... args) {
    std::string compressed_file_name(path);
    META::meta_compress_info compressInfo;

    std::vector<size_t> compressed_size;
    size_t total_compressed_size = 0;
    SZ::uchar *compressed;

    {

        std::cout << "****************** compression ****************" << std::endl;
        std::cout << "Interp Level       = " << interp_level << std::endl
                  << "Interp Op          = " << interp_op << std::endl
                  << "Direction          = " << direction_op << std::endl
                  << "SZ block size      = " << block_size << std::endl
                  << "Interp block size  = " << interp_block_size << std::endl
                  << "SZ_mode            = " << sz_op << std::endl;
        struct timespec start, end;
        clock_gettime(CLOCK_REALTIME, &start);


        auto dims = std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...};
        auto sz = SZ::SZProgressiveInterpolationCompressorV3<float, N, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer<float>(eb),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd(),
                dims,
                interp_op,
                direction_op,
                32,
                interp_block_size,
                interp_level
        );
        compressed = sz.compress(data, compressed_size);
        total_compressed_size = std::accumulate(compressed_size.begin(), compressed_size.end(), (size_t) 0);

        clock_gettime(CLOCK_REALTIME, &end);
        double compress_time =
                (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
        auto compression_ratio = num * sizeof(float) * 1.0 / total_compressed_size;
        compressInfo.compress_time = compress_time;
        std::cout << "Compression time = " << compress_time << "s" << std::endl;
        std::cout << "Compressed size = " << total_compressed_size << std::endl;
        std::cout << "Compression ratio = " << compression_ratio << std::endl << std::endl;
//        exit(0);
    }
    {
        std::cout << "****************** Decompression ****************" << std::endl;

        struct timespec start, end;
        clock_gettime(CLOCK_REALTIME, &start);
        float *dec_data;

        auto dims = std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...};
        auto sz = SZ::SZProgressiveInterpolationCompressorV3<float, N, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer<float>(eb),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd(),
                dims,
                interp_op,
                direction_op,
                32,
                interp_block_size,
                interp_level
        );
        dec_data = sz.decompress(compressed, compressed_size);


        clock_gettime(CLOCK_REALTIME, &end);
        double decompress_time =
                (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
        compressInfo.decompress_time = decompress_time;
        std::cout << "Decompression time = " << decompress_time << "s" << std::endl;

//        std::string decompressed_file_name(path);
//        std::stringstream ss;
//        ss << decompressed_file_name.substr(decompressed_file_name.rfind('/') + 1)
//           << ".sz3.out";
//        decompressed_file_name = ss.str();
//        std::cout << "DEBUG decompressed file = " << decompressed_file_name << std::endl;
//        SZ::writefile(decompressed_file_name.c_str(), dec_data.get(), num);

        size_t num1 = 0;
        auto ori_data = SZ::readfile<float>(path, num1);
        assert(num1 == num);
        double psnr, nrmse;
        SZ::verify<float>(ori_data.get(), dec_data, num, psnr, nrmse);

//        std::vector<float> error(num);
//        for (size_t i = 0; i < num; i++) {
//            error[i] = ori_data[i] - dec_data[i];
//        }
//        std::string error_file(path);
//        error_file += ".error";
//        SZ::writefile(error_file.c_str(), error.data(), num);
        auto compression_ratio = num * sizeof(float) * 1.0 / total_compressed_size;
        printf("PSNR = %f, NRMSE = %.10G, Compression Ratio = %.2f\n", psnr, nrmse, compression_ratio);
        std::cout << "FINAL " << compression_ratio
                  << " " << interp_level
                  << " " << interp_op
                  << " " << direction_op
                  << " " << block_size
                  << " " << interp_block_size
                  << " " << sz_op << std::endl;
        compressInfo.psnr = psnr;
        compressInfo.nrmse = nrmse;
        compressInfo.ratio = compression_ratio;
        return compressInfo;

    }

}


template<uint N, class ... Dims>
double interp_compress_test(float *data, size_t num, double eb, int interp_level, int interp_op, int direction_op,
                            int block_size,
                            int interp_block_size, SZ_Option sz_op, Dims ... args) {

    std::cout << "****************** compression ****************" << std::endl;
    std::cout << "Interp Level       = " << interp_level << std::endl
              << "Interp Op          = " << interp_op << std::endl
              << "Direction          = " << direction_op << std::endl
              << "SZ block size      = " << block_size << std::endl
              << "Interp block size  = " << interp_block_size << std::endl
              << "SZ_mode            = " << sz_op << std::endl;

    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);

    std::vector<float> data1(data, data + num);
    size_t compressed_size = 0;
    std::unique_ptr<unsigned char[]> compressed;

    auto dims = std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...};
    auto sz = SZ::SZInterpolationCompressor<float, N, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            SZ::LinearQuantizer<float>(eb),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd(),
            dims,
            interp_block_size,
            interp_op,
            direction_op,
            interp_level
    );
    compressed.reset(sz.compress(data1.data(), compressed_size));

    clock_gettime(CLOCK_REALTIME, &end);
    double compression_time =
            (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
    auto compression_ratio = num * sizeof(float) * 1.0 / compressed_size;
    std::cout << "Compressed size = " << compressed_size << std::endl;
    std::cout << "Compression ratio = " << compression_ratio << std::endl;
    std::cout << "TUNING Interp Compression time = " << compression_time
              << " Ratio = " << compression_ratio
              << " Params = " << interp_level
              << " " << interp_op
              << " " << direction_op
              << " " << block_size
              << " " << interp_block_size
              << " " << sz_op << std::endl;

    std::cout << "****************** end ****************" << std::endl;
    return compression_ratio;
}

template<uint N, class ... Dims>
META::meta_compress_info interp_tuning(char *path, double reb, Dims ... args) {
    size_t num = 0;
    std::string compressed_file_name(path);

    auto data = SZ::readfile<float>(path, num);
    std::cout << "Read " << num << " elements\n";
    float max = data[0];
    float min = data[0];
    for (int i = 1; i < num; i++) {
        if (max < data[i]) max = data[i];
        if (min > data[i]) min = data[i];
    }
    float eb = reb * (max - min);
    int interp_level = -1;
    int interp_op, direction_op = 0, block_size = 6, interp_block_size = 32;
    SZ_Option sz_op = SZ_LR;

    double best_ratio = 0, ratio;
    for (int i = 0; i < 2; i++) {
        ratio = interp_compress_test<N>(data.get(), num, eb, interp_level, i, direction_op, block_size,
                                        interp_block_size, sz_op,
                                        args...);
        if (ratio > best_ratio) {
            best_ratio = ratio;
            interp_op = i;
        }
    }
    std::cout << "TUNING best interp_op = " << interp_op << " , best ratio = " << best_ratio << std::endl;

    int direction_op_reverse = std::tgamma(N + 1) - 1;
    ratio = interp_compress_test<N>(data.get(), num, eb, interp_level, interp_op, direction_op_reverse, block_size,
                                    interp_block_size,
                                    sz_op,
                                    args...);
    if (ratio > best_ratio * 1.0) {
        best_ratio = ratio;
        direction_op = direction_op_reverse;
    }
    std::cout << "TUNING best direction_op = " << direction_op << " , best ratio = " << best_ratio << std::endl;

    return interp_compress_decompress<N>(path, data.get(), num, eb, interp_level, interp_op, direction_op, block_size,
                                         interp_block_size, sz_op,
                                         args...);
}


int main(int argc, char **argv) {


    int dim = atoi(argv[2] + 1);
    assert(1 <= dim && dim <= 4);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }
    float reb = atof(argv[argp++]);
    if (argp >= argc) {
        if (dim == 1) {
            interp_tuning<1>(argv[1], reb, dims[0]);
        } else if (dim == 2) {
            interp_tuning<2>(argv[1], reb, dims[0], dims[1]);
        } else if (dim == 3) {
            interp_tuning<3>(argv[1], reb, dims[0], dims[1], dims[2]);
        } else if (dim == 4) {
            interp_tuning<4>(argv[1], reb, dims[0], dims[1], dims[2], dims[3]);
        }
        return 0;
    }

    int interp_level = 1;
    int interp_op = 0; // linear
    int direction_op = 0; // dimension high -> low
    if (argp < argc) {
        interp_level = atoi(argv[argp++]);
    }
    if (argp < argc) {
        interp_op = atoi(argv[argp++]);
    }
    if (argp < argc) {
        direction_op = atoi(argv[argp++]);
    }

    int block_size = 6;
    int interp_block_size = 32;
    SZ_Option sz_op = SZ_LR;
    if (argp < argc) {
        block_size = atoi(argv[argp++]);
    }
    if (argp < argc) {
        interp_block_size = atoi(argv[argp++]);
    }
    if (argp < argc) {
        sz_op = (SZ_Option) atoi(argv[argp++]);
    }
    size_t num = 0;

    auto data = SZ::readfile<float>(argv[1], num);
    float max = data[0];
    float min = data[0];
    for (int i = 1; i < num; i++) {
        if (max < data[i]) max = data[i];
        if (min > data[i]) min = data[i];
    }
    float eb = reb * (max - min);


    if (dim == 1) {
        interp_compress_decompress<1>(argv[1], data.get(), num, eb, interp_level, interp_op, direction_op, block_size,
                                      interp_block_size, sz_op, dims[0]);
    } else if (dim == 2) {
        interp_compress_decompress<2>(argv[1], data.get(), num, eb, interp_level, interp_op, direction_op, block_size,
                                      interp_block_size, sz_op, dims[0], dims[1]);
    } else if (dim == 3) {
        interp_compress_decompress<3>(argv[1], data.get(), num, eb, interp_level, interp_op, direction_op, block_size,
                                      interp_block_size, sz_op, dims[0], dims[1], dims[2]);
    } else if (dim == 4) {
        interp_compress_decompress<4>(argv[1], data.release(), num, eb, interp_level, interp_op, direction_op,
                                      block_size, interp_block_size, sz_op, dims[0], dims[1], dims[2], dims[3]);
    }


    return 0;
}
