#include <compressor/SZProgressiveInterpolationCompressorV4.hpp>
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
