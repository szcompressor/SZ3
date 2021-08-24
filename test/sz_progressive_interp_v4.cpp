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
void interp_compress_decompress(const char *path, double eb, int interp_level, int interp_op,
                                int direction_op, int block_size, int interp_block_size, Dims ... args) {
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
    size_t num = std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<size_t>());
    {

        std::cout << "****************** compression ****************" << std::endl;
        std::cout << "Interp Level       = " << interp_level << std::endl
                  << "Interp Op          = " << interp_op << std::endl
                  << "Direction          = " << direction_op << std::endl
                  << "SZ block size      = " << block_size << std::endl
                  << "Interp block size  = " << interp_block_size << std::endl;
        struct timespec start, end;
        clock_gettime(CLOCK_REALTIME, &start);


        auto sz = SZ::SZProgressiveInterpolationCompressorV4<float, N, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
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
                dims,
                interp_op,
                direction_op,
                32,
                interp_block_size,
                interp_level
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

    int interp_level = 1;
    int interp_op = 0; // linear
    int direction_op = 0; // dimension high -> low
    int block_size = 6;
    int interp_block_size = 32;
    if (argp < argc) {
        interp_level = atoi(argv[argp++]);
    }
    if (argp < argc) {
        interp_op = atoi(argv[argp++]);
    }
    if (argp < argc) {
        direction_op = atoi(argv[argp++]);
    }

    if (argp < argc) {
        block_size = atoi(argv[argp++]);
    }
    if (argp < argc) {
        interp_block_size = atoi(argv[argp++]);
    }

    if (dim == 1) {
        interp_compress_decompress<1>(argv[1], eb, interp_level, interp_op, direction_op, block_size,
                                      interp_block_size, dims[0]);
    } else if (dim == 2) {
        interp_compress_decompress<2>(argv[1], eb, interp_level, interp_op, direction_op, block_size,
                                      interp_block_size, dims[0], dims[1]);
    } else if (dim == 3) {
        interp_compress_decompress<3>(argv[1], eb, interp_level, interp_op, direction_op, block_size,
                                      interp_block_size, dims[0], dims[1], dims[2]);
    } else if (dim == 4) {
        interp_compress_decompress<4>(argv[1], eb, interp_level, interp_op, direction_op,
                                      block_size, interp_block_size, dims[0], dims[1], dims[2], dims[3]);
    }


    return 0;
}
