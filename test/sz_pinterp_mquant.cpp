//#include <compressor/SZProgressiveIndependentBlock.hpp>
//#include <compressor/SZProgressive.hpp>
#include <compressor/SZProgressiveMQuant.hpp>
#include <quantizer/IntegerQuantizer2.hpp>
#include <predictor/ComposedPredictor.hpp>
#include <lossless/Lossless_zstd.hpp>
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
    std::vector<size_t> compressed_size;
    size_t total_compressed_size = 0;
    SZ::uchar *compressed;

    size_t num = 0;
    auto data = SZ::readfile<float>(path, num);

    {

        std::cout << "****************** compression ****************" << std::endl;
        std::cout << "Interp op          = " << interp_op << std::endl
                  << "Direction          = " << direction_op << std::endl
                  << "Progressive Level  = " << level_independent << std::endl
                  << "Block size         = " << block_size << std::endl
                  << "Level fill         = " << level_fill << std::endl;

        SZ::Timer timer(true);
        auto dims = std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...};
        auto sz = SZ::SZProgressiveMQuant<float, N, SZ::LinearQuantizer2<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer2<float>(eb),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd(2),
                dims, interp_op, direction_op, 50000, level_independent, block_size, level_fill
        );
        compressed = sz.compress(data.get(), compressed_size);
        timer.stop("Compression");

        total_compressed_size = std::accumulate(compressed_size.begin(), compressed_size.end(), (size_t) 0);
        auto compression_ratio = num * sizeof(float) * 1.0 / total_compressed_size;
        std::cout << "Compressed size = " << total_compressed_size << std::endl;
        std::cout << "Compression ratio = " << compression_ratio << std::endl << std::endl;
    }
//    std::string file = std::string(path).substr(std::string(path).rfind('/') + 1) + ".sz3";
//    std::cout << "compressed file = " << file << std::endl;
//    SZ::writefile(file.c_str(), compressed, total_compressed_size);

    data = SZ::readfile<float>(path, num);
    {
        std::cout << "****************** Decompression ****************" << std::endl;

        float *dec_data;
        SZ::Timer timer(true);
        auto dims = std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...};
        auto sz = SZ::SZProgressiveMQuant<float, N, SZ::LinearQuantizer2<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer2<float>(eb),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd(),
                dims, interp_op, direction_op, 50000, level_independent, block_size, level_fill
        );
        dec_data = sz.decompress(compressed, compressed_size, data.get());


        timer.stop("Decompression");

//        std::string file = std::string(path).substr(std::string(path).rfind('/') + 1) + ".sz3.out";
//        std::cout << "decompressed file = " << file << std::endl;
//        SZ::writefile(file.c_str(), dec_data.get(), num);

        if (level_independent <= 0) {
            size_t num1 = 0;
            auto ori_data = SZ::readfile<float>(path, num1);
            assert(num1 == num);
            double psnr, nrmse;
            SZ::verify<float>(ori_data.get(), dec_data, num, psnr, nrmse);
            delete[]dec_data;
            delete[]compressed;
        }
//        std::vector<float> error(num);
//        for (size_t i = 0; i < num; i++) {
//            error[i] = ori_data[i] - dec_data[i];
//        }
//        std::string error_file(path);
//        error_file += ".error";
//        SZ::writefile(error_file.c_str(), error.data(), num);
//        auto compression_ratio = num * sizeof(float) * 1.0 / total_compressed_size;
//        printf("PSNR = %f, NRMSE = %.10G, Compression Ratio = %.2f\n", psnr, nrmse, compression_ratio);
    }

}


int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "psz usage: " << argv[0] <<
                  " data_file -num_dim dim0 .. dimn abs_eb [interp_op direction_op p_level block_size]"
                  << std::endl
                  << "example: " << argv[0] <<
                  " qmcpack.dat -3 33120 69 69 1e-3 [1 0 0 128]" << std::endl;
        return 0;
    }

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
        std::cout << "Tuning not support.\n";
        return 0;
    }

    int level_progressive = 0;
    if (argp < argc) {
        level_progressive = atoi(argv[argp++]);
    }

    int block_size = 128;
    if (argp < argc) {
        block_size = atoi(argv[argp++]);
    }

    int level_fill = 0;
    if (argp < argc) {
        level_fill = atoi(argv[argp++]);
    }


    if (dim == 1) {
        interp_compress_decompress<1>(argv[1], eb, interp_op, direction_op, level_progressive,
                                      block_size, level_fill, dims[0]);
    } else if (dim == 2) {
        interp_compress_decompress<2>(argv[1], eb, interp_op, direction_op, level_progressive,
                                      block_size, level_fill, dims[0], dims[1]);
    } else if (dim == 3) {
        interp_compress_decompress<3>(argv[1], eb, interp_op, direction_op, level_progressive,
                                      block_size, level_fill, dims[0], dims[1], dims[2]);
    } else if (dim == 4) {
        interp_compress_decompress<4>(argv[1], eb, interp_op, direction_op, level_progressive,
                                      block_size, level_fill, dims[0], dims[1], dims[2], dims[3]);
    }


    return 0;
}
