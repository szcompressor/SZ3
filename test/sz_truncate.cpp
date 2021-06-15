#include "compressor/SZTruncateCompressor.hpp"
#include "lossless/Lossless_zstd.hpp"
#include "lossless/Lossless_bypass.hpp"
#include "utils/FileUtil.hpp"
#include "utils/Config.hpp"
#include "utils/Verification.hpp"
#include "utils/Timer.hpp"
#include "def.hpp"
#include <cstdio>
#include <iostream>
#include <memory>
#include <sstream>
#include <random>

std::string src_file_name;
float relative_error_bound = 0;
float compression_time = 0;
int byteLen = 2;

template<typename T, class Lossless, uint N>
float SZ_compress(std::unique_ptr<T[]> const &data,
                  const SZ::Config<T, N> &conf, Lossless lossless) {

    std::cout << "****************** Options ********************" << std::endl;
    std::cout << "dimension = " << N
              << ", error bound = " << conf.eb
              << ", byteLen = " << byteLen
              << ", lossless = " << conf.lossless_op
              << std::endl;

    std::vector<T> data_ = std::vector<T>(data.get(), data.get() + conf.num);

    auto sz = make_sz_truncate_compressor(conf, lossless, byteLen);

    SZ::Timer timer;
    timer.start();
    std::cout << "****************** Compression ******************" << std::endl;


    size_t compressed_size = 0;
    std::unique_ptr<SZ::uchar[]> compressed;
    compressed.reset(sz.compress(data.get(), compressed_size));

    compression_time = timer.stop("Compression");

    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compression Ratio = " << ratio << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, 10000);
    std::stringstream ss;
    ss << src_file_name.substr(src_file_name.rfind('/') + 1)
       << "." << relative_error_bound << "." << dis(gen) << ".sz3";
    auto compressed_file_name = ss.str();
    SZ::writefile(compressed_file_name.c_str(), compressed.get(), compressed_size);
    std::cout << "Compressed file = " << compressed_file_name << std::endl;

    std::cout << "****************** Decompression ****************" << std::endl;
    compressed = SZ::readfile<SZ::uchar>(compressed_file_name.c_str(), compressed_size);

    timer.start();
    T *dec_data = sz.decompress(compressed.get(), compressed_size);
    timer.stop("Decompression");

    SZ::verify<T>(data_.data(), dec_data, conf.num);

    remove(compressed_file_name.c_str());
//    auto decompressed_file_name = compressed_file_name + ".out";
//    SZ::writefile(decompressed_file_name.c_str(), dec_data, conf.num);
//    std::cout << "Decompressed file = " << decompressed_file_name << std::endl;

    delete[] dec_data;
    return ratio;
}

template<class T, uint N>
float SZ_compress_parse_args(int argc, char **argv, int argp, std::unique_ptr<T[]> &data, float eb,
                             std::array<size_t, N> dims) {
    SZ::Config<float, N> conf(eb, dims);
    if (argp < argc) {
        conf.lossless_op = atoi(argv[argp++]);
    }
    if (conf.lossless_op > 0) {
        return SZ_compress<T>(data, conf, SZ::Lossless_zstd(conf.lossless_op));
    } else {
        return SZ_compress<T>(data, conf, SZ::Lossless_bypass());
    }
}


int main(int argc, char **argv) {
    if (argc < 4) {
        std::cout << "usage: " << argv[0] << " data_file -num_dim dim0 .. dimn byteLens " << std::endl;
        std::cout << "example: " << argv[0] << " qmcpack.dat -3 33120 69 69 2" << std::endl;
        return 0;
    }

    size_t num = 0;
    auto data = SZ::readfile<float>(argv[1], num);
    src_file_name = argv[1];
    std::cout << "Read " << num << " elements\n";

    int dim = atoi(argv[2] + 1);
    assert(1 <= dim && dim <= 4);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }
    if (argp < argc) {
        byteLen = atoi(argv[argp++]);
    }

    float eb = 0;


    if (dim == 1) {
        SZ_compress_parse_args<float, 1>(argc, argv, argp, data, eb, std::array<size_t, 1>{dims[0]});
    } else if (dim == 2) {
        SZ_compress_parse_args<float, 2>(argc, argv, argp, data, eb, std::array<size_t, 2>{dims[0], dims[1]});
    } else if (dim == 3) {
        SZ_compress_parse_args<float, 3>(argc, argv, argp, data, eb,
                                         std::array<size_t, 3>{dims[0], dims[1], dims[2]});
    } else if (dim == 4) {
        SZ_compress_parse_args<float, 4>(argc, argv, argp, data, eb,
                                         std::array<size_t, 4>{dims[0], dims[1], dims[2], dims[3]});
    }


    return 0;
}