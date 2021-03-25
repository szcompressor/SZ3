#include <vector>
#include <sstream>
#include "utils/FileUtil.h"
#include "utils/Timer.hpp"
#include "utils/Verification.hpp"
#include "def.hpp"
#include "quantizer/IntegerQuantizer.hpp"
#include "lossless/Lossless_zstd.hpp"
#include "lossless/Lossless_bypass.hpp"
#include "encoder/HuffmanEncoder.hpp"
#include "encoder/BypassEncoder.hpp"
#include "compressor/SZExaaltCompressor.hpp"
#include "utils/KmeansUtil.h"
#include "utils/QuantOptimizatioin.hpp"

std::string src_file_name;
float relative_error_bound = 0;


template<typename T, uint N>
float
SZ_Compress(std::unique_ptr<T[]> const &data, SZ::Config<T, N> conf) {

    float level_start, level_offset;
    int level_num;
    SZ::get_cluster(data.get(), conf.num, level_start, level_offset, level_num, 0.001);
    printf("start = %.3f , level_offset = %.3f, nlevel=%d\n", level_start, level_offset, level_num);

    //    level_start = -58.291; //trinity-110x
//    level_offset = 2.241; //trinity-110x
//    level_start = 0;
//    level_offset = 1.961;
    std::vector<T> data_(data.get(), data.get() + conf.num);

    SZ::Timer timer;
    std::cout << "****************** Compression ******************" << std::endl;

    timer.start();
    conf.quant_state_num = 1024;
    auto sz = SZ::SZ_Exaalt_Compressor(conf, SZ::LinearQuantizer<T>(conf.eb, conf.quant_state_num / 2),
                                       SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd(), 0);
    sz.set_level(level_start, level_offset, level_num);

    size_t compressed_size = 0;
    std::unique_ptr<SZ::uchar[]> compressed;
    compressed.reset(sz.compress(data.get(), compressed_size));
    timer.stop("Compression");

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
    std::unique_ptr<T[]> dec_data;
    dec_data.reset(sz.decompress(compressed.get(), compressed_size));
    timer.stop("Decompression");
    SZ::verify<T>(data_.data(), dec_data.get(), conf.num);

    auto decompressed_file_name = compressed_file_name + ".out";
//    SZ::writefile(decompressed_file_name.c_str(), dec_data.get(), conf.num);
//    std::cout << "Decompressed file = " << decompressed_file_name << std::endl;

    return ratio;
}

int main(int argc, char **argv) {


    size_t num = 0;
    auto data = SZ::readfile<float>(argv[1], num);
    src_file_name = argv[1];
    std::cout << "Read " << num << " elements\n";

    int dim = atoi(argv[2] + 1);
    assert(1 <= dim && dim <= 2);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }

    float max = data[0];
    float min = data[0];
    for (int i = 1; i < num; i++) {
        if (max < data[i]) max = data[i];
        if (min > data[i]) min = data[i];
    }
    char *eb_op = argv[argp++] + 1;
    float eb = 0;
    if (*eb_op == 'a') {
        eb = atof(argv[argp++]);
        relative_error_bound = eb / (max - min);
    } else {
        relative_error_bound = atof(argv[argp++]);
        eb = relative_error_bound * (max - min);
    }

    if (dim == 1) {
        SZ_Compress(data, SZ::Config<float, 1>(eb, {dims[0]}));
    } else {
        SZ_Compress(data, SZ::Config<float, 2>(eb, {dims[0], dims[1]}));
    }

}