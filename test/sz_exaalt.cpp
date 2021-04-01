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
#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>
#include <libgen.h>

/**
 *
 * @tparam T
 * @tparam N
 * @param data
 * @param conf
 * @param level_start
 * @param level_offset
 * @param level_num
 * @param timestep_op 0:no timestep, 1:timestep, 2: timestep+levels
 * @return
 */
template<typename T, uint N>
float SZ_Compress(SZ::Config<T, N> conf) {
    assert(N == 2);
    conf.quant_state_num = 1024;
    if (conf.timestep_batch == 0) {
        conf.timestep_batch = conf.dims[0];
    }
    std::cout << "****************** Options ********************" << std::endl;
    std::cout << "dimension = " << N
              << ", error bound = " << conf.eb
              << ", timestep_op = " << conf.timestep_op
              << ", timestep_batch = " << conf.timestep_batch
              << ", quan_state_num = " << conf.quant_state_num
              << ", encoder = " << conf.encoder_op
              << ", lossless = " << conf.lossless_op
              << std::endl;


    float level_start, level_offset;
    int level_num;
//    SZ::get_cluster(data.get(), conf.num, level_start, level_offset, level_num, 0.001);
    auto sample_rate = (conf.timestep_batch == conf.dims[0]) ? 0.001 : 1;
    auto data_all = SZ::readfile<T>(conf.src_file_name.data(), 0, conf.num);
    SZ::get_cluster(data_all.get(), conf.dims[1] * conf.timestep_batch, level_start, level_offset, level_num,
                    sample_rate);
    //    level_start = -58.291; //trinity-110x
//    level_offset = 2.241; //trinity-110x
//    level_start = 0;
//    level_offset = 1.961;

    printf("start = %.3f , level_offset = %.3f, nlevel=%d\n", level_start, level_offset, level_num);

    double total_compressed_size = 0;
    double total_compress_time = 0;
    double total_decompress_time = 0;
    auto dims = conf.dims;
    auto total_num = conf.num;
    std::vector<T> dec_data(total_num);

    for (size_t ts = 0; ts < dims[0]; ts += conf.timestep_batch) {
        conf.dims[0] = (ts + conf.timestep_batch > dims[0] ? dims[0] - ts : conf.timestep_batch);
        conf.num = conf.dims[0] * conf.dims[1];

        auto data = SZ::readfile<T>(conf.src_file_name.data(), ts, conf.num);
        T max = *std::max_element(data.get(), data.get() + conf.num);
        T min = *std::min_element(data.get(), data.get() + conf.num);
        if (conf.eb_mode == 0) {
            conf.relative_eb = conf.eb / (max - min);
        } else if (conf.eb_mode == 1) {
            conf.eb = conf.relative_eb * (max - min);
        }

        std::cout << "****************** Compression From " << ts << " to " << ts + conf.dims[0] - 1
                  << " ******************" << std::endl;
        SZ::Timer timer(true);
        auto sz = SZ::SZ_Exaalt_Compressor(conf, SZ::LinearQuantizer<float>(conf.eb, conf.quant_state_num / 2),
                                           SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd(), conf.timestep_op);
        sz.set_level(level_start, level_offset, level_num);

        size_t compressed_size = 0;
        std::unique_ptr<SZ::uchar[]> compressed;
        compressed.reset(sz.compress(data.get(), compressed_size));
        total_compress_time += timer.stop("Compression");

        auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
        std::cout << "Compression Ratio = " << ratio << std::endl;
        std::cout << "Compressed size = " << compressed_size << std::endl;

        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> dis(0, 10000);
        std::stringstream ss;
        ss << conf.src_file_name.substr(conf.src_file_name.rfind('/') + 1)
           << "." << conf.relative_eb << "." << dis(gen) << ".sz3";
        auto compressed_file_name = ss.str();
        SZ::writefile(compressed_file_name.c_str(), compressed.get(), compressed_size);
        std::cout << "Compressed file = " << compressed_file_name << std::endl;

        std::cout << "****************** Decompression ****************" << std::endl;
        compressed = SZ::readfile<SZ::uchar>(compressed_file_name.c_str(), compressed_size);

        timer.start();
        auto ts_dec_data = sz.decompress(compressed.get(), compressed_size);
        total_decompress_time += timer.stop("Decompression");

        auto decompressed_file_name = compressed_file_name + ".out";
//    SZ::writefile(decompressed_file_name.c_str(), dec_data.get(), conf.num);
//    std::cout << "Decompressed file = " << decompressed_file_name << std::endl;
        remove(compressed_file_name.c_str());
        memcpy(&dec_data[ts * conf.dims[1]], ts_dec_data, conf.num * sizeof(T));
        total_compressed_size += compressed_size;
    }

    std::cout << "****************** Final ****************" << std::endl;
    float ratio = total_num * sizeof(T) / total_compressed_size;
    auto data = SZ::readfile<T>(conf.src_file_name.data(), 0, total_num);

    std::stringstream ss;
    ss << dirname(strdup(conf.src_file_name.c_str()))
       << "/exaalt/";
    struct stat st = {0};
    if (stat(ss.str().c_str(), &st) == -1) {
        mkdir(ss.str().c_str(), 0700);
    }
    ss << basename(strdup(conf.src_file_name.c_str()))
       << ".b" << conf.timestep_batch
       << "." << conf.relative_eb << ".out";
    std::cout << "Decompressed file = " << ss.str() << std::endl;
    SZ::writefile(ss.str().data(), dec_data.data(), total_num);

    double max_diff, psnr, nrmse;
    SZ::verify<T>(data.get(), dec_data.data(), total_num, max_diff, psnr, nrmse);

    printf("method=exaalt, file=%s, block=%lu, compression_ratio=%.3f, reb=%.1e, eb=%.6f, psnr=%.3f, nsmse=%e, compress_time=%.3f, decompress_time=%.3f, timestep_op=%d\n",
           conf.src_file_name.data(), conf.timestep_batch,
           ratio,
           conf.relative_eb,
           max_diff, psnr, nrmse,
           total_compress_time, total_decompress_time,
           conf.timestep_op);

    return ratio;
}

int main(int argc, char **argv) {


    int dim = atoi(argv[2] + 1);
    assert(1 <= dim && dim <= 2);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }

    SZ::Config<float, 2> conf({1, dims[0]});
    if (dim == 2) {
        conf = SZ::Config<float, 2>({dims[0], dims[1]});
    }
    conf.src_file_name = argv[1];

    char *eb_op = argv[argp++] + 1;
    if (*eb_op == 'a') {
        conf.eb_mode = 0;
        conf.eb = atof(argv[argp++]);
    } else {
        conf.eb_mode = 1;
        conf.relative_eb = atof(argv[argp++]);
    }

    if (argp < argc) {
        conf.timestep_batch = atoi(argv[argp++]);
    }
    if (argp < argc) {
        conf.timestep_op = atoi(argv[argp++]);
    }
    SZ_Compress(conf);

}