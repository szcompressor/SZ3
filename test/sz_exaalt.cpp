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
float SZ_Compress(std::unique_ptr<T[]> const &data, SZ::Config<T, N> conf, int timestep_op, size_t timestep_batch) {
    assert(N == 2);
    conf.quant_state_num = 1024;
    std::cout << "****************** Options ********************" << std::endl;
    std::cout << "dimension = " << N
              << ", error bound = " << conf.eb
              << ", timestep_op = " << timestep_op
              << ", timestep_batch = " << timestep_batch
              << ", quan_state_num = " << conf.quant_state_num
              << ", encoder = " << conf.encoder_op
              << ", lossless = " << conf.lossless_op
              << std::endl;

    std::vector<T> data_(data.get(), data.get() + conf.num);
    std::vector<T> dec_data(conf.num);

    float level_start, level_offset;
    int level_num;
    SZ::get_cluster(data.get(), conf.num, level_start, level_offset, level_num, 0.001);
    //    level_start = -58.291; //trinity-110x
//    level_offset = 2.241; //trinity-110x
//    level_start = 0;
//    level_offset = 1.961;

    printf("start = %.3f , level_offset = %.3f, nlevel=%d\n", level_start, level_offset, level_num);
    double total_compressed_size = 0;
    auto dims = conf.dims;
    auto num = conf.num;
    if (timestep_batch == 0) {
        timestep_batch = dims[0];
    }
    for (size_t ts = 0; ts < dims[0]; ts += timestep_batch) {
        conf.dims[0] = (ts + timestep_batch - 1 > dims[0] ? dims[0] - ts + 1 : timestep_batch);
        conf.num = conf.dims[0] * conf.dims[1];

        std::cout << "****************** Compression ******************" << std::endl;
        SZ::Timer timer(true);
        auto sz = SZ::SZ_Exaalt_Compressor(conf, SZ::LinearQuantizer<float>(conf.eb, conf.quant_state_num / 2),
                                           SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd(), timestep_op);
        sz.set_level(level_start, level_offset, level_num);

        size_t compressed_size = 0;
        std::unique_ptr<SZ::uchar[]> compressed;
        auto ts_data = data.get() + ts * dims[1];
        compressed.reset(sz.compress(ts_data, compressed_size));
        timer.stop("Compression");

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
        timer.stop("Decompression");
        std::copy(ts_dec_data, ts_dec_data + conf.num, dec_data.begin() + ts * dims[1]);

        auto decompressed_file_name = compressed_file_name + ".out";
//    SZ::writefile(decompressed_file_name.c_str(), dec_data.get(), conf.num);
//    std::cout << "Decompressed file = " << decompressed_file_name << std::endl;
        remove(compressed_file_name.c_str());
        total_compressed_size += compressed_size;
    }
    SZ::verify<T>(data_.data(), dec_data.data(), num);

    float ratio = num * sizeof(T) / total_compressed_size;
    std::cout << "Total Compression Size = " << total_compressed_size << std::endl;
    std::cout << "Total Compression Ratio = " << ratio << std::endl;
    return ratio;
}

int main(int argc, char **argv) {


    size_t num = 0;
    auto data = SZ::readfile<float>(argv[1], num);
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
    float eb = 0, relative_eb = 0;
    if (*eb_op == 'a') {
        eb = atof(argv[argp++]);
        relative_eb = eb / (max - min);
    } else {
        relative_eb = atof(argv[argp++]);
        eb = relative_eb * (max - min);
    }
    int timestep_op = 0;
    if (argp < argc) {
        timestep_op = atoi(argv[argp++]);
    }
    size_t timestep_batch = 0;
    if (argp < argc) {
        timestep_batch = atoi(argv[argp++]);
    }
    if (dim == 1) {
        auto conf = SZ::Config<float, 2>(eb, {1, dims[0]});
        conf.relative_eb = relative_eb;
        conf.src_file_name = argv[1];
        SZ_Compress(data, conf, timestep_op, timestep_batch);
    } else {
        auto conf = SZ::Config<float, 2>(eb, {dims[0], dims[1]});
        conf.relative_eb = relative_eb;
        conf.src_file_name = argv[1];
        SZ_Compress(data, conf, timestep_op, timestep_batch);
    }

}