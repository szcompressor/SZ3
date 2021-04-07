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
#include "compressor/SZGeneralCompressor.hpp"
#include "frontend/SZ3TimeBasedFrontend.hpp"
#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "predictor/RegressionPredictor.hpp"
#include "predictor/ComposedPredictor.hpp"
#include "utils/KmeansUtil.h"
#include "utils/QuantOptimizatioin.hpp"
#include <cstdio>

double total_compress_time = 0;
double total_decompress_time = 0;

template<typename T, uint N, class Predictor>
SZ::concepts::CompressorInterface<T> *
make_sz_timebased2(const SZ::Config<T, N> &conf, Predictor predictor, T *data_ts0) {

    return new SZ::SZGeneralCompressor<T, N, SZ::SZ3TimeBasedFrontend<T, N, Predictor, SZ::LinearQuantizer<T>>,
            SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            conf,
            make_sz3_timebased_frontend(
                    conf, predictor, SZ::LinearQuantizer<T>(conf.eb, conf.quant_state_num / 2), data_ts0),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd());
}

template<typename T, uint N>
SZ::concepts::CompressorInterface<T> *
make_sz_timebased(const SZ::Config<T, N> &conf, T *data_ts0) {
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N - 1>>> predictors;

    int use_single_predictor =
            (conf.enable_lorenzo + conf.enable_2ndlorenzo + conf.enable_regression) == 1;
    if (conf.enable_lorenzo) {
        if (use_single_predictor) {
            return make_sz_timebased2(conf, SZ::LorenzoPredictor<T, N - 1, 1>(conf.eb), data_ts0);
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N - 1, 1>>(conf.eb));
        }
    }
    if (conf.enable_2ndlorenzo) {
        if (use_single_predictor) {
            return make_sz_timebased2(conf, SZ::LorenzoPredictor<T, N - 1, 2>(conf.eb), data_ts0);
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N - 1, 2>>(conf.eb));
        }
    }
    if (conf.enable_regression) {
        if (use_single_predictor) {
            return make_sz_timebased2(conf, SZ::RegressionPredictor<T, N - 1>(conf.block_size, conf.eb), data_ts0);
        } else {
            predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N - 1>>(conf.block_size, conf.eb));
        }
    }

    return make_sz_timebased2(conf, SZ::ComposedPredictor<T, N - 1>(predictors), data_ts0);
}

template<typename T, uint N>
float *
VQ(SZ::Config<T, N> conf, size_t ts, T *data, size_t &compressed_size, bool decom,
   int timestep_op, float level_start, float level_offset, int level_num) {
    auto sz = SZ::SZ_Exaalt_Compressor(conf, SZ::LinearQuantizer<float>(conf.eb, conf.quant_state_num / 2),
                                       SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd(), timestep_op);
    sz.set_level(level_start, level_offset, level_num);

    std::unique_ptr<SZ::uchar[]> compressed;
    SZ::Timer timer(true);
    compressed.reset(sz.compress(data, compressed_size));
    total_compress_time += timer.stop("Compression");
    if (!decom) {
        return nullptr;
    }
    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compression Ratio = " << ratio << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;


    std::stringstream ss;
    ss << conf.src_file_name.substr(conf.src_file_name.rfind('/') + 1)
       << ".b" << conf.timestep_batch
       << "." << conf.relative_eb
       << ".md-" << timestep_op
       << ".t" << ts;
    auto compressed_file_name = ss.str();
    SZ::writefile(compressed_file_name.c_str(), compressed.get(), compressed_size);
    std::cout << "Compressed file = " << compressed_file_name << std::endl;

    std::cout << "****************** Decompression ****************" << std::endl;
    compressed = SZ::readfile<SZ::uchar>(compressed_file_name.c_str(), compressed_size);

    timer.start();
    auto ts_dec_data = sz.decompress(compressed.get(), compressed_size);
    total_decompress_time += timer.stop("Decompression");

//    auto decompressed_file_name = compressed_file_name + ".out";
//        SZ::writefile(decompressed_file_name.c_str(), dec_data.get(), conf.num);
//        std::cout << "Decompressed file = " << decompressed_file_name << std::endl;
    remove(compressed_file_name.c_str());

    return ts_dec_data;
}


template<typename T, uint N>
float *
MT(SZ::Config<T, N> conf, size_t ts, T *data, size_t &compressed_size, bool decom,
   int timestep_op, T *data_ts0) {
    T *ts0 = timestep_op == 0 ? nullptr : data_ts0;
    auto sz = make_sz_timebased(conf, ts0);

    SZ::Timer timer(true);
    std::unique_ptr<SZ::uchar[]> compressed;
    compressed.reset(sz->compress(data, compressed_size));
    total_compress_time += timer.stop("Compression");
    if (!decom) {
        return nullptr;
    }

    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compression Ratio = " << ratio << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;


    std::stringstream ss;
    ss << conf.src_file_name.substr(conf.src_file_name.rfind('/') + 1)
       << ".b" << conf.timestep_batch
       << "." << conf.relative_eb
       << ".md-" << timestep_op
       << ".t" << ts;
    auto compressed_file_name = ss.str();
    SZ::writefile(compressed_file_name.c_str(), compressed.get(), compressed_size);
    std::cout << "Compressed file = " << compressed_file_name << std::endl;

    std::cout << "****************** Decompression ****************" << std::endl;
    compressed = SZ::readfile<SZ::uchar>(compressed_file_name.c_str(), compressed_size);

    timer.start();
    auto ts_dec_data = sz->decompress(compressed.get(), compressed_size);
    total_decompress_time += timer.stop("Decompression");

//    auto decompressed_file_name = compressed_file_name + ".out";
//        SZ::writefile(decompressed_file_name.c_str(), dec_data.get(), conf.num);
//        std::cout << "Decompressed file = " << decompressed_file_name << std::endl;

    remove(compressed_file_name.c_str());
    delete sz;
    return ts_dec_data;
}

template<typename T, uint N>
float SZ_Compress(SZ::Config<T, N> conf, int method_op, int method_block) {
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
    auto ts0_name_str = conf.src_file_name + ".ts0";
    auto ts0_name = ts0_name_str.data();
    std::unique_ptr<T[]> data_ts0;
    assert(SZ::file_exist(ts0_name) && "ts0 is required");
    size_t num_ts0;
    data_ts0 = SZ::readfile<T>(ts0_name, num_ts0);


    float level_start, level_offset;
    int level_num;
//    SZ::get_cluster(data.get(), conf.num, level_start, level_offset, level_num, 0.001);
    size_t sample_num = 0.1 * conf.dims[1];
    sample_num = std::min(sample_num, (size_t) 20000);
    sample_num = std::max(sample_num, std::min((size_t) 5000, conf.dims[1]));
    auto data_all = SZ::readfile<T>(conf.src_file_name.data(), 0, conf.num);
    SZ::get_cluster(data_all.get(), conf.dims[1] * conf.timestep_batch, level_start, level_offset, level_num,
                    sample_num);
    //    level_start = -58.291; //trinity-110x
//    level_offset = 2.241; //trinity-110x
//    level_start = 0;
//    level_offset = 1.961;

    printf("start = %.3f , level_offset = %.3f, nlevel=%d\n", level_start, level_offset, level_num);


    auto dims = conf.dims;
    auto total_num = conf.num;
    std::vector<T> dec_data(total_num);
    double total_compressed_size = 0;
    std::vector<std::string> compressor_names = {"VQ", "VQT", "MT"};

    for (size_t ts = 0; ts < dims[0]; ts += conf.timestep_batch) {
        conf.dims[0] = (ts + conf.timestep_batch > dims[0] ? dims[0] - ts : conf.timestep_batch);
        conf.num = conf.dims[0] * conf.dims[1];

//        auto data = SZ::readfile<T>(conf.src_file_name.data(), ts * conf.dims[1], conf.num);
        T *data = &data_all[ts * conf.dims[1]];

        T max = *std::max_element(data, data + conf.num);
        T min = *std::min_element(data, data + conf.num);
        if (conf.eb_mode == 0) {
            conf.relative_eb = conf.eb / (max - min);
        } else if (conf.eb_mode == 1) {
            conf.eb = conf.relative_eb * (max - min);
        }

        std::cout << "****************** Compression From " << ts << " to " << ts + conf.dims[0] - 1
                  << " ******************" << std::endl;
        if (level_num > 0 && method_block > 0 && ts / conf.timestep_batch % method_block == 0) {
            std::vector<size_t> compressed_size = {0, 0, 0};
            std::vector<T> data1(&data_all[ts * conf.dims[1]], &data_all[ts * conf.dims[1]] + conf.num);
            VQ(conf, ts, data1.data(), compressed_size[0], false, 0, level_start, level_offset, level_num);

            data1 = std::vector(&data_all[ts * conf.dims[1]], &data_all[ts * conf.dims[1]] + conf.num);
            VQ(conf, ts, data1.data(), compressed_size[1], false, 1, level_start, level_offset, level_num);

            data1 = std::vector(&data_all[ts * conf.dims[1]], &data_all[ts * conf.dims[1]] + conf.num);
            MT(conf, ts, data1.data(), compressed_size[2], false, 1, data_ts0.get());
            method_op = std::distance(compressed_size.begin(),
                                      std::min_element(compressed_size.begin(), compressed_size.end()));
            printf("Select %s as Compressor, timestep=%lu, %d %lu %lu %lu\n", compressor_names[method_op].data(),
                   ts, method_op, compressed_size[0], compressed_size[1], compressed_size[2]);
            std::cout << "****************** END Selection ****************" << std::endl;
        }
        T *ts_dec_data;
        size_t compressed_size;
        if (method_op == 0) {
            ts_dec_data = VQ(conf, ts, data, compressed_size, true, 0, level_start, level_offset, level_num);
        } else if (method_op == 1) {
            ts_dec_data = VQ(conf, ts, data, compressed_size, true, 1, level_start, level_offset, level_num);
        } else {
            ts_dec_data = MT(conf, ts, data, compressed_size, true, 1, data_ts0.get());
        }
        total_compressed_size += compressed_size;
        memcpy(&dec_data[ts * conf.dims[1]], ts_dec_data, conf.num * sizeof(T));
    }

    std::cout << "****************** Final ****************" << std::endl;
    float ratio = total_num * sizeof(T) / total_compressed_size;
    auto data = SZ::readfile<T>(conf.src_file_name.data(), 0, total_num);

    std::stringstream ss;
    ss << conf.src_file_name.substr(conf.src_file_name.rfind('/') + 1)
       << ".b" << conf.timestep_batch
       << "." << conf.relative_eb << ".md-" << conf.timestep_op << ".out";
    std::cout << "Decompressed file = " << ss.str() << std::endl;
    SZ::writefile(ss.str().data(), dec_data.data(), total_num);

    double max_diff, psnr, nrmse;
    SZ::verify<T>(data.get(), dec_data.data(), total_num, max_diff, psnr, nrmse);

    printf("method=md, file=%s, block=%lu, compression_ratio=%.3f, reb=%.1e, eb=%.6f, psnr=%.3f, nsmse=%e, compress_time=%.3f, decompress_time=%.3f, timestep_op=%d\n",
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
    int method_op = 2; //default compressor: MT
    if (argp < argc) {
        method_op = atoi(argv[argp++]);
    }
    int method_batch = 50;
    if (argp < argc) {
        method_batch = atoi(argv[argp++]);
    }
    conf.block_size = 128;
    conf.stride = 128;
    SZ_Compress(conf, method_op, method_batch);

}