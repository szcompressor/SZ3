#include <vector>
#include <sstream>
#include "mdz.h"
#include "utils/FileUtil.h"
#include "utils/Timer.hpp"
#include "utils/Verification.hpp"
#include "def.hpp"
#include "quantizer/IntegerQuantizer.hpp"
#include "lossless/Lossless_zstd.hpp"
#include "encoder/HuffmanEncoder.hpp"
#include "frontend/SZ3Frontend.hpp"
#include "utils/KmeansUtil.h"
#include "utils/QuantOptimizatioin.hpp"
#include <cstdio>

double total_compress_time = 0;
double total_decompress_time = 0;

const char *compressor_names[] = {"VQ", "VQT", "MT", "LR", "TS"};

template<typename T, uint N>
float *
VQ(SZ::Config<T, N> conf, size_t ts, T *data, size_t &compressed_size, bool decom,
   int method, float level_start, float level_offset, int level_num) {
    if (level_num == 0) {
        printf("VQ/VQT not availble on current dataset, please use ADP or MT\n");
        exit(0);
    }
    SZ::Timer timer(true);

    auto sz = SZ::SZ_Exaalt_Compressor(conf, SZ::LinearQuantizer<float>(conf.eb, conf.quant_state_num / 2),
                                       SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd(), method);
    sz.set_level(level_start, level_offset, level_num);

    std::unique_ptr<SZ::uchar[]> compressed;
    compressed.reset(sz.compress(data, compressed_size));
    total_compress_time += timer.stop("Compression");
    if (!decom) {
        return nullptr;
    }
    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compression Ratio = " << ratio << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;


//    std::stringstream ss;
//    ss << conf.src_file_name.substr(conf.src_file_name.rfind('/') + 1)
//       << ".b" << conf.timestep_batch
//       << "." << conf.relative_eb
//       << ".md-" << method
//       << ".t" << ts;
//    auto compressed_file_name = ss.str();
//    SZ::writefile(compressed_file_name.c_str(), compressed.get(), compressed_size);
//    std::cout << "Compressed file = " << compressed_file_name << std::endl;
//
//    std::cout << "****************** Decompression ****************" << std::endl;
//    compressed = SZ::readfile<SZ::uchar>(compressed_file_name.c_str(), compressed_size);

    timer.start();
    auto ts_dec_data = sz.decompress(compressed.get(), compressed_size);
    total_decompress_time += timer.stop("Decompression");

//    auto decompressed_file_name = compressed_file_name + ".out";
//        SZ::writefile(decompressed_file_name.c_str(), dec_data.get(), conf.num);
//        std::cout << "Decompressed file = " << decompressed_file_name << std::endl;
//    remove(compressed_file_name.c_str());

    return ts_dec_data;
}


template<typename T, uint N>
float *
MT(SZ::Config<T, N> conf, size_t ts, T *data, size_t &compressed_size, bool decom, T *ts0) {
    SZ::Timer timer(true);

    auto sz = make_sz_timebased(conf, ts0);

    std::unique_ptr<SZ::uchar[]> compressed;
    compressed.reset(sz->compress(data, compressed_size));
    total_compress_time += timer.stop("Compression");
    if (!decom) {
        return nullptr;
    }

    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compression Ratio = " << ratio << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;


//    std::stringstream ss;
//    ss << conf.src_file_name.substr(conf.src_file_name.rfind('/') + 1)
//       << ".b" << conf.timestep_batch
//       << "." << conf.relative_eb
//       << ".md-2"
//       << ".t" << ts;
//    auto compressed_file_name = ss.str();
//    SZ::writefile(compressed_file_name.c_str(), compressed.get(), compressed_size);
//    std::cout << "Compressed file = " << compressed_file_name << std::endl;
//
//    std::cout << "****************** Decompression ****************" << std::endl;
//    compressed = SZ::readfile<SZ::uchar>(compressed_file_name.c_str(), compressed_size);

    timer.start();
    auto ts_dec_data = sz->decompress(compressed.get(), compressed_size);
    total_decompress_time += timer.stop("Decompression");

//    auto decompressed_file_name = compressed_file_name + ".out";
//        SZ::writefile(decompressed_file_name.c_str(), dec_data.get(), conf.num);
//        std::cout << "Decompressed file = " << decompressed_file_name << std::endl;

//    remove(compressed_file_name.c_str());
    delete sz;
    return ts_dec_data;
}

template<typename T, uint N>
float *
SZ2(SZ::Config<T, N> conf, size_t ts, T *data, size_t &compressed_size, bool decom) {
    SZ::Timer timer(true);
//    conf.block_size = 16;
//    conf.enable_regression=false;
//    conf.stride = conf.block_size;
    auto sz = make_sz(conf);

    std::unique_ptr<SZ::uchar[]> compressed;
    compressed.reset(sz->compress(data, compressed_size));
    total_compress_time += timer.stop("Compression");
    if (!decom) {
        return nullptr;
    }

    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compression Ratio = " << ratio << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;


//    std::stringstream ss;
//    ss << conf.src_file_name.substr(conf.src_file_name.rfind('/') + 1)
//       << ".b" << conf.timestep_batch
//       << "." << conf.relative_eb
//       << ".md-3"
//       << ".t" << ts;
//    auto compressed_file_name = ss.str();
//    SZ::writefile(compressed_file_name.c_str(), compressed.get(), compressed_size);
//    std::cout << "Compressed file = " << compressed_file_name << std::endl;
//
//    std::cout << "****************** Decompression ****************" << std::endl;
//    compressed = SZ::readfile<SZ::uchar>(compressed_file_name.c_str(), compressed_size);

    timer.start();
    auto ts_dec_data = sz->decompress(compressed.get(), compressed_size);
    total_decompress_time += timer.stop("Decompression");

//    auto decompressed_file_name = compressed_file_name + ".out";
//        SZ::writefile(decompressed_file_name.c_str(), dec_data.get(), conf.num);
//        std::cout << "Decompressed file = " << decompressed_file_name << std::endl;

//    remove(compressed_file_name.c_str());
    delete sz;
    return ts_dec_data;
}

int ts_last_select = -1;
int method_batch = 0;

template<typename T, uint N>
void select(SZ::Config<T, N> conf, int &method, size_t ts, T *data_all,
            float level_start, float level_offset, int level_num, T *data_ts0) {
    if (method_batch <= 0) {
        return;
    }
//        && (ts_last_select == -1 || t - ts_last_select >= conf.timestep_batch * 10)) {
    std::cout << "****************** BEGIN Selection ****************" << std::endl;
//        ts_last_select = ts;
    std::vector<size_t> compressed_size(10, std::numeric_limits<size_t>::max());
    std::vector<T> data1;
    if (conf.timestep_batch > 10) {
        conf.dims[0] = 10;
        conf.num = conf.dims[0] * conf.dims[1];
    }
    size_t t = ts;
    if (ts == 0) {
        t = conf.dims[0] / 2;
        conf.dims[0] /= 2;
        conf.num = conf.dims[0] * conf.dims[1];
    }
    std::cout << conf.dims[0] << " " << conf.dims[1] << " " << t << std::endl;

    if (level_num > 0) {
        data1 = std::vector(&data_all[t * conf.dims[1]], &data_all[t * conf.dims[1]] + conf.num);
        VQ(conf, t, data1.data(), compressed_size[0], false, 0, level_start, level_offset, level_num);

        data1 = std::vector(&data_all[t * conf.dims[1]], &data_all[t * conf.dims[1]] + conf.num);
        VQ(conf, t, data1.data(), compressed_size[1], false, 1, level_start, level_offset, level_num);
    } else {
        data1 = std::vector(&data_all[t * conf.dims[1]], &data_all[t * conf.dims[1]] + conf.num);
        SZ2(conf, t, data1.data(), compressed_size[3], false);
    }

    data1 = std::vector(&data_all[t * conf.dims[1]], &data_all[t * conf.dims[1]] + conf.num);
    MT(conf, t, data1.data(), compressed_size[2], false, data_ts0);

//    data1 = std::vector(&data_all[t * conf.dims[1]], &data_all[t * conf.dims[1]] + conf.num);
//    MT(conf, t, data1.data(), compressed_size[4], false, (T *) nullptr);

    method = std::distance(compressed_size.begin(),
                           std::min_element(compressed_size.begin(), compressed_size.end()));
    printf("Select %s as Compressor, timestep=%lu, method=%d, %lu %lu %lu %lu %lu\n",
           compressor_names[method],
           ts, method, compressed_size[0], compressed_size[1], compressed_size[2], compressed_size[3],
           compressed_size[4]);
    std::cout << "****************** END Selection ****************" << std::endl;
}

template<typename T, uint N>
float SZ_Compress(SZ::Config<T, N> conf, int method) {
    assert(N == 2);
    if (conf.timestep_batch == 0) {
        conf.timestep_batch = conf.dims[0];
    }
    std::cout << "****************** Options ********************" << std::endl;
    std::cout << "dimension = " << N
              << ", error bound = " << conf.eb
              << ", method = " << method
              << ", method_update_batch = " << method_batch
              << ", timestep_batch = " << conf.timestep_batch
              << ", quan_state_num = " << conf.quant_state_num
              << ", encoder = " << conf.encoder_op
              << ", lossless = " << conf.lossless_op
              << std::endl;
//    auto ts0_name_str = conf.src_file_name + ".ts0";
//    auto ts0_name = ts0_name_str.data();
//    std::unique_ptr<T[]> data_ts0;
//    assert(SZ::file_exist(ts0_name) && "ts0 is required");
//    size_t num_ts0;
//    data_ts0 = SZ::readfile<T>(ts0_name, num_ts0);
    auto data_all = SZ::readfile<T>(conf.src_file_name.data(), 0, conf.num);
    auto data_ts0 = std::vector<T>(data_all.get(), data_all.get() + conf.dims[1]);

    float level_start, level_offset;
    int level_num;
    size_t sample_num = 0.1 * conf.dims[1];
    sample_num = std::min(sample_num, (size_t) 20000);
    sample_num = std::max(sample_num, std::min((size_t) 5000, conf.dims[1]));
    SZ::get_cluster(data_all.get(), conf.dims[1], level_start, level_offset, level_num,
                    sample_num);
    if (level_num > conf.dims[1] * 0.25) {
        level_num = 0;
    }
    if (level_num != 0) {
        printf("start = %.3f , level_offset = %.3f, nlevel=%d\n", level_start, level_offset, level_num);
    }

    auto dims = conf.dims;
    auto total_num = conf.num;
    std::vector<T> dec_data(total_num);
    double total_compressed_size = 0;
    double compressed_size_pre = total_num * sizeof(T);
    int current_method = method;

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
//        std::cout<<method_batch<<" "<<ts<<" "<<conf.timestep_batch<<" "<<method_batch<<std::endl;
        if (method_batch > 0 && ts / conf.timestep_batch % method_batch == 0) {
            select(conf, current_method, ts, data_all.get(), level_start, level_offset, level_num, data_ts0.data());
        }
        printf("Compressor = %s\n", compressor_names[current_method]);

        T *ts_dec_data;
        size_t compressed_size;
        if (current_method == 0) {
            ts_dec_data = VQ(conf, ts, data, compressed_size, true, current_method, level_start, level_offset,
                             level_num);
        } else if (current_method == 1) {
            ts_dec_data = VQ(conf, ts, data, compressed_size, true, current_method, level_start, level_offset,
                             level_num);
        } else if (current_method == 2) {
            ts_dec_data = MT(conf, ts, data, compressed_size, true, data_ts0.data());
        } else if (current_method == 4) {
            ts_dec_data = MT(conf, ts, data, compressed_size, true, (T *) nullptr);
        } else {
            ts_dec_data = SZ2(conf, ts, data, compressed_size, true);
        }
        total_compressed_size += compressed_size;
//        if (compressed_size > 4.0 * compressed_size_pre) {
//            select(conf, current_method, ts, data_all.get(), level_start, level_offset, level_num, data_ts0.get());
//        }
        compressed_size_pre = compressed_size;
        memcpy(&dec_data[ts * conf.dims[1]], ts_dec_data, conf.num * sizeof(T));
    }

    std::cout << "****************** Final ****************" << std::endl;

    float ratio = total_num * sizeof(T) / total_compressed_size;
    auto data = SZ::readfile<T>(conf.src_file_name.data(), 0, total_num);

    std::stringstream ss;
    ss << conf.src_file_name.substr(conf.src_file_name.rfind('/') + 1)
       << ".b" << conf.timestep_batch
       << "." << conf.relative_eb << ".md-" << method << ".out";
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
           method);

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
    int method = 9;
    method_batch = 50; //method_batch >0 indicates ADP
    if (argp < argc) {
        int tmp = atoi(argv[argp++]);
        if (tmp <= 0) {
            method = -tmp;
            method_batch = 0;
        } else {
            method_batch = tmp;
        }
    }

    conf.block_size = 128;
    conf.stride = 128;

    conf.quant_state_num = 1024;
    if (argp < argc) {
        conf.quant_state_num = atoi(argv[argp++]);
    }
    SZ_Compress(conf, method);

}