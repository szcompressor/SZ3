#include <compressor/SZInterpolationCompressor.hpp>
#include <compressor/SZBlockInterpolationCompressor.hpp>
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
interp_compress_decompress_nowritefile(char *path, float *data, size_t num, double eb, int interp_level,
                                       int interp_op, int direction_op,
                                       int block_size, int interp_block_size, SZ_Option sz_op, Dims ... args) {
    std::string compressed_file_name(path);
    META::meta_compress_info compressInfo;


    std::cout << "****************** compression ****************" << std::endl;
    std::cout << "Interp Level       = " << interp_level << std::endl
              << "Interp Op          = " << interp_op << std::endl
              << "Direction          = " << direction_op << std::endl
              << "SZ block size      = " << block_size << std::endl
              << "Interp block size  = " << interp_block_size << std::endl
              << "SZ_mode            = " << sz_op << std::endl;

    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);
    size_t compressed_size = 0;
    std::unique_ptr<unsigned char[]> compressed;

    auto dims = std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...};
    {
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
        compressed.reset(sz.compress(data, compressed_size, true));
    }

    clock_gettime(CLOCK_REALTIME, &end);
    double compress_time =
            (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
    auto compression_ratio = num * sizeof(float) * 1.0 / compressed_size;
    compressInfo.compress_time = compress_time;
    std::cout << "LEVEL2 Interp compression time = " << compress_time << "s" << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;
    std::cout << "Compression ratio = " << compression_ratio << std::endl;

    std::cout << "***************** Decompression ****************" << std::endl;

    clock_gettime(CLOCK_REALTIME, &start);

    std::unique_ptr<float[]> dec_data;
    {
        auto sz_dec = SZ::SZInterpolationCompressor<float, N, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer<float>(eb),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd(),
                dims,
                interp_block_size,
                interp_op,
                direction_op,
                interp_level
        );
        dec_data.reset(sz_dec.decompress(compressed.get(), compressed_size));
    }
    compressed.reset();

    clock_gettime(CLOCK_REALTIME, &end);
    double decompress_time =
            (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
    compressInfo.decompress_time = decompress_time;
    std::cout << "LEVEL2 Interp decompression time = " << decompress_time << "s" << std::endl;

    size_t num1 = 0;
    auto ori_data = SZ::readfile<float>(path, num1);
    assert(num1 == num);
    double psnr, nrmse;
    SZ::verify<float>(ori_data.get(), dec_data.get(), num, psnr, nrmse);

    compression_ratio = num * sizeof(float) * 1.0 / compressed_size;
    printf("PSNR = %f, NRMSE = %.10G, Compression Ratio = %.2f\n", psnr, nrmse, compression_ratio);
    std::cout << "LEVEL2 " << compression_ratio
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

template<uint N, class ... Dims>
META::meta_compress_info
interp_compress_decompress(char *path, float *data, size_t num, double eb, int interp_level, int interp_op, int direction_op,
                           int block_size, int interp_block_size, SZ_Option sz_op, Dims ... args) {
    std::string compressed_file_name(path);
    META::meta_compress_info compressInfo;

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
        compressed.reset(sz.compress(data, compressed_size));


        clock_gettime(CLOCK_REALTIME, &end);
        double compress_time =
                (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
        auto compression_ratio = num * sizeof(float) * 1.0 / compressed_size;
        compressInfo.compress_time = compress_time;
        std::cout << "LEVEL2 Interp compression time = " << compress_time << "s" << std::endl;
        std::cout << "Compressed size = " << compressed_size << std::endl;
        std::cout << "Compression ratio = " << compression_ratio << std::endl;

        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(0, 100000000);
        std::stringstream ss;
        ss << compressed_file_name.substr(compressed_file_name.rfind('/') + 1)
           << "_" << eb << "_" << dis(gen) << ".sz";
        compressed_file_name = ss.str();
        std::cout << "CompressFileName = " << compressed_file_name << std::endl;
        fflush(stdout);
        SZ::writefile(compressed_file_name.c_str(), compressed.get(), compressed_size);
    }
    {
        std::cout << "***************** Decompression ****************" << std::endl;
        size_t compressed_size = 0;
        auto compressed = SZ::readfile<SZ::uchar>(compressed_file_name.c_str(), compressed_size);
        remove(compressed_file_name.c_str());

        struct timespec start, end;
        clock_gettime(CLOCK_REALTIME, &start);
        std::unique_ptr<float[]> dec_data;

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
        dec_data.reset(sz.decompress(compressed.get(), compressed_size));


        clock_gettime(CLOCK_REALTIME, &end);
        double decompress_time =
                (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
        compressInfo.decompress_time = decompress_time;
        std::cout << "LEVEL2 Interp decompression time = " << decompress_time << "s" << std::endl;

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
        SZ::verify<float>(ori_data.get(), dec_data.get(), num, psnr, nrmse);

//        std::vector<float> error(num);
//        for (size_t i = 0; i < num; i++) {
//            error[i] = ori_data[i] - dec_data[i];
//        }
//        std::string error_file(path);
//        error_file += ".error";
//        SZ::writefile(error_file.c_str(), error.data(), num);
        auto compression_ratio = num * sizeof(float) * 1.0 / compressed_size;
        printf("PSNR = %f, NRMSE = %.10G, Compression Ratio = %.2f\n", psnr, nrmse, compression_ratio);
        std::cout << "LEVEL2 " << compression_ratio
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
double
interp_compress_test_block(float *data, size_t num, double eb, int interp_level, int interp_op, int direction_op, int block_size,
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

    SZ::Config<float, N> conf(eb, std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...});
    conf.block_size = block_size;
    conf.stride = conf.block_size;
    auto sz = SZ::make_sz_fast_block_interpolation_compressor(
            conf,
            SZ::SimplePredictor<float, N>(eb),
            SZ::LinearQuantizer<float>(eb),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd(),
            interp_op,
            direction_op,
            interp_level
    );

//    auto dims = std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...};
//    auto sz = SZ::SZBlockInterpalationCompressor<float, N, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
//            conf,
//            SZ::LinearQuantizer<float>(eb),
//            SZ::HuffmanEncoder<int>(),
//            SZ::Lossless_zstd(),
//            interp_op,
//            direction_op,
//            interp_level
//    );

    sz.compress(data1.data(), compressed_size);
    clock_gettime(CLOCK_REALTIME, &end);
    double compression_time = (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
    auto compression_ratio = num * sizeof(float) * 1.0 / compressed_size;
    std::cout << "Compressed size = " << compressed_size << std::endl;
    std::cout << "Compression ratio = " << compression_ratio << std::endl;
    std::cout << "LEVEL1 Interp Compression time = " << compression_time
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
double interp_compress_test(float *data, size_t num, double eb, int interp_level, int interp_op, int direction_op, int block_size,
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
    double compression_time = (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
    auto compression_ratio = num * sizeof(float) * 1.0 / compressed_size;
    std::cout << "Compressed size = " << compressed_size << std::endl;
    std::cout << "Compression ratio = " << compression_ratio << std::endl;
    std::cout << "LEVEL1 Interp Compression time = " << compression_time
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
        ratio = interp_compress_test<N>(data.get(), num, eb, interp_level, i, direction_op, block_size, interp_block_size, sz_op,
                                        args...);
        if (ratio > best_ratio) {
            best_ratio = ratio;
            interp_op = i;
        }
    }
    std::cout << "LEVEL1 best interp_op = " << interp_op << " , best ratio = " << best_ratio << std::endl;

    ratio = interp_compress_test<N>(data.get(), num, eb, interp_level, interp_op, 5, block_size, interp_block_size, sz_op,
                                    args...);
    if (ratio > best_ratio * 1.0) {
        best_ratio = ratio;
        direction_op = 5;
    }
    std::cout << "LEVEL1 best direction_op = " << direction_op << " , best ratio = " << best_ratio << std::endl;

    return interp_compress_decompress<N>(path, data.get(), num, eb, interp_level, interp_op, direction_op, block_size,
                                         interp_block_size, sz_op,
                                         args...);
}

template<typename T>
META::meta_compress_info meta_compress_3d(T *data, size_t num_elements, int r1, int r2, int r3, float precision,
                                          META::meta_params params) {
    size_t result_size = 0;
    META::meta_compress_info compressInfo;

    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);

    unsigned char *result = meta_compress_3d<T>(data, r1, r2, r3, precision, result_size, params, compressInfo);
    unsigned char *result_after_lossless = NULL;
    size_t lossless_outsize = META::meta_lossless_compress(ZSTD_COMPRESSOR, 3, result, result_size, &result_after_lossless);
    free(result);
    free(result_after_lossless);

    clock_gettime(CLOCK_REALTIME, &end);
    compressInfo.compress_time = (float) (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / (double) 1000000000;

    float ratio = (num_elements * sizeof(T)) * 1.0 / lossless_outsize;
    compressInfo.ori_bytes = num_elements * sizeof(T);
    compressInfo.compress_bytes = lossless_outsize;
    compressInfo.ratio = compressInfo.ori_bytes * 1.0 / compressInfo.compress_bytes;

    return compressInfo;
}

template<uint N>
float cal_sampling_ratio(size_t block, size_t n, size_t dmin, std::array<size_t, N> dims) {
    size_t sample_n = 1;
    for (auto dim:dims) {
        sample_n *= dim / dmin * 2 * block;
    }
    return sample_n * 1.0 / n;
}

template<uint N, class ... Dims>
void interp_meta_tuning(char *path, double reb, Dims ... args) {
    size_t num = 0;
    auto data = SZ::readfile<float>(path, num);
    std::cout << "Read " << num << " elements\n";
    auto dims = std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...};
    float max = data[0];
    float min = data[0];
    for (int i = 1; i < num; i++) {
        if (max < data[i]) max = data[i];
        if (min > data[i]) min = data[i];
    }
    float eb = reb * (max - min);

    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);

    bool meta = true;
    size_t dmin = *std::min_element(dims.begin(), dims.end());
    int sampling_block = dmin;
    while (cal_sampling_ratio<3>(sampling_block, num, dmin, dims) > 0.035) {
        sampling_block--;
    }
    if (sampling_block * 2 > dmin) {
        sampling_block = dmin / 2;
    }
    size_t b1 = dims[0] / dmin;
    size_t b2 = dims[1] / dmin;
    size_t b3 = dims[2] / dmin;
    size_t r1 = b1 * 2 * sampling_block;
    size_t r2 = b2 * 2 * sampling_block;
    size_t r3 = b3 * 2 * sampling_block;
    size_t di, dj, dk;
    size_t sampling_num = r1 * r2 * r3;
    std::vector<float> sampling_data(sampling_num, 0);
//        std::vector<int> debug(sampling_num, 0);
    struct timespec start1, end1;
    clock_gettime(CLOCK_REALTIME, &start1);
    for (size_t bi = 0; bi < b1; bi++) {
        for (size_t bj = 0; bj < b2; bj++) {
            for (size_t bk = 0; bk < b3; bk++) {
                for (size_t i = 0; i < 2 * sampling_block; i++) {
                    for (size_t j = 0; j < 2 * sampling_block; j++) {
                        for (size_t k = 0; k < 2 * sampling_block; k++) {
                            di = i < sampling_block ? i + sampling_block : dmin - 3 * sampling_block + i;
                            dj = j < sampling_block ? j + sampling_block : dmin - 3 * sampling_block + j;
                            dk = k < sampling_block ? k + sampling_block : dmin - 3 * sampling_block + k;
                            auto d = data[(bi * dmin + di) * dims[1] * dims[2] + (bj * dmin + dj) * dims[2] + bk * dmin + dk];
                            sampling_data[(bi * 2 * sampling_block + i) * r2 * r3
                                          + (bj * 2 * sampling_block + j) * r3
                                          + bk * 2 * sampling_block + k] = d;
//                                debug[(bi * 2 * sampling_block + i) * r2 * r3
//                                      + (bj * 2 * sampling_block + j) * r3
//                                      + bk * 2 * sampling_block + k]++;
                        }
                    }
                }
            }
        }
    }
    clock_gettime(CLOCK_REALTIME, &end1);
    auto sampling_time = (double) (end1.tv_sec - start1.tv_sec) + (double) (end1.tv_nsec - start1.tv_nsec) / (double) 1000000000;
    printf("LEVEL1 sampling block = %d percent = %.3f%% Time = %.3f \n", sampling_block, sampling_num * 100.0 / num,
           sampling_time);
//        for (size_t i = 0; i < sampling_num; i++) {
//            if (debug[i] != 1) {
//                printf("%lu %d\n", i, debug[i]);
//                exit(0);
//            }
//        }

    META::meta_params meta_params(false, 5, 3, 0, true, true, false, eb);
    meta_params.capacity = 65536 * 2;
    meta_params.lossless = true;
    int interp_level = -1, interp_op, direction_op = 0, block_size = sampling_block, interp_block_size = sampling_block;
    SZ_Option sz_op = SZ_LR;

    auto result_meta = meta_compress_3d(sampling_data.data(), sampling_num, r1, r2, r3, eb, meta_params);
    printf("LEVEL1 meta ratio = %.2f, lorenzo:%d, lorenzo2:%d, pred_dim:%d compress_time:%.3f\n",
           result_meta.ratio, meta_params.use_lorenzo, meta_params.use_lorenzo_2layer, meta_params.prediction_dim,
           result_meta.compress_time);

    double best_meta_ratio = result_meta.ratio, best_interp_ratio = 0, ratio;
    for (int i = 0; i < 2; i++) {
        ratio = interp_compress_test_block<N>(sampling_data.data(), sampling_num, eb, interp_level, i, direction_op,
                                              block_size, interp_block_size, sz_op, r1, r2, r3);
        if (ratio > best_interp_ratio) {
            best_interp_ratio = ratio;
            interp_op = i;
        }
    }
    std::cout << "LEVEL1 interp best interp_op = " << interp_op << " , best ratio = " << best_interp_ratio << std::endl;

    ratio = interp_compress_test_block<N>(sampling_data.data(), sampling_num, eb, interp_level, interp_op, 5, block_size,
                                          interp_block_size, sz_op, r1, r2, r3);
    if (ratio > best_interp_ratio * 1.02) {
        best_interp_ratio = ratio;
        direction_op = 5;
    }
    std::cout << "LEVEL1 interp best direction_op = " << direction_op << " , best ratio = " << best_interp_ratio << std::endl;

    meta = result_meta.ratio > best_interp_ratio && result_meta.ratio < 80 && best_interp_ratio < 80;
    printf("LEVEL2 meta Compression Ratio = %.2f\n", result_meta.ratio);
    printf("LEVEL2 interp Compression Ratio = %.2f\n", best_interp_ratio);
    printf("LEVEL2 tuning choose %s\n", meta ? "meta" : "interp");

    if (meta) {
        int capacity = 65536 * 2;
        META::optimize_quant_invl_3d(data.get(), dims[0], dims[1], dims[2], eb, capacity);
        meta_params.capacity = capacity;
//    printf("tuning capacity: %d\n", capacity);


        meta_params.prediction_dim = 2;
        result_meta = meta_compress_3d(sampling_data.data(), sampling_num, r1, r2, r3, eb, meta_params);
        printf("LEVEL1 meta ratio = %.2f, lorenzo:%d, lorenzo2:%d, pred_dim:%d compress_time:%.3f\n",
               result_meta.ratio, meta_params.use_lorenzo, meta_params.use_lorenzo_2layer, meta_params.prediction_dim,
               result_meta.compress_time);
        if (result_meta.ratio > best_meta_ratio * 1.02) {
            best_meta_ratio = result_meta.ratio;
        } else {
            meta_params.prediction_dim = 3;
        }

        if (reb < 1.01e-6 && best_interp_ratio > 5) {
            meta_params.capacity = 16384;
            result_meta = meta_compress_3d(sampling_data.data(), sampling_num, r1, r2, r3, eb, meta_params);
            printf("LEVEL1 meta ratio = %.2f, lorenzo:%d, lorenzo2:%d, pred_dim:%d compress_time:%.3f\n",
                   result_meta.ratio, meta_params.use_lorenzo, meta_params.use_lorenzo_2layer, meta_params.prediction_dim,
                   result_meta.compress_time);
            if (result_meta.ratio > best_meta_ratio * 1.02) {
                best_meta_ratio = result_meta.ratio;
            } else {
                meta_params.capacity = capacity;
            }
        }


        clock_gettime(CLOCK_REALTIME, &end);
        auto tuning_time = (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
        std::cout << "LEVEL2 Tuning time = " << tuning_time << "s" << std::endl;

        meta_params.print();
        auto result = meta_compress_decompress_3d(data.get(), num, dims[0], dims[1], dims[2], eb, meta_params, true);
        printf("LEVEL2 TUNING meta PSNR = %f, NRMSE = %.10G, Compression Ratio = %.2f\n", result.psnr, result.nrmse,
               result.ratio);
        std::cout << "LEVEL2 total compress time = " << tuning_time + result.compress_time << std::endl;
        std::cout << "LEVEL2 total decompress time = " << tuning_time + result.decompress_time << std::endl;

    } else {
        clock_gettime(CLOCK_REALTIME, &end);
        auto tuning_time = (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000;
        std::cout << "LEVEL2 Tuning time = " << tuning_time << "s" << std::endl;

        block_size = 6;
        interp_block_size = 32;
        auto result = interp_compress_decompress<N>(path, data.get(), num, eb, interp_level, interp_op, direction_op, block_size,
                                                    interp_block_size, sz_op, args...);
        printf("LEVEL2 TUNING interp PSNR = %f, NRMSE = %.10G, Compression Ratio = %.2f\n", result.psnr, result.nrmse,
               result.ratio);
        std::cout << "LEVEL2 total compress time = " << tuning_time + result.compress_time << std::endl;
        std::cout << "LEVEL2 total decompress time = " << tuning_time + result.decompress_time << std::endl;
    }
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
            interp_meta_tuning<3>(argv[1], reb, dims[0], dims[1], dims[2]);
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
        interp_compress_decompress_nowritefile<4>(argv[1], data.release(), num, eb, interp_level, interp_op, direction_op,
                                                  block_size,
                                                  interp_block_size, sz_op, dims[0], dims[1], dims[2], dims[3]);
    }


    return 0;
}
