#include <compressor/SZInterpolationCompressor.hpp>
#include <compressor/SZBlockInterpolationCompressor.hpp>
#include <quantizer/IntegerQuantizer.hpp>
#include <predictor/ComposedPredictor.hpp>
#include <predictor/SimplePredictor.hpp>
#include <lossless/Lossless_zstd.hpp>
#include <meta/meta_compress.hpp>
#include <utils/Iterator.hpp>
#include <utils/Verification.hpp>
#include <utils/Extraction.hpp>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <memory>
#include <type_traits>
#include <sstream>


template<class T, uint N, class ... Dims>
META::meta_compress_info
interp_compress_decompress(char *path, T *data, size_t num, double eb, int interp_op,
                           int direction_op, int block_size, int interp_block_size, Dims ... args) {
    std::string compressed_file_name(path);
    META::meta_compress_info compressInfo;

    {

        std::cout << "****************** compression ****************" << std::endl;
        std::cout << "Interp Op          = " << interp_op << std::endl
                  << "Direction          = " << direction_op << std::endl
                  << "SZ block size      = " << block_size << std::endl
                  << "Interp block size  = " << interp_block_size << std::endl;
//                  << "SZ_mode            = " << sz_op << std::endl;

        SZ::Timer timer(true);
        size_t compressed_size = 0;
        std::unique_ptr<unsigned char[]> compressed;


        auto dims = std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...};
        auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer<T>(eb),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd(),
                dims,
                interp_block_size,
                interp_op,
                direction_op
        );
        compressed.reset(sz.compress(data, compressed_size));


        double compress_time = timer.stop();
        auto compression_ratio = num * sizeof(T) * 1.0 / compressed_size;
        compressInfo.compress_time = compress_time;
        std::cout << "Compression time = " << compress_time << "s" << std::endl;
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

        SZ::Timer timer(true);
        std::unique_ptr<T[]> dec_data;

        auto dims = std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...};
        auto sz = SZ::SZInterpolationCompressor<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                SZ::LinearQuantizer<T>(eb),
                SZ::HuffmanEncoder<int>(),
                SZ::Lossless_zstd(),
                dims,
                interp_block_size,
                interp_op,
                direction_op
        );
        dec_data.reset(sz.decompress(compressed.get(), compressed_size));


        compressInfo.decompress_time = timer.stop();
        std::cout << "Decompression time = " << compressInfo.decompress_time << "s" << std::endl;

//        std::string decompressed_file_name(path);
//        std::stringstream ss;
//        ss << decompressed_file_name.substr(decompressed_file_name.rfind('/') + 1)
//           << ".sz3.out";
//        decompressed_file_name = ss.str();
//        std::cout << "DEBUG decompressed file = " << decompressed_file_name << std::endl;
//        SZ::writefile(decompressed_file_name.c_str(), dec_data.get(), num);

        size_t num1 = 0;
        auto ori_data = SZ::readfile<T>(path, num1);
        assert(num1 == num);
        double psnr, nrmse;
        SZ::verify<T>(ori_data.get(), dec_data.get(), num, psnr, nrmse);

        auto compression_ratio = num * sizeof(T) * 1.0 / compressed_size;
        printf("PSNR = %f, NRMSE = %.10G, Compression Ratio = %.2f\n", psnr, nrmse, compression_ratio);
        std::cout << "Options " << compression_ratio
                  << " " << interp_op
                  << " " << direction_op
                  << " " << block_size
                  << " " << interp_block_size
                  << std::endl;
        compressInfo.psnr = psnr;
        compressInfo.nrmse = nrmse;
        compressInfo.ratio = compression_ratio;
        return compressInfo;

    }

}

template<class T, uint N>
double interp_compress_block_version(T *data, std::array<size_t, N> dims, size_t num, double eb, int interp_level, int interp_op, int direction_op,
                                     int block_size, int interp_block_size) {

    std::cout << "****************** compression ****************" << std::endl;
    std::cout << "Interp Op          = " << interp_op << std::endl
              << "Direction          = " << direction_op << std::endl
              << "SZ block size      = " << block_size << std::endl
              << "Interp block size  = " << interp_block_size << std::endl;
//              << "SZ_mode            = " << sz_op << std::endl;

    SZ::Timer timer(true);

    std::vector<T> data1(data, data + num);
    size_t compressed_size = 0;

    SZ::Config<T, N> conf(eb, dims);
    conf.block_size = block_size;
    conf.stride = conf.block_size;
    auto sz = SZ::make_sz_fast_block_interpolation_compressor(
            conf,
            SZ::SimplePredictor<T, N>(eb),
            SZ::LinearQuantizer<T>(eb),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd(),
            interp_op,
            direction_op,
            interp_level
    );

    sz.compress(data1.data(), compressed_size);

    double compression_time = timer.stop();

    auto compression_ratio = num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compressed size = " << compressed_size << std::endl;
    std::cout << "Compression ratio = " << compression_ratio << std::endl;
    std::cout << "Tuning Interp Compression time = " << compression_time
              << " Ratio = " << compression_ratio
              << " Params = " << interp_level
              << " " << interp_op
              << " " << direction_op
              << " " << block_size
              << " " << interp_block_size
              //              << " " << sz_op
              << std::endl;


    std::cout << "****************** end ****************" << std::endl;
    return compression_ratio;
}

template<typename T>
META::meta_compress_info lorenzo_compress_3D(T *data, size_t num_elements, int r1, int r2, int r3, float precision,
                                             META::meta_params params) {
    size_t result_size = 0;
    META::meta_compress_info compressInfo;

    SZ::Timer timer(true);

    unsigned char *result = META::meta_compress_3d<T>(data, r1, r2, r3, precision, result_size, params, compressInfo);
    unsigned char *result_after_lossless = NULL;
    size_t lossless_outsize = META::meta_lossless_compress(ZSTD_COMPRESSOR, 3, result, result_size,
                                                           &result_after_lossless);
    free(result);
    free(result_after_lossless);

    compressInfo.compress_time = timer.stop();
    compressInfo.ori_bytes = num_elements * sizeof(T);
    compressInfo.compress_bytes = lossless_outsize;
    compressInfo.ratio = compressInfo.ori_bytes * 1.0 / compressInfo.compress_bytes;

    return compressInfo;
}


template<class T, uint N, class ... Dims>
void interp_lorenzo_tuning(char *path, double reb, bool enable_lorenzo, Dims ... args) {
    std::cout << "================================ BEGIN TUNING ================================" << std::endl;

    size_t num = 0;
    auto data = SZ::readfile<T>(path, num);
    double eb = reb * SZ::data_range(data.get(), num);
    auto dims = std::array<size_t, N>{static_cast<size_t>(std::forward<Dims>(args))...};


    SZ::Timer timer(true);

    size_t sampling_num, sampling_block;
    std::array<size_t, N> sample_dims;
    std::vector<T> sampling_data = SZ::sampling<T, N>(data.get(), dims, sampling_num, sample_dims, sampling_block);
    printf("%lu %lu %lu %lu %lu\n", sampling_data.size(), sampling_num, sample_dims[0], sample_dims[1], sample_dims[2]);

    META::meta_compress_info result_lorenzo;
    META::meta_params lorenzo_params(false, 5, 3, 0, true, true, false, eb);

    if (enable_lorenzo) {
        if (N != 3) {
            printf("Lorenzo can only be enabled in 3D mode");
            exit(0);
        }
        lorenzo_params.capacity = 65536 * 2;
        lorenzo_params.lossless = true;

        result_lorenzo = lorenzo_compress_3D(sampling_data.data(), sampling_num, sample_dims[0], sample_dims[1], sample_dims[2], eb, lorenzo_params);
        printf("Tuning lorenzo ratio = %.2f, lorenzo:%d, lorenzo2:%d, pred_dim:%d compress_time:%.3f\n",
               result_lorenzo.ratio, lorenzo_params.use_lorenzo, lorenzo_params.use_lorenzo_2layer, lorenzo_params.prediction_dim,
               result_lorenzo.compress_time);
    }
    double best_lorenzo_ratio = result_lorenzo.ratio;
    double best_interp_ratio = 0, ratio;

    int interp_level = -1, interp_op, direction_op = 0, block_size = sampling_block, interp_block_size = sampling_block;
//    SZ_Option sz_op = SZ_LR;
    for (int i = 0; i < 2; i++) {
        ratio = interp_compress_block_version<T, N>(sampling_data.data(), sample_dims, sampling_num, eb, interp_level, i, direction_op,
                                                    block_size, interp_block_size);
        if (ratio > best_interp_ratio) {
            best_interp_ratio = ratio;
            interp_op = i;
        }
    }
    std::cout << "interp best interp_op = " << interp_op << " , best ratio = " << best_interp_ratio << std::endl;

    ratio = interp_compress_block_version<T, N>(sampling_data.data(), sample_dims, sampling_num, eb, interp_level, interp_op, 5,
                                                block_size, interp_block_size);
    if (ratio > best_interp_ratio * 1.02) {
        best_interp_ratio = ratio;
        direction_op = 5;
    }
    std::cout << "interp best direction_op = " << direction_op << " , best ratio = " << best_interp_ratio
              << std::endl;

    bool lorenzo = enable_lorenzo && result_lorenzo.ratio > best_interp_ratio && result_lorenzo.ratio < 80 && best_interp_ratio < 80;
    printf("Lorenzo compression ratio = %.2f\n", result_lorenzo.ratio);
    printf("Interp compression ratio = %.2f\n", best_interp_ratio);
    printf("choose %s\n", lorenzo ? "lorenzo" : "interp");

    if (lorenzo) {
        int capacity = 65536 * 2;
        META::optimize_quant_invl_3d(data.get(), dims[0], dims[1], dims[2], eb, capacity);
        lorenzo_params.capacity = capacity;
//    printf("tuning capacity: %d\n", capacity);

        lorenzo_params.prediction_dim = 2;
        result_lorenzo = lorenzo_compress_3D(sampling_data.data(), sampling_num, sample_dims[0], sample_dims[1], sample_dims[2], eb, lorenzo_params);
        printf("Tuning Lorenzo ratio = %.2f, lorenzo:%d, lorenzo2:%d, pred_dim:%d compress_time:%.3f\n",
               result_lorenzo.ratio, lorenzo_params.use_lorenzo, lorenzo_params.use_lorenzo_2layer, lorenzo_params.prediction_dim,
               result_lorenzo.compress_time);
        if (result_lorenzo.ratio > best_lorenzo_ratio * 1.02) {
            best_lorenzo_ratio = result_lorenzo.ratio;
        } else {
            lorenzo_params.prediction_dim = 3;
        }

        if (reb < 1.01e-6 && best_lorenzo_ratio > 5) {
            lorenzo_params.capacity = 16384;
            result_lorenzo = lorenzo_compress_3D(sampling_data.data(), sampling_num, sample_dims[0], sample_dims[1], sample_dims[2], eb,
                                                 lorenzo_params);
            printf("Tuning Lorenzo ratio = %.2f, lorenzo:%d, lorenzo2:%d, pred_dim:%d compress_time:%.3f\n",
                   result_lorenzo.ratio, lorenzo_params.use_lorenzo, lorenzo_params.use_lorenzo_2layer,
                   lorenzo_params.prediction_dim,
                   result_lorenzo.compress_time);
            if (result_lorenzo.ratio > best_lorenzo_ratio * 1.02) {
                best_lorenzo_ratio = result_lorenzo.ratio;
            } else {
                lorenzo_params.capacity = capacity;
            }
        }


        double tuning_time = timer.stop();
        std::cout << "Tuning time = " << tuning_time << "s" << std::endl;
        std::cout << "================================ END TUNING ================================" << std::endl << std::endl;
        std::cout << "================================ BEGIN SZ-Lorenzo ================================" << std::endl;

        lorenzo_params.print();
        auto result = meta_compress_decompress_3d(data.get(), num, dims[0], dims[1], dims[2], eb, lorenzo_params, true);
        printf("PSNR = %f, NRMSE = %.10G, Compression Ratio = %.2f\n", result.psnr, result.nrmse,
               result.ratio);
        std::cout << "Total compress time = " << tuning_time + result.compress_time << std::endl;
        std::cout << "Total decompress time = " << result.decompress_time << std::endl;
        std::cout << "==================================== END SZ-Lorenzo ===================================" << std::endl;

    } else {
        double tuning_time = timer.stop();
        std::cout << "Tuning time = " << tuning_time << "s" << std::endl;
        std::cout << "====================================== END TUNING ======================================" << std::endl << std::endl;
        std::cout << "==================================== BEGIN SZ-Interp ===================================" << std::endl;

        block_size = 6;
        interp_block_size = 32;
        auto result = interp_compress_decompress<T, N>(path, data.get(), num, eb, interp_op, direction_op,
                                                       block_size, interp_block_size, args...);
//        printf("PSNR = %f, NRMSE = %.10G, Compression Ratio = %.2f\n", result.psnr, result.nrmse,
//               result.ratio);
        std::cout << "Total compress time = " << tuning_time + result.compress_time << std::endl;
        std::cout << "Total decompress time = " << result.decompress_time << std::endl;
        std::cout << "==================================== END SZ-Interp ===================================" << std::endl;
    }
}