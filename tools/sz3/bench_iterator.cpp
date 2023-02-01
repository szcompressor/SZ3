//
// Created by Kai Zhao on 1/25/23.
//

#include "SZ3/api/sz.hpp"
#include "omp.h"

using namespace SZ;


template<class T>
static inline __attribute__((always_inline)) int quantize_and_overwrite(T &data, T pred, std::vector<T> &unpred, double error_bound,
                                                                        double error_bound_reciprocal, int radius) {
    T diff = data - pred;
    int quant_index = (int) (fabs(diff) * error_bound_reciprocal) + 1;
    if (quant_index < radius * 2) {
        quant_index >>= 1;
        int half_index = quant_index;
        quant_index <<= 1;
        int quant_index_shifted;
        if (diff < 0) {
            quant_index = -quant_index;
            quant_index_shifted = radius - half_index;
        } else {
            quant_index_shifted = radius + half_index;
        }
        T decompressed_data = pred + quant_index * error_bound;
        if (fabs(decompressed_data - data) > error_bound) {
            unpred.push_back(data);
            return 0;
        } else {
            data = decompressed_data;
            return quant_index_shifted;
        }
    } else {
        unpred.push_back(data);
        return 0;
    }
}

template<class T, uint N>
void estimate_compress(Config conf, T *data) {
    T *data1 = (T *) malloc(sizeof(T) * conf.num);
    T *data2 = (T *) malloc(sizeof(T) * conf.num);
    T *data3 = (T *) malloc(sizeof(T) * conf.num);
    T *data4 = (T *) malloc(sizeof(T) * conf.num);
    T *data5 = (T *) malloc(sizeof(T) * conf.num);
    T *data6 = (T *) malloc(sizeof(T) * conf.num);
    T *data7 = (T *) malloc(sizeof(T) * conf.num);
    T *data8 = (T *) malloc(sizeof(T) * conf.num);
    for (size_t i = 0; i < conf.num; i++) {
        data1[i] = data[i];
        data2[i] = data[i];
        data3[i] = data[i];
        data4[i] = data[i];
        data5[i] = data[i];
        data6[i] = data[i];
        data7[i] = data[i];
        data8[i] = data[i];
    }

    std::vector<int> quant_inds_1(conf.num);
    std::vector<int> quant_inds_2(conf.num);
    std::vector<int> quant_inds_3(conf.num);
    std::vector<int> quant_inds_4(conf.num);
    std::vector<int> quant_inds_5(conf.num);
    std::vector<int> quant_inds_6(conf.num);
    std::vector<int> quant_inds_7(conf.num);
    std::vector<int> quant_inds_8(conf.num);
    std::vector<T> unpred1;
    std::vector<T> unpred2;
    std::vector<T> unpred3;
    std::vector<T> unpred4;
    std::vector<T> unpred5;
    std::vector<T> unpred6;
    std::vector<T> unpred7;
    std::vector<T> unpred8;


    {
        LinearQuantizer<T> quantizer;
        size_t quant_count = 0;

        Timer timer(true);
        auto element_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                data1, std::begin(conf.dims), std::end(conf.dims), 1, 0);
        for (auto element = element_range->begin(); element != element_range->end(); ++element) {
            quant_inds_1[quant_count++] = quantizer.quantize_and_overwrite(*element, 0);
        }
        timer.stop("Baseline SZ3 style (iterator, inline member function)");
    }

    {
        LinearQuantizer<T> quantizer;

        Timer timer(true);
        unpred2.reserve(conf.num);
        double error_bound = quantizer.get_eb();
        double error_bound_reciprocal = 1 / quantizer.get_eb();
        int radius = quantizer.get_radius();

        for (size_t i = 0; i < conf.dims[0]; i++) {
            for (size_t j = 0; j < conf.dims[1]; j++) {
                for (size_t k = 0; k < conf.dims[2]; k++) {
                    size_t offset = i * conf.dims[1] * conf.dims[2] + j * conf.dims[2] + k;
                    T diff = data2[offset] - 0;
                    int quant_index = (int) (fabs(diff) * error_bound_reciprocal) + 1;
                    if (quant_index < radius * 2) {
                        quant_index >>= 1;
                        int half_index = quant_index;
                        quant_index <<= 1;
                        int quant_index_shifted;
                        if (diff < 0) {
                            quant_index = -quant_index;
                            quant_index_shifted = radius - half_index;
                        } else {
                            quant_index_shifted = radius + half_index;
                        }
                        T decompressed_data = 0 + quant_index * error_bound;
                        if (fabs(decompressed_data - data2[offset]) > error_bound) {
                            unpred2.push_back(data2[offset]);
                            quant_inds_2[offset] = 0;
                        } else {
                            data2[offset] = decompressed_data;
                            quant_inds_2[offset] = quant_index_shifted;
                        }
                    } else {
                        unpred2.push_back(data2[offset]);
                        quant_inds_2[offset] = 0;
                    }
                }
            }
        }
        timer.stop("Baseline SZ2 style (nested loop, function substituted)");
    }

    {
        LinearQuantizer<T> quantizer;

        Timer timer(true);
        unpred3.reserve(conf.num);

        size_t bsize = 6;
        double error_bound = quantizer.get_eb();
        double error_bound_reciprocal = 1 / quantizer.get_eb();
        int radius = quantizer.get_radius();
        unpred3.reserve(conf.num);
        auto blocks = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                data3, std::begin(conf.dims), std::end(conf.dims), bsize, 0);
        for (auto block = blocks->begin(); block != blocks->end(); ++block) {
            auto idx = block.get_global_index();
            for (size_t i = idx[0]; i < ((idx[0] + bsize >= conf.dims[0]) ? conf.dims[0] : idx[0] + bsize); i++) {
                for (size_t j = idx[1]; j < ((idx[1] + bsize >= conf.dims[1]) ? conf.dims[1] : idx[1] + bsize); j++) {
                    for (size_t k = idx[2]; k < ((idx[2] + bsize >= conf.dims[2]) ? conf.dims[2] : idx[2] + bsize); k++) {
                        size_t offset = i * conf.dims[1] * conf.dims[2] + j * conf.dims[2] + k;
//                        quant_inds_3[offset] = quantizer.quantize_and_overwrite(data[offset], 0);
                        T diff = data3[offset] - 0;
                        int quant_index = (int) (fabs(diff) * error_bound_reciprocal) + 1;
                        if (quant_index < radius * 2) {
                            quant_index >>= 1;
                            int half_index = quant_index;
                            quant_index <<= 1;
                            int quant_index_shifted;
                            if (diff < 0) {
                                quant_index = -quant_index;
                                quant_index_shifted = radius - half_index;
                            } else {
                                quant_index_shifted = radius + half_index;
                            }
                            T decompressed_data = 0 + quant_index * error_bound;
                            if (fabs(decompressed_data - data3[offset]) > error_bound) {
                                unpred3.push_back(data3[offset]);
                                quant_inds_3[offset] = 0;
                            } else {
                                data3[offset] = decompressed_data;
                                quant_inds_3[offset] = quant_index_shifted;
                            }
                        } else {
                            unpred3.push_back(data3[offset]);
                            quant_inds_3[offset] = 0;
                        }
                    }

                }
            }
        }

        timer.stop("V1 Hybrid (block iterator, function substituted)");
    }

    {
        LinearQuantizer<T> quantizer;
        unpred4.reserve(conf.num);

        Timer timer(true);
        size_t bsize = 6;
        double error_bound = quantizer.get_eb();
        double error_bound_reciprocal = 1 / quantizer.get_eb();
        int radius = quantizer.get_radius();
        auto blocks = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                data4, std::begin(conf.dims), std::end(conf.dims), bsize, 0);
        for (auto block = blocks->begin(); block != blocks->end(); ++block) {
            auto idx = block.get_global_index();
            for (size_t i = idx[0]; i < ((idx[0] + bsize >= conf.dims[0]) ? conf.dims[0] : idx[0] + bsize); i++) {
                for (size_t j = idx[1]; j < ((idx[1] + bsize >= conf.dims[1]) ? conf.dims[1] : idx[1] + bsize); j++) {
                    for (size_t k = idx[2]; k < ((idx[2] + bsize >= conf.dims[2]) ? conf.dims[2] : idx[2] + bsize); k++) {
                        size_t offset = i * conf.dims[1] * conf.dims[2] + j * conf.dims[2] + k;
                        //TODO force substitution for the function call, make it as fast as Hybrid (block iterator, function substituted)
//                        quant_inds_3[offset] = quantizer.quantize_and_overwrite(data[offset], 0);
                        quant_inds_4[offset] = quantize_and_overwrite<T>(data4[offset], 0, unpred4, error_bound, error_bound_reciprocal, radius);
                    }
                }
            }
        }

        timer.stop("V2 Hybrid (block iterator, inline function)");
    }

    {
        LinearQuantizer<T> quantizer;

        Timer timer(true);
        size_t bsize = 6;
        auto blocks = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                data5, std::begin(conf.dims), std::end(conf.dims), bsize, 0);
        for (auto block = blocks->begin(); block != blocks->end(); ++block) {
            auto idx = block.get_global_index();
            for (size_t i = idx[0]; i < ((idx[0] + bsize >= conf.dims[0]) ? conf.dims[0] : idx[0] + bsize); i++) {
                for (size_t j = idx[1]; j < ((idx[1] + bsize >= conf.dims[1]) ? conf.dims[1] : idx[1] + bsize); j++) {
                    for (size_t k = idx[2]; k < ((idx[2] + bsize >= conf.dims[2]) ? conf.dims[2] : idx[2] + bsize); k++) {
                        size_t offset = i * conf.dims[1] * conf.dims[2] + j * conf.dims[2] + k;
                        //TODO force substitution for the function call, make it as fast as Hybrid (block iterator, function substituted)
                        quant_inds_5[offset] = quantizer.quantize_and_overwrite(data5[offset], 0);
//                        quant_inds_3[offset] = quantize_and_overwrite<T>(data[offset], 0, unpred, error_bound, error_bound_reciprocal, radius);
                    }
                }
            }
        }

        timer.stop("V3 Hybrid (block iterator, inline member function)");
    }

    {
        LinearQuantizer<T> quantizer;

        Timer timer(true);
        unpred6.reserve(conf.num);

        size_t bsize = 6;
        auto blocks = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                data6, std::begin(conf.dims), std::end(conf.dims), bsize, 0);
        for (auto block = blocks->begin(); block != blocks->end(); ++block) {
            auto idx = block.get_global_index();
            for (size_t i = idx[0]; i < ((idx[0] + bsize >= conf.dims[0]) ? conf.dims[0] : idx[0] + bsize); i++) {
                for (size_t j = idx[1]; j < ((idx[1] + bsize >= conf.dims[1]) ? conf.dims[1] : idx[1] + bsize); j++) {
                    for (size_t k = idx[2]; k < ((idx[2] + bsize >= conf.dims[2]) ? conf.dims[2] : idx[2] + bsize); k++) {
                        size_t offset = i * conf.dims[1] * conf.dims[2] + j * conf.dims[2] + k;
                        //TODO force substitution for the function call, make it as fast as Hybrid (block iterator, function substituted)
                        quant_inds_6[offset] = quantizer.quantize_and_overwrite_no_this(data6[offset], 0, unpred6);
//                        quant_inds_3[offset] = quantize_and_overwrite<T>(data[offset], 0, unpred, error_bound, error_bound_reciprocal, radius);
                    }
                }
            }
        }

        timer.stop("V4 Hybrid (block iterator, inline member function no this )");
    }

    {
        LinearQuantizer<T> quantizer;

        Timer timer(true);
        size_t bsize = 6;
        unpred7.reserve(conf.num);

        auto blocks = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                data7, std::begin(conf.dims), std::end(conf.dims), bsize, 0);
        double error_bound = quantizer.get_eb();
        double error_bound_reciprocal = 1 / quantizer.get_eb();
        int radius = quantizer.get_radius();
        for (auto block = blocks->begin(); block != blocks->end(); ++block) {
            auto idx = block.get_global_index();
            for (size_t i = idx[0]; i < ((idx[0] + bsize >= conf.dims[0]) ? conf.dims[0] : idx[0] + bsize); i++) {
                for (size_t j = idx[1]; j < ((idx[1] + bsize >= conf.dims[1]) ? conf.dims[1] : idx[1] + bsize); j++) {
                    for (size_t k = idx[2]; k < ((idx[2] + bsize >= conf.dims[2]) ? conf.dims[2] : idx[2] + bsize); k++) {
                        size_t offset = i * conf.dims[1] * conf.dims[2] + j * conf.dims[2] + k;
                        //TODO force substitution for the function call, make it as fast as Hybrid (block iterator, function substituted)
                        quant_inds_7[offset] = quantizer.quantize_and_overwrite_no_this2(data7[offset], 0, unpred7, error_bound,
                                                                                         error_bound_reciprocal,
                                                                                         radius);
//                        quant_inds_3[offset] = quantize_and_overwrite<T>(data[offset], 0, unpred, error_bound, error_bound_reciprocal, radius);
                    }
                }
            }
        }

        timer.stop("V5 Hybrid (block iterator, inline member function no this 2)");
    }

    {
        LinearQuantizer<T> quantizer;

        Timer timer(true);
        unpred8.reserve(conf.num);

        size_t bsize = 6;
        auto blocks = std::make_shared<SZ::multi_dimensional_range<T, N>>(
                data8, std::begin(conf.dims), std::end(conf.dims), bsize, 0);
        double error_bound = quantizer.get_eb();
        double error_bound_reciprocal = 1 / quantizer.get_eb();
        int radius = quantizer.get_radius();
        for (auto block = blocks->begin(); block != blocks->end(); ++block) {
            auto idx = block.get_global_index();
            for (size_t i = idx[0]; i < ((idx[0] + bsize >= conf.dims[0]) ? conf.dims[0] : idx[0] + bsize); i++) {
                for (size_t j = idx[1]; j < ((idx[1] + bsize >= conf.dims[1]) ? conf.dims[1] : idx[1] + bsize); j++) {
                    for (size_t k = idx[2]; k < ((idx[2] + bsize >= conf.dims[2]) ? conf.dims[2] : idx[2] + bsize); k++) {
                        size_t offset = i * conf.dims[1] * conf.dims[2] + j * conf.dims[2] + k;
                        //TODO force substitution for the function call, make it as fast as Hybrid (block iterator, function substituted)
                        quant_inds_8[offset] = quantizer.quantize_and_overwrite_no_this2_static(data8[offset], 0, unpred8, error_bound,
                                                                                                error_bound_reciprocal,
                                                                                                radius);
//                        quant_inds_3[offset] = quantize_and_overwrite<T>(data[offset], 0, unpred, error_bound, error_bound_reciprocal, radius);
                    }
                }
            }
        }

        timer.stop("V6 Hybrid (block iterator, inline member function no this 2 static)");
    }


    for (size_t i = 0; i < conf.num; i++) {
        if (quant_inds_1[i] != quant_inds_2[i] || quant_inds_2[i] != quant_inds_3[i] || quant_inds_3[i] != quant_inds_4[i] ||
            quant_inds_4[i] != quant_inds_5[i] || quant_inds_5[i] != quant_inds_6[i] || quant_inds_6[i] != quant_inds_7[i] ||
            quant_inds_7[i] != quant_inds_8[i]) {
            printf("Mismatch at %lu\n", i);
            exit(0);
        }
    }
    {
        conf.cmprAlgo = ALGO_LORENZO_REG;
        conf.errorBoundMode = SZ::EB_ABS;
        conf.absErrorBound = 1e-2;

        Timer timer(true);
        size_t cmpSize = 0;
        auto cmpdata = SZ_compress(conf, data, cmpSize);
        timer.stop("Compression (SZ3 fast Lorenzo+Reg, ABS 1e-2)");
        delete[] cmpdata;

    }

    free(data1);
    free(data2);
    free(data3);
    free(data4);
    free(data5);
    free(data6);
    free(data7);
    free(data8);

}

int main(int argc, char *argv[]) {

    if (argc < 4) {
        std::cout << "usage: " << argv[0] << " data_file dim dim0 .. dimn" << std::endl;
        std::cout << "example: " << argv[0] << " qmcpack.dat 3 33120 69 69" << std::endl;
        return 0;
    }

    size_t num = 0;
    auto data = SZ::readfile<float>(argv[1], num);

    int dim = atoi(argv[2]);
    if (dim == 1) {
        Config config(atoi(argv[3]));
        estimate_compress<float, 1>(config, data.get());
    } else if (dim == 2) {
        Config config(atoi(argv[3]), atoi(argv[4]));
        estimate_compress<float, 2>(config, data.get());
    } else if (dim == 3) {
        Config config(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
        estimate_compress<float, 3>(config, data.get());
    } else if (dim == 4) {
        Config config(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
        estimate_compress<float, 4>(config, data.get());
    }

}