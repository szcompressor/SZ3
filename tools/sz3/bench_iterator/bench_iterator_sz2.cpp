//
// Created by Kai Zhao on 1/25/23.
//

#include "SZ3/api/sz.hpp"
#include "omp.h"

using namespace SZ;

template<class T, uint N>
void estimate_compress(Config conf, T *data) {
    std::vector<int> quant_inds(conf.num);
    LinearQuantizer<T> quantizer;

    Timer timer(true);
    std::vector<T> unpred;
    unpred.reserve(conf.num);
    double error_bound = quantizer.get_eb();
    double error_bound_reciprocal = 1 / quantizer.get_eb();
    int radius = quantizer.get_radius();

    for (size_t i = 0; i < conf.dims[0]; i++) {
        for (size_t j = 0; j < conf.dims[1]; j++) {
            for (size_t k = 0; k < conf.dims[2]; k++) {
                size_t offset = i * conf.dims[1] * conf.dims[2] + j * conf.dims[2] + k;
                T diff = data[offset] - 0;
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
                    if (fabs(decompressed_data - data[offset]) > error_bound) {
                        unpred.push_back(data[offset]);
                        quant_inds[offset] = 0;
                    } else {
                        data[offset] = decompressed_data;
                        quant_inds[offset] = quant_index_shifted;
                    }
                } else {
                    unpred.push_back(data[offset]);
                    quant_inds[offset] = 0;
                }
            }
        }
    }
    timer.stop("Baseline SZ2 style (nested loop, function substituted)");

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