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
    size_t bsize = 6;
    auto blocks = std::make_shared<SZ::multi_dimensional_range<T, N>>(
            data, std::begin(conf.dims), std::end(conf.dims), bsize, 0);
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
        auto idx = block.get_global_index();
        for (size_t i = idx[0]; i < ((idx[0] + bsize >= conf.dims[0]) ? conf.dims[0] : idx[0] + bsize); i++) {
            for (size_t j = idx[1]; j < ((idx[1] + bsize >= conf.dims[1]) ? conf.dims[1] : idx[1] + bsize); j++) {
                for (size_t k = idx[2]; k < ((idx[2] + bsize >= conf.dims[2]) ? conf.dims[2] : idx[2] + bsize); k++) {
                    size_t offset = i * conf.dims[1] * conf.dims[2] + j * conf.dims[2] + k;
                    //TODO force substitution for the function call, make it as fast as Hybrid (block iterator, function substituted)
                    quant_inds[offset] = quantizer.quantize_and_overwrite(data[offset], 0);
                }
            }
        }
    }

    timer.stop("V3 Hybrid (block iterator, inline member function)");


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