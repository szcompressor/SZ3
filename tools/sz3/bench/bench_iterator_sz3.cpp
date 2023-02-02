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
    size_t quant_count = 0;

    Timer timer(true);
    auto element_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(
            data, std::begin(conf.dims), std::end(conf.dims), 1, 0);
    for (auto element = element_range->begin(); element != element_range->end(); ++element) {
        quant_inds[quant_count++] = quantizer.quantize_and_overwrite(*element, 0);
    }
    timer.stop("Baseline SZ3 style (iterator, inline member function)");
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