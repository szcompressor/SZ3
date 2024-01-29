#include "SZ3/api/sz.hpp"

/**
 * Fast Lorenzo
 * Only support 1D and 3D
 */
using namespace SZ3;

template<class T, uint N>
void SZ3_lorenzo_v1_compress(Config &conf, T *data, char *dst, size_t &outSize) {

    calAbsErrorBound(conf, data);
    if (N != 1 && N != 3) {
        throw std::invalid_argument("Lorenzo v1 only support 1D and 3D data input");
    }
    conf.lorenzo = true;
    conf.lorenzo2 = false;
    conf.regression = false;
    conf.regression2 = false;
    auto sz = make_compressor_sz_generic<T, N>(
            make_decomposition_lorenzo_regression<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
            HuffmanEncoder<int>(),
            Lossless_zstd());

    sz->compress(conf, data, (uchar *) dst, outSize);
}

template<class T, uint N>
void SZ3_lorenzo_v1_decompress(const Config &conf, const char *cmpData, size_t cmpSize, T *decData) {
    uchar const *cmpDataPos = (uchar *) cmpData;
    auto sz = make_compressor_sz_generic<T, N>(
            make_decomposition_lorenzo_regression<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
            HuffmanEncoder<int>(),
            Lossless_zstd());

    sz->decompress(conf, cmpDataPos, cmpSize, decData);
}

int main(int argc, char **argv) {


    std::vector<size_t> dims({100, 200, 300});
    Config conf({dims[0], dims[1], dims[2]});
    conf.cmprAlgo = ALGO_INTERP_LORENZO;
    conf.errorBoundMode = EB_ABS; // refer to def.hpp for all supported error bound mode
    conf.absErrorBound = 1E-3; // absolute error bound 1e-3


    std::vector<float> input_data(conf.num);
    std::vector<float> dec_data(conf.num);
    std::vector<size_t> stride({dims[1] * dims[2], dims[2], 1});

    for (size_t i = 0; i < dims[0]; ++i) {
        for (size_t j = 0; j < dims[1]; ++j) {
            for (size_t k = 0; k < dims[2]; ++k) {
                double x = static_cast<double>(i) - static_cast<double>(dims[0]) / 2.0;
                double y = static_cast<double>(j) - static_cast<double>(dims[1]) / 2.0;
                double z = static_cast<double>(k) - static_cast<double>(dims[2]) / 2.0;
                input_data[i * stride[0] + j * stride[1] + k] = static_cast<float>(.0001 * y * sin(y) + .0005 * cos(pow(x, 2) + x) + z);
            }
        }
    }

    std::vector<float> input_data_copy(input_data);

    auto dec_data_p = dec_data.data();

    std::vector<char> cmpData(conf.num);
    size_t cmpSize = cmpData.size();
    SZ3_lorenzo_v1_compress<float, 3>(conf, input_data.data(), cmpData.data(), cmpSize);
    SZ3_lorenzo_v1_decompress<float, 3>(conf, cmpData.data(), cmpSize, dec_data_p);

    double max_err = 0.0;
    for (size_t i = 0; i < conf.num; i++) {
        if (fabs(dec_data[i] - input_data_copy[i]) > max_err) {
            max_err = fabs(dec_data[i] - input_data_copy[i]);
        }
    }
    printf("Smoke test %s", max_err <= conf.absErrorBound ? "passed" : "failed");
//    printf("%lu ", conf.num);
    return 0;


}
