#include "SZ3/api/sz.hpp"
/**
 * Slow Lorenzo
 * support any data dimension
 */
using namespace SZ3;

template<class T, uint N>
void SZ3_lorenzo_v2_compress(Config &conf, T *data, char *dst, size_t &outSize) {

    calAbsErrorBound(conf, data);

    auto sz = make_compressor_sz_iterate<T, N>( //compressor
            conf,
            LorenzoPredictor<T, N, 1>(conf.absErrorBound), //predictor
            LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2), //quantizer
            HuffmanEncoder<int>(),//encoder
            Lossless_zstd()); //lossless
    sz->compress(conf, data, (uchar *) dst, outSize);
}

template<class T, uint N>
void SZ3_lorenzo_v2_decompress(const Config &conf, const char *cmpData, size_t cmpSize, T *decData) {
    auto cmpDataPos = (uchar *) cmpData;
    auto sz = make_compressor_sz_iterate<T, N>( //compressor
            conf,
            LorenzoPredictor<T, N, 1>(conf.absErrorBound), //predictor
            LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2), //quantizer
            HuffmanEncoder<int>(),//encoder
            Lossless_zstd()); //lossless

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
    SZ3_lorenzo_v2_compress<float, 3>(conf, input_data.data(), cmpData.data(), cmpSize);
    SZ3_lorenzo_v2_decompress<float, 3>(conf, cmpData.data(), cmpSize, dec_data_p);

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
