/**
 */

#include "SZ3/api/sz.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/compressor/ZFPCompressor.hpp"
#include "SZ3/decomposition/ZFPDecomposition.hpp"
#include "SZ3/encoder/ZFPEncoder.hpp"
#include "SZ3/utils/Statistic.hpp"

using namespace SZ3;

template <class T, class Int, uint N>
void test_zfp(const Config &conf, std::vector<T> &input_data) {
    // prepare compress and decompress buffers
    std::vector<T> dec_data(conf.num);
    std::vector<uchar> cmp_data(conf.num * sizeof(float) * 2);
    size_t cmp_cap = cmp_data.size();
    Timer timer;

    {
        printf("\n====== Testing ZFP standalone in FZ =========\n");
        auto input_data_copy = input_data;
        auto conf_copy = conf;

        auto zfp = make_compressor_zfp<T, N>();
        timer.start();
        auto cmpSize = zfp->compress(conf_copy, input_data_copy.data(), cmp_data.data(), cmp_cap);
        timer.stop("compression");

        timer.start();
        zfp->decompress(conf_copy, cmp_data.data(), cmpSize, dec_data.data());
        timer.stop("decompression");
        SZ3::verify<T>(input_data.data(), dec_data.data(), conf.num);
        printf("compression ratio = %.3f\n", static_cast<float>(conf.num * sizeof(float)) / cmpSize);
    }

    {
        printf("\n====== Testing ZFP modules in FZ =========\n");

        auto input_data_copy(input_data);
        Config conf_copy = conf;
        auto cmp_data_pos = cmp_data.data();

        ZFPDecomposition<T, Int, N> zfp_transform;
        ZFPEncoder<Int, N> zfp_encoder(conf);

        timer.start();
        auto transformed = zfp_transform.compress(conf_copy, input_data_copy.data());
        zfp_encoder.encode(transformed, cmp_data_pos);
        auto cmpSize = cmp_data_pos - cmp_data.data();
        timer.stop("compression");

        timer.start();
        const uchar *cmp_data_pos1 = cmp_data.data();
        transformed = zfp_encoder.decode(cmp_data_pos1, cmpSize);
        zfp_transform.decompress(conf_copy, transformed, dec_data.data());
        timer.stop("decompression");

        SZ3::verify<T>(input_data.data(), dec_data.data(), conf.num);
        printf("compression ratio = %.3f\n", static_cast<float>(conf.num * sizeof(T)) / cmpSize);
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "SZ v" << SZ3_VER << std::endl;
        std::cout << "usage: " << argv[0] << " data_file -num_dim dim0 .. dimn ABS" << std::endl;
        std::cout << "example: " << argv[0] << " qmcpack.dat -3 33120 69 69 1e-3" << std::endl;
        return 0;
    }

    std::string src_file_name = argv[1];

    // read dimensions
    int dim = atoi(argv[2] + 1);
    assert(1 <= dim && dim <= 4);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }

    // create config with data dimensions
    Config conf;
    conf.setDims(dims.begin(), dims.end());

    // set error bound (you can also use relative error bound EB_REL etc.)
    conf.errorBoundMode = EB_ABS;
    conf.absErrorBound = atof(argv[argp++]);

    // read data and convert error bound
    std::vector<float> input_data(conf.num);
    SZ3::readfile<float>(argv[1], conf.num, input_data.data());
    calAbsErrorBound(conf, input_data.data());

    if (dim == 1) {
        test_zfp<float, int32_t, 1>(conf, input_data);
    } else if (dim == 2) {
        test_zfp<float, int32_t, 2>(conf, input_data);
    } else if (dim == 3) {
        test_zfp<float, int32_t, 3>(conf, input_data);
    }
    return 0;
}
