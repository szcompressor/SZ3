#include "SZ3/api/sz.hpp"

using namespace SZ3;

template <class T, uint N>
void SZ3_interpolation_compress(Config &conf, T *data, char *dst, size_t &outSize) {
    calAbsErrorBound(conf, data);

    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());
    sz->compress(conf, data, reinterpret_cast<uchar *>(dst), outSize);
}

template <class T, uint N>
void SZ3_interpolation_decompress(const Config &conf, const char *cmpData, size_t cmpSize, T *decData) {
    auto cmpDataPos = reinterpret_cast<const uchar *>(cmpData);
    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());
    sz->decompress(conf, cmpDataPos, cmpSize, decData);
}

template <class T, uint N>
class MyDecomposition : public concepts::DecompositionInterface<T, int, N> {
   public:
    MyDecomposition(const Config &conf) : num(conf.num), quantizer(conf.absErrorBound) {}

    std::vector<int> compress(const Config &conf, T *data) override {
        // write your own logic here
        std::vector<int> output(num);
        for (size_t i = 0; i < num; i++) {
            output[i] = quantizer.quantize_and_overwrite(data[i], 0);  // this will replace the value of data[i]
        }
        return output;
    }

    T *decompress(const Config &conf, std::vector<int> &quant_inds, T *dec_data) override {
        // write your own logic here
        for (size_t i = 0; i < num; i++) {
            dec_data[i] = quantizer.recover(0, quant_inds[i]);
        }
        return dec_data;
    }

    void save(uchar *&c) override {
        write(num, c);
        quantizer.save(c);
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        read(num, c, remaining_length);
        quantizer.load(c, remaining_length);
    }

    std::pair<int, int> get_out_range() override { return std::make_pair(0, 0); }

   private:
    size_t num;
    LinearQuantizer<T> quantizer;
};

template <class T, uint N>
void SZ3_customized_compress(Config &conf, T *data, char *dst, size_t &outSize) {
    calAbsErrorBound(conf, data);

    auto sz = make_compressor_sz_generic<T, N>(MyDecomposition<T, N>(conf), HuffmanEncoder<int>(), Lossless_zstd());

    outSize = sz->compress(conf, data, reinterpret_cast<uchar *>(dst), outSize);
}

template <class T, uint N>
void SZ3_customized_decompress(const Config &conf, const char *cmpData, size_t cmpSize, T *decData) {
    auto cmpDataPos = reinterpret_cast<const uchar *>(cmpData);
    auto sz = make_compressor_sz_generic<T, N>(MyDecomposition<T, N>(conf), HuffmanEncoder<int>(), Lossless_zstd());

    sz->decompress(conf, cmpDataPos, cmpSize, decData);
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

    Config conf;
    conf.setDims(dims.begin(), dims.end());

    std::vector<float> input_data(conf.num);
    SZ3::readfile<float>(argv[1], conf.num, input_data.data());

    std::vector<float> input_data_copy(input_data);
    std::vector<float> dec_data(conf.num);
    auto dec_data_pos = dec_data.data();
    std::vector<char> cmpData(conf.num * sizeof(float) * 2);
    size_t cmpSize = cmpData.size();

    // set error bound
    conf.errorBoundMode = EB_ABS;
    conf.absErrorBound = atof(argv[argp++]);

    // Use this code to set relative error bound
    // conf.errorBoundMode=EB_REL;
    // conf.relErrorBound = atof(argv[argp++]);

    // compress, decompress, and verify
    if (dim == 3) {
        SZ_compress<float>(conf, input_data.data(), cmpData.data(), cmpSize);
        SZ_decompress<float>(conf, cmpData.data(), cmpSize, dec_data_pos);

        SZ3_interpolation_compress<float, 3>(conf, input_data.data(), cmpData.data(), cmpSize);
        SZ3_interpolation_decompress<float, 3>(conf, cmpData.data(), cmpSize, dec_data_pos);

        SZ3_customized_compress<float, 3>(conf, input_data.data(), cmpData.data(), cmpSize);
        SZ3_customized_decompress<float, 3>(conf, cmpData.data(), cmpSize, dec_data_pos);

        SZ3::verify<float>(input_data_copy.data(), dec_data_pos, conf.num);
    } else {
        printf("Data dimension not equal to 3 is not supported.\n");
        return 0;
    }

    // write compressed data and decompressed data to the same folder as original data
    writefile((src_file_name + ".demo").data(), cmpData.data(), cmpSize);
    writefile((src_file_name + ".demo.out").data(), dec_data.data(), conf.num);

    return 0;
}
