/**
 * This is a demo for creating your own compressor in SZ3.
 * Essentially, there are two main steps: first, you need to implement your own compressor.
 * Please check out the 4 examples provided in the main function for implementing your own compressor.
 *
 * Second, you need to run your compressor with SZ3 executable or a new executable.
 * This demon code is a new executable, so you have control over the parameter parsing, IO, etc.
 * If you want to use the SZ3 executable to run your own compressor, please follow the steps below:
 * 1. Add a new ALGO in SZ3 Config
 * 2. Add a new hpp file in "include/SZ3/api/impl/" to assemble your compressor, one example easy to follow is
 * "include/SZ3/api/impl/SZAlgoNopred.hpp"
 * 3. Add the corresponding code in "include/SZ3/api/impl/SZDispatcher.hpp" to dispatch the new compressor.
 * 4. When executing the SZ3 executable, use -c to specify the config file. The config file should use the new ALGO.
 *  Example config file is "tools/sz3/sz3.config".
 */

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
    if (dim != 3) {
        printf("This demo only supports 3D data.\n");
        return 0;
    }
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }

    // create config with data dimensions
    Config conf;
    conf.setDims(dims.begin(), dims.end());

    // prepare input and output buffers
    std::vector<float> input_data(conf.num);
    SZ3::readfile<float>(argv[1], conf.num, input_data.data());
    std::vector<float> input_data_copy(input_data);
    std::vector<float> dec_data(conf.num);
    auto dec_data_pos = dec_data.data();
    std::vector<char> cmpData(conf.num * sizeof(float) * 2);
    size_t cmpSize = cmpData.size();

    // set error bound (you can also use relative error bound EB_REL etc.)
    conf.errorBoundMode = EB_ABS;
    conf.absErrorBound = atof(argv[argp++]);

    // ================================================================================
    // There are 4 customization examples, you just need to pick up one.

    // Example 1: using SZ3 API defined in SZ3/API/sz.hpp
    // SZ3 has some built-in compressors with fixed modules, for example, the interpolation compressor uses
    // interpolation decomposition, linear quantizer, and huffman encoder.
    // The API takes the SZ3::Config as parameter which allows some basic customization, such as interpolation
    // method (linear or cubic) and the quantization levels.
    SZ_compress<float>(conf, input_data.data(), cmpData.data(), cmpSize);
    SZ_decompress<float>(conf, cmpData.data(), cmpSize, dec_data_pos);

    // Example 2: using provided modules to assemble a generic compressor
    // The generic compressor follow the decomposition -> encoding -> lossless pipeline.
    // You can use any decomposition, encoder, and lossless modules to form such a compressor, including bypass
    // modules. This example assembles the interpolation compressor with interpolation decomposition, linear
    // quantizer, and huffman encoder. You can change any of the modules to get a new compressor that is not
    // built-in in SZ3.
    SZ3_interpolation_compress<float, 3>(conf, input_data.data(), cmpData.data(), cmpSize);
    SZ3_interpolation_decompress<float, 3>(conf, cmpData.data(), cmpSize, dec_data_pos);

    // Example 3: using new modules to assemble a generic compressor
    // This is an extension of Example 2, which uses a custom decomposition module instead of built-in ones.
    // MyDecomposition is the new decomposition module, which implements the DecompositionInterface.
    // We put the source code of MyDecomposition here for demonstration purpose.
    // You should put your new modules in the corresponding module folders (include/SZ3/decomposition for Decomposition
    // modules).
    SZ3_customized_compress<float, 3>(conf, input_data.data(), cmpData.data(), cmpSize);
    SZ3_customized_decompress<float, 3>(conf, cmpData.data(), cmpSize, dec_data_pos);

    // Example 4: assemble a specialized compressor
    // If your compressor doesn't follow the decomposition -> encoding -> lossless pipeline,
    // you can implement your own compressor by implementing the CompressorInterface
    // Please check out "SZ3/compressor/specialized/SZExaaltCompressor.hpp". It contains two separate encoding process.

    // ================================================================================

    // In the end, you can use the SZ3 API to verify the decompressed data.
    SZ3::verify<float>(input_data_copy.data(), dec_data_pos, conf.num);

    // write compressed data and decompressed data to the same folder as original data
    writefile((src_file_name + ".demo").data(), cmpData.data(), cmpSize);
    writefile((src_file_name + ".demo.out").data(), dec_data.data(), conf.num);

    return 0;
}
