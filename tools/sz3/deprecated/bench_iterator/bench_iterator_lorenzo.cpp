//
// Created by Kai Zhao on 1/25/23.
//

#include "SZ3/api/sz.hpp"

using namespace SZ;

/*
 * global padding
 */

template<class T, uint N>
uchar *compress(Config &conf, T *data, size_t &compressed_size) {

    std::vector<int> quant_inds;
    quant_inds.reserve(conf.num);
    LinearQuantizer<T> quantizer(conf.absErrorBound);

    //====================== create a copy of the input data with padding ==================
    int padding = 2;
    size_t ds0_ = (conf.dims[2] + padding) * (conf.dims[1] + padding);
    size_t ds1_ = (conf.dims[2] + padding);
    size_t ds0 = conf.dims[2] * conf.dims[1];
    size_t ds1 = conf.dims[2];

    size_t num_ = 1;
    for (auto &dim: conf.dims) {
        num_ *= (dim + padding);
    }
    std::vector<T> data_(num_);
    for (size_t i = 0; i < conf.dims[0]; i++) {
        for (size_t j = 0; j < conf.dims[1]; j++) {
            memcpy(&data_[(i + padding) * ds0_ + (j + padding) * ds1_ + padding], &data[i * ds0 + j * ds1], ds1 * sizeof(T));
        }
    }

    T *datap = &data_[padding * ds0_ + padding * ds1_ + padding];

    //=====================================================================================

    size_t bsize = 6;
    auto blocks = std::make_shared<SZ::multi_dimensional_range<T, N>>(datap, std::begin(conf.dims), std::end(conf.dims), bsize, 0);
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
        {
            // ================================================================
            //          This part will be in the Lorenzo predictor
            auto idx = block.get_global_index();
            for (size_t i = idx[0]; i < std::min(conf.dims[0], idx[0] + bsize); i++) {
                for (size_t j = idx[1]; j < std::min(conf.dims[1], idx[1] + bsize); j++) {
                    for (size_t k = idx[2]; k < std::min(conf.dims[2], idx[2] + bsize); k++) {
                        size_t offset_ = i * ds0_ + j * ds1_ + k;
                        T *data_pos = &datap[offset_];
                        T pred = data_pos[-1] + data_pos[-ds1_] + data_pos[-ds0_]
                                 - data_pos[-ds1_ - 1] - data_pos[-ds0_ - 1]
                                 - data_pos[-ds0_ - ds1_] + data_pos[-ds0_ - ds1_ - 1];
                        quant_inds.push_back(quantizer.quantize_and_overwrite(datap[offset_], pred));
                    }
                }
            }
            // ================================================================
        }
    }

    HuffmanEncoder<int> encoder;
    encoder.preprocess_encode(quant_inds, 0);
    size_t bufferSize = 1.2 * (encoder.size_est() + sizeof(T) * quant_inds.size() + quantizer.size_est());

    uchar *buffer = new uchar[bufferSize];
    uchar *buffer_pos = buffer;

    quantizer.save(buffer_pos);

    encoder.save(buffer_pos);
    encoder.encode(quant_inds, buffer_pos);
    encoder.postprocess_encode();
    assert(buffer_pos - buffer < bufferSize);

    Lossless_zstd lossless;
    uchar *lossless_data = lossless.compress(buffer, buffer_pos - buffer, compressed_size);
    lossless.postcompress_data(buffer);

    return lossless_data;
}

template<class T, uint N>
void decompress(Config &conf, uchar const *cmpData, const size_t &cmpSize, T *decData) {
    size_t remaining_length = cmpSize;
    LinearQuantizer<T> quantizer(conf.absErrorBound);

    Lossless_zstd lossless;
    auto compressed_data = lossless.decompress(cmpData, remaining_length);
    uchar const *compressed_data_pos = compressed_data;

    quantizer.load(compressed_data_pos, remaining_length);

    HuffmanEncoder<int> encoder;
    encoder.load(compressed_data_pos, remaining_length);
    auto quant_inds = encoder.decode(compressed_data_pos, conf.num);
    encoder.postprocess_decode();

    lossless.postdecompress_data(compressed_data);

    //====================== create a copy of the input data with padding ==================
    int padding = 2;
    size_t ds0_ = (conf.dims[2] + padding) * (conf.dims[1] + padding);
    size_t ds1_ = (conf.dims[2] + padding);
    size_t ds0 = conf.dims[2] * conf.dims[1];
    size_t ds1 = conf.dims[2];

    size_t num_ = 1;
    for (auto &dim: conf.dims) {
        num_ *= (dim + padding);
    }

    std::vector<T> data_(num_);
    T *datap = &data_[padding * ds0_ + padding * ds1_ + padding];
    int *quant_inds_pos = &quant_inds[0];
    //=======================================================================================

    size_t bsize = 6;
    auto blocks = std::make_shared<SZ::multi_dimensional_range<T, N>>(datap, std::begin(conf.dims), std::end(conf.dims), bsize, 0);
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
        {
            // ================================================================
            //          This part will be in the Lorenzo predictor
            auto idx = block.get_global_index();
            for (size_t i = idx[0]; i < std::min(conf.dims[0], idx[0] + bsize); i++) {
                for (size_t j = idx[1]; j < std::min(conf.dims[1], idx[1] + bsize); j++) {
                    for (size_t k = idx[2]; k < std::min(conf.dims[2], idx[2] + bsize); k++) {
                        size_t offset_ = i * ds0_ + j * ds1_ + k;
                        T *data_pos = &datap[offset_];
                        T pred = data_pos[-1] + data_pos[-ds1_] + data_pos[-ds0_]
                                 - data_pos[-ds1_ - 1] - data_pos[-ds0_ - 1]
                                 - data_pos[-ds0_ - ds1_] + data_pos[-ds0_ - ds1_ - 1];
                        *data_pos = quantizer.recover(pred, *(quant_inds_pos++));
                    }
                }
            }
            // ================================================================
        }
    }

    for (size_t i = 0; i < conf.dims[0]; i++) {
        for (size_t j = 0; j < conf.dims[1]; j++) {
            memcpy(&decData[i * ds0 + j * ds1], &data_[(i + padding) * ds0_ + (j + padding) * ds1_ + padding], ds1 * sizeof(T));
        }
    }
}

template<class T, uint N>
void estimate_compress(Config &conf, T *data) {
    conf.absErrorBound = 1e-2;

    Timer timer(true);
    size_t cmpr_size;
    auto cmpr_data = compress<T, N>(conf, data, cmpr_size);
    timer.stop("compress");

    T *decData = new T[conf.num];
    timer.start();
    decompress<T, N>(conf, cmpr_data, cmpr_size, decData);
    timer.stop("decompress");

    printf("CR= %.3f\n", conf.num * sizeof(T) * 1.0 / cmpr_size);
    SZ::verify(data, decData, conf.num);
};

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