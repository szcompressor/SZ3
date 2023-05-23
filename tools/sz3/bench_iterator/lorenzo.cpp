//
// Created by Kai Zhao on 1/25/23.
//

#include "SZ3/api/sz.hpp"
#include "SZ3/utils/mddata.hpp"

using namespace SZ;

/*
 * global padding
 */

template<class T, uint N>
uchar *compress(Config &conf, T *data, size_t &compressed_size) {

    std::vector<int> quant_inds;
    quant_inds.reserve(conf.num);
    LinearQuantizer<T> quantizer(conf.absErrorBound);

    //=====================================================================================

    size_t bsize = 6;
    size_t padding = 2;
    auto mddata = std::make_shared<SZ::multi_dimensional_data<T, N>>(data, conf.dims, true, padding);
    auto block = mddata->block_iter(bsize);
    do {
        auto range = block->get_block_range();
        auto d = block->mddata;
        auto ds = mddata->get_dim_strides();
        if (N == 3) {
            for (size_t i = range[0].first; i < range[0].second; i++) {
                for (size_t j = range[1].first; j < range[1].second; j++) {
                    for (size_t k = range[2].first; k < range[2].second; k++) {
                        T *c = d->get_data(i, j, k);
                        T pred = c[-1] + c[-ds[1]] + c[-ds[0]]
                                 - c[-ds[1] - 1] - c[-ds[0] - 1]
                                 - c[-ds[0] - ds[1]] + c[-ds[0] - ds[1] - 1];
                        quant_inds.push_back(quantizer.quantize_and_overwrite(*c, pred));
                    }
                }
            }
        }
    } while (block->next());


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

    int *quant_inds_pos = &quant_inds[0];

    size_t bsize = 6;
    size_t padding = 2;
    auto mddata = std::make_shared<SZ::multi_dimensional_data<T, N>>(decData, conf.dims, false, padding);
    auto block = mddata->block_iter(bsize);
    do {
        auto range = block->get_block_range();
        auto d = block->mddata;
        auto ds = mddata->get_dim_strides();
        if (N == 3) {
            for (size_t i = range[0].first; i < range[0].second; i++) {
                for (size_t j = range[1].first; j < range[1].second; j++) {
                    for (size_t k = range[2].first; k < range[2].second; k++) {
                        T *c = d->get_data(i, j, k);
                        T pred = c[-1] + c[-ds[1]] + c[-ds[0]]
                                 - c[-ds[1] - 1] - c[-ds[0] - 1]
                                 - c[-ds[0] - ds[1]] + c[-ds[0] - ds[1] - 1];
                        *c = quantizer.recover(pred, *(quant_inds_pos++));
                    }
                }
            }
        }
    } while (block->next());

    mddata->copy_data_out(decData);
}

template<class T, uint N>
void estimate_compress(Config &conf, T *data) {
    conf.absErrorBound = 1e-2;

    Timer timer(true);
    size_t cmpr_size;
    auto cmpr_data = compress<T, N>(conf, data, cmpr_size);
    timer.stop("compress");
    fflush(stdout);

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