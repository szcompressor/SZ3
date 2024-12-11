/*
 * Based on SZBioMDDecomposition.hpp
 * \author: Magnus Lundborg
 */

#ifndef SZ3_SZBIOMDXTCBASED_FRONTEND
#define SZ3_SZBIOMDXTCBASED_FRONTEND

#include <list>

#include "Decomposition.hpp"
#include "SZ3/utils/Config.hpp"

namespace SZ3 {

template <class T, uint N, class Quantizer>
class SZBioMDXtcDecomposition : public concepts::DecompositionInterface<T, int, N> {
   public:
    SZBioMDXtcDecomposition(const Config &conf, Quantizer quantizer) : quantizer(quantizer), conf(conf) {
        if (N != 1 && N != 2 && N != 3) {
            throw std::invalid_argument("SZBioMDXtcDecomposition only support 1D, 2D or 3D data");
        }
    }

    std::vector<int> compress(const Config &conf, T *data) override {
        if (N <= 2) {
            return compressSingleFrame(data);
        } else {
            return compressMultiFrame(data);
        }
    }

    T *decompress(const Config &conf, std::vector<int> &quantData, T *decData) override {
        if (N <= 2) {
            return decompressSingleFrame(quantData, decData);
        } else {
            return decompressMultiFrame(quantData, decData);
        }
    }

    void save(uchar *&c) override {
        write(firstFillFrame_, c);
        write(fillValue_, c);
        quantizer.save(c);
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        read(firstFillFrame_, c, remaining_length);
        read(fillValue_, c, remaining_length);
        quantizer.load(c, remaining_length);
    }

    size_t size_est() override { return 0; }

    std::pair<int, int> get_out_range() override { return quantizer.get_out_range(); }

    size_t get_num_elements() const {
        if (N == 3) {
            return firstFillFrame_ * conf.dims[1] * conf.dims[2];
        }
        return conf.num;
    }

   private:
    std::vector<int> compressSingleFrame(T *data) {
        std::vector<int> quantData(conf.num);
        for (size_t i = 0; i < conf.num; i++) {
            quantData[i] = quantizer.quantize_and_overwrite(data[i], 0);
        }
        quantizer.postcompress_data();
        return quantData;
    }

    T *decompressSingleFrame(std::vector<int> &quantData, T *decData) {
        for (size_t i = 0; i < conf.num; i++) {
            decData[i] = quantizer.recover(0, quantData[i]);
        }
        quantizer.postdecompress_data();
        return decData;
    }

    /* Start from the last frame and look for frames filled with the same value.
     * Those frames do not need to be compressed - they can all just be filled. */
    std::tuple<size_t, T> findFillValueAndFirstFilledFrame(T *data, std::vector<size_t> dims) {
        size_t numDims = dims.size();
        if (numDims < 3) {
            return std::make_tuple(dims[0], 0);
        }
        size_t frameStride = dims[1] * dims[2];
        size_t firstFillFrame = dims[0];

        if (firstFillFrame == 0) {
            return std::make_tuple(firstFillFrame, 0);
        }
        /* Assume that the first value of the last frame is the potential fill value */
        T fillFrameValue = data[(dims[0] - 1) * frameStride];
        /* To simplify compression/decompression below, assume that the first frame needs
         * to be compressed. */
        for (size_t i = dims[0] - 1; i > 0; i--) {
            size_t idx = i * frameStride;
            bool allFrameValuesAreFillValue = true;

            for (size_t j = 0; j < dims[1] * dims[2]; j++) {
                size_t idy = idx + j;
                if (data[idy] != fillFrameValue) {
                    allFrameValuesAreFillValue = false;
                    break;
                }
            }
            if (allFrameValuesAreFillValue) {
                firstFillFrame = i;
            } else {
                break;
            }
        }
        return std::make_tuple(firstFillFrame, fillFrameValue);
    }

    /* This just converts float to integer based on the absolute error. */
    std::vector<int> compressMultiFrame(T *data) {
        auto dims = conf.dims;
        std::vector<size_t> stride({dims[1] * dims[2], dims[2], 1});

        /* Find out if the last frames are all filled with the same value. */
        std::tuple<size_t, T> fillValueSettings = findFillValueAndFirstFilledFrame(data, dims);
        firstFillFrame_ = std::get<0>(fillValueSettings);
        fillValue_ = std::get<1>(fillValueSettings);
        size_t lastFrame = std::min(dims[0], firstFillFrame_);
        std::vector<int> quantData(lastFrame * dims[1] * dims[2]);

        for (size_t i = 0; i < lastFrame; i++) {        // time
            for (size_t j = 0; j < dims[1]; j++) {      // atoms
                for (size_t k = 0; k < dims[2]; k++) {  // xyz
                    size_t idx = i * stride[0] + j * stride[1] + k;
                    quantData[idx] = quantizer.quantize_and_overwrite(data[idx], 0);
                }
            }
        }
        quantizer.postcompress_data();

        return quantData;
    }

    /* This just converts integer to float based on the absolute error. */
    T *decompressMultiFrame(std::vector<int> &quantData, T *decData) {
        // printf("Decompressing 3D.\n");
        auto dims = conf.dims;
        std::vector<size_t> stride({dims[1] * dims[2], dims[2], 1});

        size_t lastFrame = std::min(dims[0], firstFillFrame_);

        for (size_t i = 0; i < lastFrame; i++) {        // time
            for (size_t j = 0; j < dims[1]; j++) {      // atoms
                for (size_t k = 0; k < dims[2]; k++) {  // xyz
                    size_t idx = i * stride[0] + j * stride[1] + k;
                    decData[idx] = quantizer.recover(0, quantData[idx]);
                }
            }
        }
        quantizer.postdecompress_data();

        /* Fill frames at the end with the fill value. */
        for (size_t i = firstFillFrame_; i < dims[0]; i++) {
            size_t idx = i * stride[0];
            for (size_t j = 0; j < dims[1] * dims[2]; j++) {
                size_t idy = idx + j;
                decData[idy] = fillValue_;
            }
        }
        return decData;
    }

    Config conf;
    size_t firstFillFrame_;
    T fillValue_;
    Quantizer quantizer;
};

template <class T, uint N, class Quantizer>
SZBioMDXtcDecomposition<T, N, Quantizer> make_decomposition_biomdxtc(const Config &conf, Quantizer quantizer) {
    return SZBioMDXtcDecomposition<T, N, Quantizer>(conf, quantizer);
}

}  // namespace SZ3

#endif
