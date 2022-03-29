#ifndef _SZ_QOI_ENCODER_HPP
#define _SZ_QOI_ENCODER_HPP

#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/utils/ska_hash/unordered_map.hpp"
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>


namespace SZ {


    template<class T>
    class QoIEncoder : public concepts::EncoderInterface<T> {

    public:

        QoIEncoder() : data_encoder(), eb_encoder() {

        }

        ~QoIEncoder() {}

        /**
         * build huffman tree using bins
         * @param bins
         * @param stateNum is no longer needed
         */
        void preprocess_encode(const std::vector<T> &bins, int stateNum) {
            size_t num_elements = bins.size() / 2;
            eb_encoder.preprocess_encode(bins.data(), num_elements, stateNum);
            data_encoder.preprocess_encode(bins.data() + num_elements, num_elements, stateNum);
        }

        //save the huffman Tree in the compressed data
        uint save(uchar *&c) {
            auto s1 = eb_encoder.save(c);
            auto s2 = data_encoder.save(c);
            return s1 + s2;
        }

        //perform encoding
        size_t encode(const std::vector<T> &bins, uchar *&bytes) {
            size_t num_elements = bins.size() / 2;
            // printf("%d %d %d %d\n", bins[0], bins[1], bins[2], bins[3]);
            size_t s1 = eb_encoder.encode(bins.data(), num_elements, bytes);
            size_t s2 = data_encoder.encode(bins.data() + num_elements, num_elements, bytes);
            // printf("eb_size = %ld, data_size = %ld\n", s1, s2);
            return s1 + s2;
        }

        void postprocess_encode() {
            eb_encoder.postprocess_encode();
            data_encoder.postprocess_encode();
        }

        void preprocess_decode() {}

        //perform decoding
        std::vector<T> decode(const uchar *&bytes, size_t targetLength) {
            std::vector<T> eb_inds = eb_encoder.decode(bytes, targetLength);
            std::vector<T> data_inds = data_encoder.decode(bytes, targetLength);
            eb_inds.insert(eb_inds.end(), data_inds.begin(), data_inds.end());
            return eb_inds;
        }

        void postprocess_decode() {
            eb_encoder.postprocess_decode();
            data_encoder.postprocess_decode();
        }

        //load Huffman tree
        void load(const uchar *&c, size_t &remaining_length) {
            eb_encoder.load(c, remaining_length);
            data_encoder.load(c, remaining_length);
            loaded = true;
        }

        bool isLoaded() { return loaded; }

    private:
        HuffmanEncoder<T> eb_encoder;
        HuffmanEncoder<T> data_encoder;
        bool loaded = false;

    };
}

#endif