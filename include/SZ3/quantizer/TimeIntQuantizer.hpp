#ifndef SZ3_TIME_INT_QUANTIZER_HPP
#define SZ3_TIME_INT_QUANTIZER_HPP

#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>

namespace SZ3 {

template<class T>
class TimeIntQuantizer : public concepts::QuantizerInterface<T, int> {
public:

    TimeIntQuantizer(int pred_dim) : pred_dim(pred_dim) {
        quant_table.resize(2 * pred_dim + 1);
        for (int i = 0; i < quant_table.size(); i++) {
            quant_table[i] = i - pred_dim;
        }
    }

    TimeIntQuantizer() : pred_dim(4) {
        quant_table.resize(2 * pred_dim + 1);
        for (int i = 0; i < quant_table.size(); i++) {
            quant_table[i] = i - pred_dim;
        }
    }

    std::pair<int, int> get_out_range() const override {
        return std::make_pair(0, 2 * pred_dim);
    }

    // quantize the data with a prediction value, and returns the quantization index and the decompressed data
    ALWAYS_INLINE int quantize_and_overwrite(T &data, T pred) override {
        T diff = data - pred;
        int quant_index = diff + pred_dim;
        if (quant_index < 0) {
            quant_index = 0;
        } else if (quant_index > 2 * pred_dim) {
            quant_index = 2 * pred_dim;
        }
        data = pred + quant_table[quant_index];
        return quant_index;
    }

    // recover the data using the quantization index
    ALWAYS_INLINE T recover(T pred, int quant_index) override {
        return pred + quant_table[quant_index];
    }

    void save(unsigned char *&c) const override {
        write(uid, c);
        write(this->pred_dim, c);
    }

    void load(const unsigned char *&c, size_t &remaining_length) override {
        uchar uid_read;
        read(uid_read, c, remaining_length);
        if (uid_read != uid) {
            throw std::invalid_argument("TimeIntQuantizer uid mismatch");
        }
        read(this->pred_dim, c, remaining_length);
        quant_table.resize(2 * pred_dim + 1);
        for (int i = 0; i < quant_table.size(); i++) {
            quant_table[i] = i - pred_dim;
        }
    }

    int force_save_unpred(T data) override {
        return 0;
    }

    void print() override {
        printf("TimeIntQuantizer, pred_dim = %d\n", pred_dim);
    }

private:
    int pred_dim;
    std::vector<T> quant_table;
    uchar uid = 0b100;

};

}
#endif
