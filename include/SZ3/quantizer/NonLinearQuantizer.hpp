/**
 * @file NonLinearQuantizer.hpp
 * @ingroup Quantizer
 */

#ifndef SZ3_NONLINEAR_QUANTIZER_HPP
#define SZ3_NONLINEAR_QUANTIZER_HPP

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
class NonLinearQuantizer : public concepts::QuantizerInterface<T, int> {
public:

    NonLinearQuantizer(double eb, int r = 32768) : error_bound(eb), radius(r) {
        assert(eb != 0);
        build_quant_table();
    }

    NonLinearQuantizer() : error_bound(1e-4), radius(32768) {
        build_quant_table();
    }

    double get_eb() const { return error_bound; }

    void set_eb(double eb) {
        error_bound = eb;
        build_quant_table();
    }

    std::pair<int, int> get_out_range() const override {
        return std::make_pair(0, 2 * radius);
    }

    void build_quant_table() {
        quant_table.resize(2 * radius);
        std::vector<double> pos_levels(radius);
        pos_levels[0] = error_bound; // Smallest positive level

        for (int i = 1; i < radius; i++) {
            double r = (double) i / (radius - 1);
            pos_levels[i] = error_bound + (2 * error_bound * (radius - 1)) * pow(r, 2);
        }

        for (int i = 1; i < radius; i++) {
            if (pos_levels[i] - pos_levels[i-1] > 2 * error_bound) {
                pos_levels[i] = pos_levels[i-1] + 2 * error_bound;
            }
        }

        for (int i = 0; i < radius; i++) {
            quant_table[radius + i] = pos_levels[i];
            quant_table[radius - 1 - i] = -pos_levels[i];
        }
    }

    // quantize the data with a prediction value, and returns the quantization index and the decompressed data
    // int quantize(T data, T pred, T& dec_data);
    ALWAYS_INLINE int quantize_and_overwrite(T &data, T pred) override {
        T diff = data - pred;
        int quant_index = search_best_quant_interval(diff);
        data = pred + quant_table[quant_index];
        return quant_index;
    }

    // recover the data using the quantization index
    ALWAYS_INLINE T recover(T pred, int quant_index) override {
        return pred + quant_table[quant_index];
    }

    int search_best_quant_interval(T diff) {
        int index = 0;
        if (diff > 0) {
            index = radius;
            while (index < 2 * radius && quant_table[index] < diff) {
                index++;
            }
            if (index == 2 * radius) {
                index--;
            } else if (index > radius) {
                if (std::abs(quant_table[index] - diff) > std::abs(quant_table[index - 1] - diff)) {
                    index--;
                }
            }
        } else {
            index = radius;
            while (index > 0 && quant_table[index] > diff) {
                index--;
            }
            if (index == 0) {
                // do nothing
            } else if (index < radius) {
                if (std::abs(quant_table[index] - diff) > std::abs(quant_table[index + 1] - diff)) {
                    index++;
                }
            }
        }
        return index;
    }

    void save(unsigned char *&c) const override {
        write(uid, c);
        write(this->error_bound, c);
        write(this->radius, c);
    }

    void load(const unsigned char *&c, size_t &remaining_length) override {
        uchar uid_read;
        read(uid_read, c, remaining_length);
        if (uid_read != uid) {
            throw std::invalid_argument("NonLinearQuantizer uid mismatch");
        }
        read(this->error_bound, c, remaining_length);
        read(this->radius, c, remaining_length);
        build_quant_table();
    }

    int force_save_unpred(T data) override {
        return 0;
    }

    void print() override {
        printf("NonLinearQuantizer, error_bound = %f, radius = %d\n", error_bound, radius);
    }

private:
    double error_bound;
    int radius;
    std::vector<T> quant_table;
    uchar uid = 0b11;

};

}
#endif
