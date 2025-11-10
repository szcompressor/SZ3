#ifndef SZ3_LINEAR_QUANTIZER_HPP
#define SZ3_LINEAR_QUANTIZER_HPP

#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {
template <class T>
class LinearQuantizer : public concepts::QuantizerInterface<T, int> {
public:
    LinearQuantizer()
        : error_bound(1),
          error_bound_reciprocal(1),
          radius(32768) {
    }

    LinearQuantizer(double eb, int r = 32768, bool _strict_eb = true)
        : error_bound(eb),
          error_bound_reciprocal(1.0 / eb),
          radius(r),
          strict_eb(_strict_eb) {
        assert(eb != 0);
    }


    double get_eb() const { return error_bound; }

    void set_eb(double eb) {
        error_bound = eb;
        error_bound_reciprocal = 1.0 / eb;
    }

    std::pair<int, int> get_out_range() const override { return std::make_pair(0, radius * 2); }

    // quantize the data with a prediction value, and returns the quantization index and the decompressed data
    // int quantize(T data, T pred, T& dec_data);
    ALWAYS_INLINE int quantize_and_overwrite(T& data, T pred) override {
        T diff = data - pred;
        auto quant_index = static_cast<int64_t>(fabs(diff) * this->error_bound_reciprocal) + 1;
        if (quant_index < this->radius * 2) {
            quant_index >>= 1;
            int half_index = quant_index;
            quant_index <<= 1;
            int quant_index_shifted;
            if (diff < 0) {
                quant_index = -quant_index;
                quant_index_shifted = this->radius - half_index;
            } else {
                quant_index_shifted = this->radius + half_index;
            }
            T decompressed_data = pred + quant_index * this->error_bound;
            // if data is NaN, the diff is NaN, and NaN <= 0 is false
            diff = fabs(decompressed_data - data) - this->error_bound;
            if (diff <= 0 || (!strict_eb && diff <= get_floating_point_threshold(data))) {
                data = decompressed_data;
                return quant_index_shifted;
            } else {
                unpred.push_back(data);
                return 0;
            }
        } else {
            unpred.push_back(data);
            return 0;
        }
    }

    // recover the data using the quantization index
    ALWAYS_INLINE T recover(T pred, int quant_index) override {
        if (quant_index) {
            return recover_pred(pred, quant_index);
        } else {
            return recover_unpred();
        }
    }

    ALWAYS_INLINE T recover_pred(T pred, int quant_index) {
        return pred + 2 * (quant_index - this->radius) * this->error_bound;
    }

    ALWAYS_INLINE T recover_unpred() { return unpred[index++]; }

    ALWAYS_INLINE int force_save_unpred(T ori) override {
        unpred.push_back(ori);
        return 0;
    }

    size_t size_est() { return unpred.size() * sizeof(T); }

    void save(unsigned char*& c) const override {
        write(uid, c);
        write(this->error_bound, c);
        write(this->radius, c);
        size_t unpred_size = unpred.size();
        write(unpred_size, c);
        if (unpred_size > 0) {
            write(unpred.data(), unpred.size(), c);
        }
    }

    void load(const unsigned char*& c, size_t& remaining_length) override {
        uchar uid_read;
        read(uid_read, c, remaining_length);
        if (uid_read != uid) {
            throw std::invalid_argument("LinearQuantizer uid mismatch");
        }
        read(this->error_bound, c, remaining_length);
        this->error_bound_reciprocal = 1.0 / this->error_bound;
        read(this->radius, c, remaining_length);
        size_t unpred_size = 0;
        read(unpred_size, c, remaining_length);
        if (unpred_size > 0) {
            unpred.resize(unpred_size);
            read(unpred.data(), unpred_size, c, remaining_length);
        }
        index = 0;
    }

    void print() override {
        printf("[LinearQuantizer] error_bound = %.8G, radius = %d, unpred = %zu\n", error_bound, radius, unpred.size());
    }

private:
    /**
     * Get the floating-point precision threshold (2 * ULP) for the given value.
     * 
     * This function calculates 2 times the Unit in the Last Place (ULP) for the input value.
     * ULP is the smallest representable difference between two floating-point numbers
     * in the vicinity of the given value. 
     * This threshold is used to determine if small deviations from the error bound
     * are acceptable due to floating-point precision limitations.
     */
    T get_floating_point_threshold(T value) {
        if constexpr (std::is_same<T, float>::value) {
            union {
                float f;
                uint32_t u;
            } conv = {value};
            int biased_exp = (conv.u >> 23) & 0xFF;

            if (biased_exp == 0) {
                return 2.0f * 1.401298464324817e-45f; // 2 * 2^(-149) for denormals
            } else if (biased_exp == 255) {
                return 0.0f; // Infinity/NaN
            } else {
                int unbiased_exp = biased_exp - 127;
                int shift = unbiased_exp - 23 + 1; // +1 for *2
                return (shift >= 0) ? static_cast<float>(1U << shift) : 1.0f / static_cast<float>(1U << (-shift));
            }
        } else if constexpr (std::is_same<T, double>::value) {
            union {
                double d;
                uint64_t u;
            } conv = {value};
            int biased_exp = (conv.u >> 52) & 0x7FF;

            if (biased_exp == 0) {
                return 2.0 * 4.9406564584124654e-324; // 2 * 2^(-1074) for denormals
            } else if (biased_exp == 2047) {
                return 0.0; // Infinity/NaN
            } else {
                int unbiased_exp = biased_exp - 1023;
                int shift = unbiased_exp - 52 + 1; // +1 for *2
                return (shift >= 0) ? static_cast<double>(1ULL << shift) : 1.0 / static_cast<double>(1ULL << (-shift));
            }
        } else {
            return 0; // Non-floating-point types
        }
    }


    std::vector<T> unpred;
    size_t index = 0; // used in decompression only
    uchar uid = 0b10;

    double error_bound;
    double error_bound_reciprocal;
    int radius; // quantization interval radius
    bool strict_eb = true;
};
} // namespace SZ3
#endif