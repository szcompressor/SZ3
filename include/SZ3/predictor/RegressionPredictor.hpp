#ifndef _SZ_REGRESSION_PREDICTOR_HPP
#define _SZ_REGRESSION_PREDICTOR_HPP

#include <iostream>

#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"

namespace SZ3 {

// N-d regression predictor
template <class T, uint N>
class RegressionPredictor : public concepts::PredictorInterface<T, N> {
   public:
    using block_iter = typename block_data<T, N>::block_iterator;

    static const uint8_t predictor_id = 0b00000010;

    RegressionPredictor() : quantizer_independent(0), quantizer_liner(0), prev_coeffs{0}, current_coeffs{0} {}

    RegressionPredictor(uint block_size, double eb)
        : quantizer_independent(eb / (N + 1)),
          quantizer_liner(eb / (N + 1) / block_size),
          prev_coeffs{0},
          current_coeffs{0} {}

    bool precompress(const block_iter &block) override {
        auto range = block.get_block_range();

        std::array<double, N> dims{0};
        double num_elements = 1;
        for (int i = 0; i < N; i++) {
            dims[i] = range[i].second - range[i].first;
            if (dims[i] <= 1) {
                return false;
            }
            num_elements *= dims[i];
        }

        std::array<double, N + 1> sum{0};
        block_iter::foreach (block, [&](T *c, const std::array<size_t, N> &index) {
            for (int i = 0; i < N; i++) {
                sum[i] += index[i] * (*c);
            }
            sum[N] += *c;
        });
        std::fill(current_coeffs.begin(), current_coeffs.end(), 0);
        current_coeffs[N] = sum[N] / num_elements;
        for (int i = 0; i < N; i++) {
            current_coeffs[i] = (2 * sum[i] / (dims[i] - 1) - sum[N]) * 6 / num_elements / (dims[i] + 1);
            current_coeffs[N] -= (dims[i] - 1) * current_coeffs[i] / 2;
        }
        return true;
    }

    void precompress_block_commit() noexcept override {
        pred_and_quantize_coefficients();
        std::copy(current_coeffs.begin(), current_coeffs.end(), prev_coeffs.begin());
    }

    bool predecompress(const block_iter &block) override {
        auto range = block.get_block_range();
        for (const auto &r : range) {
            if (r.second - r.first <= 1) {
                return false;
            }
        }
        pred_and_recover_coefficients();
        return true;
    }

    ALWAYS_INLINE T estimate_error(const block_iter &block, T *d, const std::array<size_t, N> &index) override {
        return fabs(*d - predict(block, d, index));
    }

    ALWAYS_INLINE T predict(const block_iter &block, T *d, const std::array<size_t, N> &index) override {
        if constexpr (N == 1) {
            return current_coeffs[0] * index[0] + current_coeffs[1];
        } else if constexpr (N == 2) {
            return current_coeffs[0] * index[0] + current_coeffs[1] * index[1] + current_coeffs[2];
        } else if constexpr (N == 3) {
            return current_coeffs[0] * index[0] + current_coeffs[1] * index[1] + current_coeffs[2] * index[2] +
                   current_coeffs[3];
        } else if constexpr (N == 4) {
            return current_coeffs[0] * index[0] + current_coeffs[1] * index[1] + current_coeffs[2] * index[2] +
                   current_coeffs[3] * index[3] + current_coeffs[4];
        } else {
            static_assert(N <= 4, "Unsupported dimension or layer configuration");
            return T(0);
        }
    }

    void save(uchar *&c) override {
        write(regression_coeff_quant_inds.size(), c);
        if (!regression_coeff_quant_inds.empty()) {
            quantizer_independent.save(c);
            quantizer_liner.save(c);
            HuffmanEncoder<int> encoder = HuffmanEncoder<int>();
            encoder.preprocess_encode(
                regression_coeff_quant_inds,
                std::max(quantizer_independent.get_out_range().second, quantizer_liner.get_out_range().second));
            encoder.save(c);
            encoder.encode(regression_coeff_quant_inds, c);
            encoder.postprocess_encode();
        }
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        size_t coeff_size = 0;
        read(coeff_size, c, remaining_length);
        if (coeff_size > 0) {
            quantizer_independent.load(c, remaining_length);
            quantizer_liner.load(c, remaining_length);
            HuffmanEncoder<int> encoder = HuffmanEncoder<int>();
            encoder.load(c, remaining_length);
            regression_coeff_quant_inds = encoder.decode(c, coeff_size);
            encoder.postprocess_decode();
            remaining_length -= coeff_size * sizeof(int);
            std::fill(current_coeffs.begin(), current_coeffs.end(), 0);
            regression_coeff_index = 0;
        }
    }

    void print() const override {
        std::cout << "Regression predictor, indendent term eb = " << quantizer_independent.get_eb() << "\n";
        std::cout << "Regression predictor, linear term eb = " << quantizer_liner.get_eb() << "\n";
        // int count = 0;
        // int ind = regression_coeff_index ? regression_coeff_index : regression_coeff_quant_inds.size();
        std::cout << "Prev coeffs: ";
        for (const auto &c : prev_coeffs) {
            std::cout << c << " ";
        }
        std::cout << "\nCurrent coeffs: ";
        for (const auto &c : current_coeffs) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
    }

   private:
    LinearQuantizer<T> quantizer_liner, quantizer_independent;
    std::vector<int> regression_coeff_quant_inds;
    size_t regression_coeff_index = 0;
    std::array<T, N + 1> current_coeffs;
    std::array<T, N + 1> prev_coeffs;

    void pred_and_quantize_coefficients() {
        for (int i = 0; i < static_cast<int>(N); i++) {
            regression_coeff_quant_inds.push_back(
                quantizer_liner.quantize_and_overwrite(current_coeffs[i], prev_coeffs[i]));
        }
        regression_coeff_quant_inds.push_back(
            quantizer_independent.quantize_and_overwrite(current_coeffs[N], prev_coeffs[N]));
    }

    void pred_and_recover_coefficients() {
        for (int i = 0; i < static_cast<int>(N); i++) {
            current_coeffs[i] =
                quantizer_liner.recover(current_coeffs[i], regression_coeff_quant_inds[regression_coeff_index++]);
        }
        current_coeffs[N] =
            quantizer_independent.recover(current_coeffs[N], regression_coeff_quant_inds[regression_coeff_index++]);
    }
};

}  // namespace SZ3

#endif
