#ifndef _SZ_REGRESSION_PREDICTOR_HPP
#define _SZ_REGRESSION_PREDICTOR_HPP

#include "def.hpp"
#include "utils/Iterator.hpp"
#include "quantizer/Quantizer.hpp"
#include "encoder/Encoder.hpp"
#include <cstring>
#include <iostream>
#include <fstream>

namespace SZ {

// N-d regression predictor
    template<class T, uint N>
    class RegressionPredictor {
    public:
        static const uint8_t predictor_id = 0b00000010;

        RegressionPredictor() : quantizer(0), prev_coeffs{0}, current_coeffs{0} {}

        RegressionPredictor(T eb) : quantizer(eb), prev_coeffs{0}, current_coeffs{0} {}

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        void precompress_data(const iterator &) const noexcept {}

        void postcompress_data(const iterator &) const noexcept {}

        void predecompress_data(const iterator &) const noexcept {}

        void postdecompress_data(const iterator &) const noexcept {}

        inline T estimate_error(const iterator &iter) const noexcept {
            return fabs(*iter - predict(iter));
        }

        void precompress_block(const std::shared_ptr<Range> &range) noexcept {
            // std::cout << "precompress_block" << std::endl;
            std::array<size_t, N> dims;
            for (int i = 0; i < N; i++) {
                dims[i] = range->get_dimensions(i);
            }
            current_coeffs = compute_regression_coefficients(range, dims);
        }

        void precompress_block_commit() noexcept {
            pred_and_quantize_coefficients();
            std::copy(current_coeffs.begin(), current_coeffs.end(), prev_coeffs.begin());
        }

        inline T predict(const iterator &iter) const noexcept {
            T pred = 0;
            for (int i = 0; i < N; i++) {
                pred += iter.get_current_index(i) * current_coeffs[i];
            }
            pred += current_coeffs[N];
            return pred;
        }

        void save(uchar *&c) const {
            std::cout << "save regression predictor" << std::endl;
            c[0] = 0b00000010;
            c += sizeof(uint8_t);
            quantizer.save(c);
            *reinterpret_cast<size_t *>(c) = regression_coeff_quant_inds.size();
            c += sizeof(size_t);
            HuffmanEncoder<int> encoder = HuffmanEncoder<int>();
            encoder.preprocess_encode(regression_coeff_quant_inds, 4 * quantizer.get_radius());
            encoder.save(c);
            encoder.encode(regression_coeff_quant_inds, c);
            encoder.postprocess_encode();
        }

        void predecompress_block(const std::shared_ptr<Range> &range) noexcept {
            pred_and_recover_coefficients();
        }

        void load(const uchar *&c, size_t &remaining_length) {
            //TODO: adjust remaining_length
            std::cout << "load regression predictor" << std::endl;
            c += sizeof(uint8_t);
            remaining_length -= sizeof(uint8_t);
            quantizer.load(c, remaining_length);
            size_t coeff_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            remaining_length -= sizeof(size_t);
            HuffmanEncoder<int> encoder = HuffmanEncoder<int>();
            encoder.load(c, remaining_length);
            regression_coeff_quant_inds = encoder.decode(c, coeff_size);
            encoder.postprocess_decode();
            remaining_length -= coeff_size * sizeof(int);
            std::fill(current_coeffs.begin(), current_coeffs.end(), 0);
            regression_coeff_index = 0;
        }

        void print() const {
            std::cout << "Regression predictor, eb = " << quantizer.get_eb() << "\n";
            int count = 0;
            int ind = regression_coeff_index ? regression_coeff_index : regression_coeff_quant_inds.size();
            std::cout << "\nPrev coeffs: ";
            for (const auto &c:prev_coeffs) {
                std::cout << c << " ";
            }
            std::cout << "\nCurrent coeffs: ";
            for (const auto &c:current_coeffs) {
                std::cout << c << " ";
            }
            std::cout << std::endl;
        }

    private:
        LinearQuantizer<T> quantizer;
        std::vector<int> regression_coeff_quant_inds;
        size_t regression_coeff_index = 0;
        std::array<T, N + 1> current_coeffs;
        std::array<T, N + 1> prev_coeffs;

        std::array<T, N + 1>
        compute_regression_coefficients(const std::shared_ptr<Range> &range, const std::array<size_t, N> &dims) const {
            std::array<double, N + 1> sum{0};
            size_t num_elements = 1;
            for (const auto &dim : dims) {
                num_elements *= dim;
            }
            T num_elements_recip = 1.0 / num_elements;
            {
                auto range_begin = range->begin();
                auto range_end = range->end();
                for (auto iter = range_begin; iter != range_end; ++iter) {
                    T data = *iter;
                    for (int i = 0; i < N; i++) {
                        sum[i] += iter.get_current_index(i) * data;
                    }
                    sum[N] += data;
                }
            }
            std::array<T, N + 1> coeffs;
            coeffs[N] = sum[N] * num_elements_recip;
            for (int i = 0; i < N; i++) {
                coeffs[i] = (2 * sum[i] / (dims[i] - 1) - sum[N]) * 6 * num_elements_recip / (dims[i] + 1);
                coeffs[N] -= (dims[i] - 1) * coeffs[i] / 2;
            }
            return coeffs;
        }

        void pred_and_quantize_coefficients() {
            for (int i = 0; i <= N; i++) {
                regression_coeff_quant_inds.push_back(quantizer.quantize_and_overwrite(current_coeffs[i], prev_coeffs[i]));
            }
        }

        void pred_and_recover_coefficients() {
            for (auto &coeffs : current_coeffs) {
                coeffs = quantizer.recover(coeffs, regression_coeff_quant_inds[regression_coeff_index++]);
            }
        }
    };

}

#endif
