#ifndef _SZ_POLY_REGRESSION_PREDICTOR_HPP
#define _SZ_POLY_REGRESSION_PREDICTOR_HPP

#include "def.hpp"
#include "utils/Iterator.hpp"
#include "utils/Compat.hpp"
#include "utils/fileUtil.h"
#include "quantizer/Quantizer.hpp"
#include "encoder/Encoder.hpp"
#include <cstring>
#include <iostream>

namespace SZ {

#define COEF_N 10
#define COEF_AUX_DIM 10
#define COEF_AUX_DIM_START 3

    // N-d regression predictor
    template<class T, uint N>
    class PolyRegressionPredictor {
    public:
        static const uint8_t predictor_id = 0b00000011;

        PolyRegressionPredictor() : quantizer(0), current_coeffs{0} {
            init_poly();
        }

        PolyRegressionPredictor(T eb) : quantizer(eb), current_coeffs{0} {
            init_poly();
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        void precompress_data(const iterator &) const noexcept {}

        void precompress_block_commit() noexcept {}

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
            std::array<T, COEF_N> coeffs = compute_regression_coefficients(range, dims);
            pred_and_quantize_coefficients(coeffs);
            std::copy(coeffs.begin(), coeffs.end(), current_coeffs.begin());
        }

        inline std::array<T, COEF_N> get_poly_index(const iterator &iter) const {
            T i = iter.get_current_index(0);
            T j = iter.get_current_index(1);
            T k = iter.get_current_index(2);

            return std::array<T, COEF_N>{1.0, i, j, k, i * i, i * j, i * k, j * j, j * k, k * k};
        }

        inline T predict(const iterator &iter) const noexcept {
            T pred = 0;
            auto poly_index = get_poly_index(iter);
            for (int i = 0; i < COEF_N; i++) {
                pred += poly_index[i] * current_coeffs[i];
            }
            return pred;
        }

        void save(uchar *&c) const {
            std::cout << "save predictor" << std::endl;
            c[0] = predictor_id;
            c += 1;
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
            std::cout << "load predictor" << std::endl;
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
        }

        void display_coef_aux(std::array<T, COEF_N * COEF_N> aux) {
            for (int i = 0; i < COEF_N; i++) {
                for (int j = 0; j < COEF_N; j++)
                    std::cout << std::setw(10) << std::setprecision(6) << aux[i * COEF_N + j] << " ";
                std::cout << std::endl;
            }
        }

        void init_poly() {
            size_t num = 0;
            std::string bin_file_path = SZ::get_cmake_project_source_dir() + "/include/predictor/PolyRegressionCoefAux.f32";
            auto data = SZ::readfile<float>(bin_file_path.data(), num);

            auto data1 = &data[0];
            coef_aux = std::vector<std::array<T, COEF_N * COEF_N>>((COEF_AUX_DIM + 1) * (COEF_AUX_DIM + 1) * (COEF_AUX_DIM + 1),
                                                                   {0});
            for (int i = COEF_AUX_DIM_START; i <= COEF_AUX_DIM; i++) {
                for (int j = COEF_AUX_DIM_START; j <= COEF_AUX_DIM; j++) {
                    for (int k = COEF_AUX_DIM_START; k <= COEF_AUX_DIM; k++) {
                        std::copy_n(data1, COEF_N * COEF_N,
                                    coef_aux[i * COEF_AUX_DIM * COEF_AUX_DIM + j * COEF_AUX_DIM + k].begin());
                        data1 += COEF_N * COEF_N;
                    }
                }
            }
//            display_coef_aux(coef_aux[6 * COEF_AUX_DIM * COEF_AUX_DIM + 6 * COEF_AUX_DIM + 6]);
        }

    private:
        LinearQuantizer<T> quantizer;
        std::vector<int> regression_coeff_quant_inds;
        size_t regression_coeff_index = 0;
        std::array<T, COEF_N> current_coeffs;
        std::vector<std::array<T, COEF_N * COEF_N>> coef_aux;

        std::array<T, COEF_N>
        compute_regression_coefficients(const std::shared_ptr<Range> &range, const std::array<size_t, N> &dims) const {
            std::array<double, COEF_N> sum{0};
            {
                for (auto iter = range->begin(); iter != range->end(); ++iter) {
                    T data = *iter;
                    auto poly_index = get_poly_index(iter);
                    for (int i = 0; i < COEF_N; i++) {
                        sum[i] += poly_index[i] * data;
                    }
                }
            }
            std::array<double, COEF_N> coeffs{0};
            auto coef_aux_index = 0;
            for (auto &dim:dims) {
                coef_aux_index = coef_aux_index * COEF_AUX_DIM + dim;
            }
            auto coef_aux1 = coef_aux[coef_aux_index];

            for (int i = 0; i < COEF_N; i++) {
                for (int j = 0; j < COEF_N; j++) {
                    coeffs[i] += coef_aux1[i * COEF_N + j] * sum[j];
                }
            }
//            std::cout << "sum ";
//            for (auto &su:sum) {
//                std::cout << su << " ";
//            }
//            std::cout << std::endl;
//
//            std::cout << "coeff ";
//            for (auto &coef:coeffs) {
//                std::cout << coef << " ";
//            }
//            std::cout << std::endl;

            std::array<T, COEF_N> coeffsT;
            for (int i = 0; i < COEF_N; i++) {
                coeffsT[i] = coeffs[i];
            }
            return coeffsT;
        }

        void pred_and_quantize_coefficients(std::array<T, COEF_N> &coeffs) {
            for (int i = 0; i <= N; i++) {
                regression_coeff_quant_inds.push_back(quantizer.quantize_and_overwrite(coeffs[i], current_coeffs[i]));
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