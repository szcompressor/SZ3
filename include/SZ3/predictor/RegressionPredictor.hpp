#ifndef _SZ_REGRESSION_PREDICTOR_HPP
#define _SZ_REGRESSION_PREDICTOR_HPP

#include "SZ3/def.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include <cassert>

namespace SZ {

    // N-dimension regression predictor
    template<class T, uint N, class Quantizer>
    class RegressionPredictor : public concepts::PredictorInterface<T, N> {
    public:
        using block_iter = typename multi_dimensional_data<T, N>::block_iterator;

        static const uint8_t predictor_id = 0b00000010;

        RegressionPredictor(Quantizer quantizer) : quantizer_independent(0),
                                                   quantizer(quantizer), quantizer_liner(0),
                                                   prev_coeffs{0}, current_coeffs{0} {}

        RegressionPredictor(Quantizer quantizer, int block_size, double eb) : quantizer(quantizer),
                                                                              quantizer_independent(eb / (N + 1)),
                                                                              quantizer_liner(eb / (N + 1) / block_size),
                                                                              prev_coeffs{0}, current_coeffs{0} {
        }

        bool precompress(const block_iter &block) {
            auto range = block.get_block_range();
            auto d = block.mddata;
            auto ds = d->get_dim_strides();

//            auto dims = mddata->get_dims();
            if (N == 3) {
                //preprocess
                T *cur_data_pos = d->get_data(range[0].first, range[1].first, range[2].first);
                double fx = 0.0, fy = 0.0, fz = 0.0, f = 0, sum_x, sum_y;
                int size_x = range[0].second - range[0].first;
                int size_y = range[1].second - range[1].first;
                int size_z = range[2].second - range[2].first;
                if (size_x <= 1 || size_y <= 1 || size_z <= 1) {
                    throw std::invalid_argument("Regression does not support block length in any dimension equal to 1");
                }
                for (int i = 0; i < size_x; i++) {
                    sum_x = 0;
                    for (int j = 0; j < size_y; j++) {
                        sum_y = 0;
                        for (int k = 0; k < size_z; k++) {
                            T curData = *cur_data_pos;
                            sum_y += curData;
                            fz += curData * k;
                            cur_data_pos++;
                        }
                        fy += sum_y * j;
                        sum_x += sum_y;
                        cur_data_pos += (ds[1] - size_z);
                    }
                    fx += sum_x * i;
                    f += sum_x;
                    cur_data_pos += (ds[0] - size_y * ds[1]);
                }
                double coeff = 1.0 / (size_x * size_y * size_z);
                current_coeffs[0] = (2 * fx / (size_x - 1) - f) * 6 * coeff / (size_x + 1);
                current_coeffs[1] = (2 * fy / (size_y - 1) - f) * 6 * coeff / (size_y + 1);
                current_coeffs[2] = (2 * fz / (size_z - 1) - f) * 6 * coeff / (size_z + 1);
                current_coeffs[3] = f * coeff - ((size_x - 1) * current_coeffs[0] / 2 + (size_y - 1) * current_coeffs[1] / 2 +
                                                 (size_z - 1) * current_coeffs[2] / 2);
            }
        }

        bool predecompress(const block_iter &) { return true; }

        void load(const uchar *&c, size_t &remaining_length) {


            c += sizeof(uint8_t);
            remaining_length -= sizeof(uint8_t);

            quantizer.load(c, remaining_length);

            size_t coeff_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            remaining_length -= sizeof(size_t);
            if (coeff_size != 0) {

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

        void save(uchar *&c) const {

            c[0] = 0b00000010;
            c += sizeof(uint8_t);

            quantizer.save(c);

            *reinterpret_cast<size_t *>(c) = regression_coeff_quant_inds.size();
            c += sizeof(size_t);
            if (!regression_coeff_quant_inds.empty()) {
                quantizer_independent.save(c);
                quantizer_liner.save(c);
                HuffmanEncoder<int> encoder = HuffmanEncoder<int>();
                encoder.preprocess_encode(regression_coeff_quant_inds, 0);
                encoder.save(c);
                encoder.encode(regression_coeff_quant_inds, c);
                encoder.postprocess_encode();
            }
        }


        void print() {
//            std::cout << L << "-Layer " << N << "D regression predictor, noise = " << noise << "\n";
        }

        inline T est_error(const block_iter &iter) {
//            return fabs(*iter - predict(iter)) + this->noise;
        }

        size_t size_est() {
            return quantizer.size_est();
        }


        size_t get_padding() {
            return 0;
        }


        inline void compress(const block_iter &block, std::vector<int> &quant_inds) {
            pred_and_quantize_coefficients();
            std::copy(current_coeffs.begin(), current_coeffs.end(), prev_coeffs.begin());

            auto range = block.get_block_range();
            auto d = block.mddata;
            for (size_t i = 0; i < range[0].second - range[0].first; i++) {
                for (size_t j = 0; j < range[1].second - range[1].first; j++) {
                    for (size_t k = 0; k < range[2].second - range[2].first; k++) {
                        T pred = current_coeffs[0] * i + current_coeffs[1] * j + current_coeffs[2] * k + current_coeffs[3];
                        T *c = d->get_data(i + range[0].first, j + range[1].first, k + range[2].first);
                        quant_inds.push_back(quantizer.quantize_and_overwrite(*c, pred));
                    }
                }
            }
        }

        inline void decompress(const block_iter &block, int *&quant_inds_pos) {
            pred_and_recover_coefficients();
            auto range = block.get_block_range();
            auto d = block.mddata;
            for (size_t i = 0; i < range[0].second - range[0].first; i++) {
                for (size_t j = 0; j < range[1].second - range[1].first; j++) {
                    for (size_t k = 0; k < range[2].second - range[2].first; k++) {
                        T pred = current_coeffs[0] * i + current_coeffs[1] * j + current_coeffs[2] * k + current_coeffs[3];
                        T *c = d->get_data(i + range[0].first, j + range[1].first, k + range[2].first);
                        *c = quantizer.recover(pred, *(quant_inds_pos++));
                    }
                }
            }
        }

        void clear() {}

    protected:
        T noise = 0;

    private:
        Quantizer quantizer;
        LinearQuantizer<T> quantizer_liner, quantizer_independent;
        std::vector<int> regression_coeff_quant_inds;
        size_t regression_coeff_index = 0;
        std::array<T, N + 1> current_coeffs;
        std::array<T, N + 1> prev_coeffs;

        void pred_and_quantize_coefficients() {
            for (int i = 0; i < N; i++) {
                regression_coeff_quant_inds.push_back(quantizer_liner.quantize_and_overwrite(current_coeffs[i], prev_coeffs[i]));
            }
            regression_coeff_quant_inds.push_back(
                    quantizer_independent.quantize_and_overwrite(current_coeffs[N], prev_coeffs[N]));
        }

        void pred_and_recover_coefficients() {
            for (int i = 0; i < N; i++) {
                current_coeffs[i] = quantizer_liner.recover(current_coeffs[i],
                                                            regression_coeff_quant_inds[regression_coeff_index++]);
            }
            current_coeffs[N] = quantizer_independent.recover(current_coeffs[N],
                                                              regression_coeff_quant_inds[regression_coeff_index++]);

        }
    };
}
#endif
