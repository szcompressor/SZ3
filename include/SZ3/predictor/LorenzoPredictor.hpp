#ifndef _SZ_LORENZO_PREDICTOR_HPP
#define _SZ_LORENZO_PREDICTOR_HPP

#include "SZ3/def.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include <cassert>

namespace SZ {

    // N-dimension L-layer lorenzo predictor
    template<class T, uint N, uint L, class Quantizer>
    class LorenzoPredictor : public concepts::PredictorInterface<T, N> {
    public:
        using block_iter = typename multi_dimensional_data<T, N>::block_iterator;

        static const uint8_t predictor_id = 0b00000001;

        LorenzoPredictor(Quantizer quantizer) : quantizer(quantizer) {
            this->noise = 0;
        }

        LorenzoPredictor(Quantizer quantizer, double eb) : quantizer(quantizer) {
            this->noise = 0;
            if (L == 1) {
                if (N == 1) {
                    this->noise = 0.5 * eb;
                } else if (N == 2) {
                    this->noise = 0.81 * eb;
                } else if (N == 3) {
                    this->noise = 1.22 * eb;
                } else if (N == 4) {
                    this->noise = 1.79 * eb;
                }
            } else if (L == 2) {
                if (N == 1) {
                    this->noise = 1.08 * eb;
                } else if (N == 2) {
                    this->noise = 2.76 * eb;
                } else if (N == 3) {
                    this->noise = 6.8 * eb;
                }
            }
        }

        inline void save(uchar *&c) const {
            c[0] = predictor_id;
            c += sizeof(uint8_t);
            quantizer.save(c);
        }

        inline void load(const uchar *&c, size_t &remaining_length) {
            c += sizeof(uint8_t);
            remaining_length -= sizeof(uint8_t);
            quantizer.load(c, remaining_length);
        }

        void print() {
            std::cout << L << "-Layer " << N << "D Lorenzo predictor, noise = " << noise << "\n";
        }

#ifndef LORENZO_FUNC
#define prev1(x) (c[ - x])
#define prev2(y, x) (c[ - y * ds[0] - x])
#define prev3(z, y, x) (c[ - z * ds[1] - y * ds[0] - x])
#define prev4(w, z, y, x) (c[ - w * ds[2] - z * ds[1] - y * ds[0] - x])
#endif

        inline T est_error(const block_iter &block) {
            auto range = block.get_block_range();
//            auto d = block.mddata;
            auto ds = block.get_dim_strides();

            size_t min_size = std::numeric_limits<size_t>::max();
            for (const auto &r: range) {
                min_size = std::min(min_size, r.second - r.first);
            }

            T err = 0;
            for (size_t i = 2; i < min_size; i++) {
                size_t bmi = min_size - i;
                if (N == 3 && L == 1) {
                    T *c = block.get_data(i, i, i);
                    T pred = prev3(0, 0, 1) + prev3(0, 1, 0) + prev3(1, 0, 0)
                             - prev3(0, 1, 1) - prev3(1, 0, 1) - prev3(1, 1, 0)
                             + prev3(1, 1, 1);
                    err += fabs(*c - pred) + this->noise;

                    c = block.get_data(i, i, bmi);
                    pred = prev3(0, 0, 1) + prev3(0, 1, 0) + prev3(1, 0, 0)
                           - prev3(0, 1, 1) - prev3(1, 0, 1) - prev3(1, 1, 0)
                           + prev3(1, 1, 1);
                    err += fabs(*c - pred) + this->noise;

                    c = block.get_data(i, bmi, i);
                    pred = prev3(0, 0, 1) + prev3(0, 1, 0) + prev3(1, 0, 0)
                           - prev3(0, 1, 1) - prev3(1, 0, 1) - prev3(1, 1, 0)
                           + prev3(1, 1, 1);
                    err += fabs(*c - pred) + this->noise;

                    c = block.get_data(i, bmi, bmi);
                    pred = prev3(0, 0, 1) + prev3(0, 1, 0) + prev3(1, 0, 0)
                           - prev3(0, 1, 1) - prev3(1, 0, 1) - prev3(1, 1, 0)
                           + prev3(1, 1, 1);
                    err += fabs(*c - pred) + this->noise;
                }
            }
            return err;
        }

        size_t size_est() {
            return quantizer.size_est();
        }


        size_t get_padding() {
            return 2;
        }


        inline void compress(const block_iter &block, std::vector<int> &quant_inds) {
            auto range = block.get_block_range();
            auto d = block.mddata;
            auto ds = d->get_dim_strides();
            if (N == 1 && L == 1) {
                for (size_t i = range[0].first; i < range[0].second; i++) {
                    T *c = d->get_data(i);
                    T pred = prev1(1);
                    quant_inds.push_back(quantizer.quantize_and_overwrite(*c, pred));
                }
            } else if (N == 2 && L == 1) {
                for (size_t i = range[0].first; i < range[0].second; i++) {
                    for (size_t j = range[1].first; j < range[1].second; j++) {
                        for (size_t k = range[2].first; k < range[2].second; k++) {
                            T *c = d->get_data(i, j, k);
                            T pred = prev2(0, 1) + prev2(1, 0) - prev2(1, 1);
                            quant_inds.push_back(quantizer.quantize_and_overwrite(*c, pred));
                        }
                    }
                }
            } else if (N == 3 && L == 1) {
                for (size_t i = range[0].first; i < range[0].second; i++) {
                    for (size_t j = range[1].first; j < range[1].second; j++) {
                        for (size_t k = range[2].first; k < range[2].second; k++) {
                            T *c = d->get_data(i, j, k);
                            T pred = prev3(0, 0, 1) + prev3(0, 1, 0) + prev3(1, 0, 0)
                                     - prev3(0, 1, 1) - prev3(1, 0, 1) - prev3(1, 1, 0)
                                     + prev3(1, 1, 1);
                            quant_inds.push_back(quantizer.quantize_and_overwrite(*c, pred));
                        }
                    }
                }
            } else if (N == 4 && L == 1) {
                for (size_t i = range[0].first; i < range[0].second; i++) {
                    for (size_t j = range[1].first; j < range[1].second; j++) {
                        for (size_t k = range[2].first; k < range[2].second; k++) {
                            for (size_t t = range[3].first; t < range[2].second; t++) {
                                T *c = d->get_data(i, j, k, t);
                                T pred = prev4(0, 0, 0, 1) + prev4(0, 0, 1, 0) - prev4(0, 0, 1, 1) + prev4(0, 1, 0, 0)
                                         - prev4(0, 1, 0, 1) - prev4(0, 1, 1, 0) + prev4(0, 1, 1, 1) + prev4(1, 0, 0, 0)
                                         - prev4(1, 0, 0, 1) - prev4(1, 0, 1, 0) + prev4(1, 0, 1, 1) - prev4(1, 1, 0, 0)
                                         + prev4(1, 1, 0, 1) + prev4(1, 1, 1, 0) - prev4(1, 1, 1, 1);
                                quant_inds.push_back(quantizer.quantize_and_overwrite(*c, pred));
                            }
                        }
                    }
                }
            } else if (N == 1 && L == 2) {
                for (size_t i = range[0].first; i < range[0].second; i++) {
                    T *c = d->get_data(i);
                    T pred = 2 * prev1(1) - prev1(2);
                    quant_inds.push_back(quantizer.quantize_and_overwrite(*c, pred));
                }
            } else if (N == 2 && L == 2) {
                for (size_t i = range[0].first; i < range[0].second; i++) {
                    for (size_t j = range[1].first; j < range[1].second; j++) {
                        T *c = d->get_data(i, j);
                        T pred = 2 * prev2(0, 1) - prev2(0, 2) + 2 * prev2(1, 0)
                                 - 4 * prev2(1, 1) + 2 * prev2(1, 2) - prev2(2, 0)
                                 + 2 * prev2(2, 1) - prev2(2, 2);
                        quant_inds.push_back(quantizer.quantize_and_overwrite(*c, pred));
                    }
                }
            } else if (N == 3 && L == 2) {
                for (size_t i = range[0].first; i < range[0].second; i++) {
                    for (size_t j = range[1].first; j < range[1].second; j++) {
                        for (size_t k = range[2].first; k < range[2].second; k++) {
                            T *c = d->get_data(i, j, k);
                            T pred = 2 * prev3(0, 0, 1) - prev3(0, 0, 2) + 2 * prev3(0, 1, 0)
                                     - 4 * prev3(0, 1, 1) + 2 * prev3(0, 1, 2) - prev3(0, 2, 0)
                                     + 2 * prev3(0, 2, 1) - prev3(0, 2, 2) + 2 * prev3(1, 0, 0)
                                     - 4 * prev3(1, 0, 1) + 2 * prev3(1, 0, 2) - 4 * prev3(1, 1, 0)
                                     + 8 * prev3(1, 1, 1) - 4 * prev3(1, 1, 2) + 2 * prev3(1, 2, 0)
                                     - 4 * prev3(1, 2, 1) + 2 * prev3(1, 2, 2) - prev3(2, 0, 0)
                                     + 2 * prev3(2, 0, 1) - prev3(2, 0, 2) + 2 * prev3(2, 1, 0)
                                     - 4 * prev3(2, 1, 1) + 2 * prev3(2, 1, 2) - prev3(2, 2, 0)
                                     + 2 * prev3(2, 2, 1) - prev3(2, 2, 2);
                            quant_inds.push_back(quantizer.quantize_and_overwrite(*c, pred));
                        }
                    }
                }
            }
        }

        inline void decompress(const block_iter &block, int *&quant_inds_pos) {
            auto range = block.get_block_range();
            auto d = block.mddata;
            auto ds = d->get_dim_strides();
            if (N == 1 && L == 1) {
                for (size_t i = range[0].first; i < range[0].second; i++) {
                    T *c = d->get_data(i);
                    T pred = prev1(1);
                    *c = quantizer.recover(pred, *(quant_inds_pos++));
                }
            } else if (N == 2 && L == 1) {
                for (size_t i = range[0].first; i < range[0].second; i++) {
                    for (size_t j = range[1].first; j < range[1].second; j++) {
                        T *c = d->get_data(i, j);
                        T pred = prev2(0, 1) + prev2(1, 0) - prev2(1, 1);
                        *c = quantizer.recover(pred, *(quant_inds_pos++));
                    }
                }
            } else if (N == 3 && L == 1) {
                for (size_t i = range[0].first; i < range[0].second; i++) {
                    for (size_t j = range[1].first; j < range[1].second; j++) {
                        for (size_t k = range[2].first; k < range[2].second; k++) {
                            T *c = d->get_data(i, j, k);
                            T pred = prev3(0, 0, 1) + prev3(0, 1, 0) + prev3(1, 0, 0)
                                     - prev3(0, 1, 1) - prev3(1, 0, 1) - prev3(1, 1, 0)
                                     + prev3(1, 1, 1);
                            *c = quantizer.recover(pred, *(quant_inds_pos++));
                        }
                    }
                }
            } else if (N == 4 && L == 1) {
                for (size_t i = range[0].first; i < range[0].second; i++) {
                    for (size_t j = range[1].first; j < range[1].second; j++) {
                        for (size_t k = range[2].first; k < range[2].second; k++) {
                            for (size_t t = range[3].first; t < range[2].second; t++) {
                                T *c = d->get_data(i, j, k, t);
                                T pred = prev4(0, 0, 0, 1) + prev4(0, 0, 1, 0) - prev4(0, 0, 1, 1) + prev4(0, 1, 0, 0)
                                         - prev4(0, 1, 0, 1) - prev4(0, 1, 1, 0) + prev4(0, 1, 1, 1) + prev4(1, 0, 0, 0)
                                         - prev4(1, 0, 0, 1) - prev4(1, 0, 1, 0) + prev4(1, 0, 1, 1) - prev4(1, 1, 0, 0)
                                         + prev4(1, 1, 0, 1) + prev4(1, 1, 1, 0) - prev4(1, 1, 1, 1);
                                *c = quantizer.recover(pred, *(quant_inds_pos++));
                            }
                        }
                    }
                }
            } else if (N == 1 && L == 2) {
                for (size_t i = range[0].first; i < range[0].second; i++) {
                    T *c = d->get_data(i);
                    T pred = 2 * prev1(1) - prev1(2);
                    *c = quantizer.recover(pred, *(quant_inds_pos++));
                }
            } else if (N == 2 && L == 2) {
                for (size_t i = range[0].first; i < range[0].second; i++) {
                    for (size_t j = range[1].first; j < range[1].second; j++) {
                        for (size_t k = range[2].first; k < range[2].second; k++) {
                            T *c = d->get_data(i, j, k);
                            T pred = 2 * prev2(0, 1) - prev2(0, 2) + 2 * prev2(1, 0)
                                     - 4 * prev2(1, 1) + 2 * prev2(1, 2) - prev2(2, 0)
                                     + 2 * prev2(2, 1) - prev2(2, 2);
                            *c = quantizer.recover(pred, *(quant_inds_pos++));
                        }
                    }
                }
            } else if (N == 3 && L == 2) {
                for (size_t i = range[0].first; i < range[0].second; i++) {
                    for (size_t j = range[1].first; j < range[1].second; j++) {
                        for (size_t k = range[2].first; k < range[2].second; k++) {
                            T *c = d->get_data(i, j, k);
                            T pred = 2 * prev3(0, 0, 1) - prev3(0, 0, 2) + 2 * prev3(0, 1, 0)
                                     - 4 * prev3(0, 1, 1) + 2 * prev3(0, 1, 2) - prev3(0, 2, 0)
                                     + 2 * prev3(0, 2, 1) - prev3(0, 2, 2) + 2 * prev3(1, 0, 0)
                                     - 4 * prev3(1, 0, 1) + 2 * prev3(1, 0, 2) - 4 * prev3(1, 1, 0)
                                     + 8 * prev3(1, 1, 1) - 4 * prev3(1, 1, 2) + 2 * prev3(1, 2, 0)
                                     - 4 * prev3(1, 2, 1) + 2 * prev3(1, 2, 2) - prev3(2, 0, 0)
                                     + 2 * prev3(2, 0, 1) - prev3(2, 0, 2) + 2 * prev3(2, 1, 0)
                                     - 4 * prev3(2, 1, 1) + 2 * prev3(2, 1, 2) - prev3(2, 2, 0)
                                     + 2 * prev3(2, 2, 1) - prev3(2, 2, 2);
                            *c = quantizer.recover(pred, *(quant_inds_pos++));
                        }
                    }
                }
            }
        }

        void clear() {}

    protected:
        T noise = 0;

    private:
        Quantizer quantizer;
    };
}
#endif
