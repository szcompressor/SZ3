#ifndef _SZ_INTEGER_QUANTIZER2_HPP
#define _SZ_INTEGER_QUANTIZER2_HPP

#include <cstring>
#include <cassert>
#include <iostream>
#include <vector>
#include "def.hpp"
#include "quantizer/Quantizer.hpp"

namespace SZ {

    template<class T>
    class LinearQuantizer2 : public concepts::QuantizerInterface<T> {
    public:
        LinearQuantizer2(size_t num_) : num(num_), error_bound(1), error_bound_reciprocal(1), radius(32768) {}

        LinearQuantizer2(size_t num_, T eb, int r = 32768) : num(num_), error_bound(eb),
                                                             error_bound_reciprocal(1.0 / eb),
                                                             radius(r) {
            assert(eb != 0);
        }

        int get_radius() const { return radius; }

        T get_eb() const { return error_bound; }

        void set_eb(T eb) {
            error_bound = eb;
            error_bound_reciprocal = 1.0 / eb;
        }


        T recover(T pred, int quant_index) {
            printf("IntegerQuantizer2.recover(T pred, int quant_index) is not supported");
            exit(0);
        }

        int quantize_and_overwrite(T &data, T pred) {
            printf("IntegerQuantizer2.quantize_and_overwrite(T &data, T pred) is not supported");
            exit(0);
        }

        // quantize the data with a prediction value, and returns the quantization index and the decompressed data
        // int quantize(T data, T pred, T& dec_data);
        int quantize_and_overwrite(size_t idx, T &data, T pred) {
            T diff = data - pred;
            int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal) + 1;
            if (quant_index < this->radius * 2) {
                quant_index >>= 1;
                int half_index = quant_index;
                quant_index <<= 1;
                int quant_index_shifted;
                if (diff < 0) {
                    quant_index = -quant_index;
                    quant_index_shifted = -half_index;
                } else {
                    quant_index_shifted = half_index;
                }
                T decompressed_data = pred + quant_index * this->error_bound;
                if (fabs(decompressed_data - data) > this->error_bound) {
                    unpred_idx.push_back(idx);
                    unpred_val.push_back(data);
                    return 0;
                } else {
                    data = decompressed_data;
                    return quant_index_shifted;
                }
            } else {
                unpred_idx.push_back(idx);
                unpred_val.push_back(data);
                return 0;
            }
        }


        void recover_and_residual(size_t idx, T &dest, T &delta, T pred, int quant_index) {
            if (!isunpred[idx]) {
                T temp = recover_pred(pred, quant_index);
                dest += temp;
                delta += temp;
            } else {
                dest = unpred_map[idx];
                delta = 0;
            }
        }

        void recover_and_residual(size_t idx, T &dest, T &delta, T pred) {
            if (!isunpred[idx]) {
                dest += pred;
                delta += pred;
            } else {
                dest = unpred_map[idx];
                delta = 0;
            }
        };

        // recover the data using the quantization index
        T recover(size_t idx, T pred, int quant_index) {
            if (!isunpred[idx]) {
                return recover_pred(pred, quant_index);
            } else {
                return unpred_map[idx];
            }
        }


        T recover_pred(T pred, int quant_index) {
            return pred + 2 * (quant_index) * this->error_bound;
        }

        void save(unsigned char *&c) const {
            // std::string serialized(sizeof(uint8_t) + sizeof(T) + sizeof(int),0);
            c[0] = 0b00000010;
            c += 1;
            // std::cout << "saving eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
            *reinterpret_cast<T *>(c) = this->error_bound;
            c += sizeof(T);
            *reinterpret_cast<int *>(c) = this->radius;
            c += sizeof(int);
            size_t unpred_size = unpred_idx.size();
            *reinterpret_cast<size_t *>(c) = unpred_size;
            c += sizeof(size_t);
            memcpy(c, unpred_idx.data(), unpred_size * sizeof(size_t));
            c += unpred_size * sizeof(size_t);
            memcpy(c, unpred_val.data(), unpred_size * sizeof(T));
            c += unpred_size * sizeof(T);
//            printf("xx %lu\n", unpred_idx.size());

        };

        void load(const unsigned char *&c, size_t &remaining_length) {
            assert(remaining_length > (sizeof(uint8_t) + sizeof(T) + sizeof(int)));
            auto c0 = c;
            c += sizeof(uint8_t);
//            remaining_length -= sizeof(uint8_t);
            this->error_bound = *reinterpret_cast<const T *>(c);
            this->error_bound_reciprocal = 1.0 / this->error_bound;
            c += sizeof(T);
            this->radius = *reinterpret_cast<const int *>(c);
            c += sizeof(int);
            size_t unpred_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            const size_t *unpred_idx_ptr = reinterpret_cast<const size_t *>(c);
            c += unpred_size * sizeof(size_t);
            const T *unpred_val_ptr = reinterpret_cast<const T *>(c);
            c += unpred_size * sizeof(T);
            isunpred.resize(num, false);
            for (size_t i = 0; i < unpred_size; i++) {
                size_t idx = unpred_idx_ptr[i];
                unpred_map[idx] = unpred_val_ptr[i];
                isunpred[idx] = true;
            }
            // std::cout << "loading: eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
            // reset index
//            printf("xx %lu\n", unpred_map.size());
            remaining_length -= (c - c0);

        }

        void clear() {
        }

        size_t get_unpred_size() {
            return unpred_idx.size();
        }

        virtual void postcompress_data() {};

        virtual void postdecompress_data() {};

        virtual void precompress_data() {};

        virtual void predecompress_data() {};

    private:
        T error_bound;
        T error_bound_reciprocal;
        int radius; // quantization interval radius
        size_t num;

        std::vector<T> unpred_val; // for compression
        std::vector<size_t> unpred_idx; // for compression
        ska::unordered_map<size_t, T> unpred_map; // for decompression
        std::vector<bool> isunpred; //for decompression

    };

}
#endif
