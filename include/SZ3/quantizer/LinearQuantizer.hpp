    #ifndef SZ3_LINEAR_QUANTIZER_HPP
    #define SZ3_LINEAR_QUANTIZER_HPP

    #include <cassert>
    #include <cstring>
    #include <iostream>
    #include <vector>

    #include "SZ3/def.hpp"
    #include "SZ3/utils/MemoryUtil.hpp"
    #include "SZ3/quantizer/Quantizer.hpp"

    namespace SZ3 {

    template <class T>
    class LinearQuantizer : public concepts::QuantizerInterface<T, int> {
       public:
        LinearQuantizer() : error_bound(1), error_bound_reciprocal(1), radius(32768) {}

        LinearQuantizer(double eb, int r = 32768) : error_bound(eb), error_bound_reciprocal(1.0 / eb), radius(r) {
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
        ALWAYS_INLINE int quantize_and_overwrite(T &data, T pred) override {
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
                if (fabs(decompressed_data - data) > this->error_bound) {
                    unpred.push_back(data);
                    return 0;
                } else {
                    data = decompressed_data;
                    return quant_index_shifted;
                }
            } else {
                unpred.push_back(data);
                return 0;
            }
        }

        /**
         * For metaLorenzo only, will be removed together with metalorenzo
         * @param ori
         * @param pred
         * @param dest
         * @return
         */
        ALWAYS_INLINE int quantize_and_overwrite(T ori, T pred, T &dest) {
            T diff = ori - pred;
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
                if (fabs(decompressed_data - ori) > this->error_bound) {
                    unpred.push_back(ori);
                    dest = ori;
                    return 0;
                } else {
                    dest = decompressed_data;
                    return quant_index_shifted;
                }
            } else {
                unpred.push_back(ori);
                dest = ori;
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

        ALWAYS_INLINE int force_save_unpred(T ori) override{
            unpred.push_back(ori);
            return 0;
        }

        size_t size_est() { return unpred.size() * sizeof(T); }

        void save(unsigned char *&c) const override {
            write(uid, c);
            write(this->error_bound, c);
            write(this->radius, c);
            size_t unpred_size = unpred.size();
            write(unpred_size, c);
            write(unpred.data(), unpred.size(), c);
        }

        void load(const unsigned char *&c, size_t &remaining_length) override {
            uchar uid_read;
            read(uid_read, c, remaining_length);
            if (uid_read != uid) {
                fprintf(stderr, "LinearQuantizer uid mismatch\n");
                throw std::invalid_argument("LinearQuantizer uid mismatch");
            }
            read(this->error_bound, c, remaining_length);
            this->error_bound_reciprocal = 1.0 / this->error_bound;
            read(this->radius, c, remaining_length);
            size_t unpred_size = 0;
            read(unpred_size, c, remaining_length);
            unpred.resize(unpred_size);
            read(unpred.data(), unpred_size, c, remaining_length);

            index = 0;
        }

        void print() override {
            printf("[LinearQuantizer] error_bound = %.8G, radius = %d, unpred = %lu\n", error_bound, radius, unpred.size());
        }

       private:
        std::vector<T> unpred;
        size_t index = 0;  // used in decompression only
        uchar uid = 0b10;

        double error_bound;
        double error_bound_reciprocal;
        int radius;  // quantization interval radius
    };

    }  // namespace SZ3
    #endif