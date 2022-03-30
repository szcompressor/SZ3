#ifndef _SZ_QOI_INTEGER_QUANTIZER_HPP
#define _SZ_QOI_INTEGER_QUANTIZER_HPP

#include <cstring>
#include <cassert>
#include <algorithm>
#include <vector>
#include "SZ3/quantizer/Quantizer.hpp"

#include <iostream>

namespace SZ {

    template<class T, class T_eb>
    class VariableEBLinearQuantizer : public concepts::QuantizerInterface<T> {
    public:
        VariableEBLinearQuantizer(int r = 32768) : radius(r) {}

        int get_radius() const { return radius; }

        // quantize the data with a prediction value, and returns the quantization index
        int quantize(T data, T pred, T_eb eb) {
            if(eb == 0) return 0;
            T diff = data - pred;
            int quant_index = (int) (fabs(diff) / eb) + 1;
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
                T decompressed_data = pred + quant_index * eb;
                if (fabs(decompressed_data - data) > eb) {
                    return 0;
                } else {
                    return quant_index_shifted;
                }
            } else {
                return 0;
            }
        }

        // quantize the data with a prediction value, and returns the quantization index and the decompressed data
        int quantize_and_overwrite(T &data, T pred, T_eb eb) {
            if(eb == 0){
                unpred.push_back(data);
                return 0;
            }
            // if(fabs(data + 14.3927) < 0.0001){
            //     std::cout << data << " " << pred << " " << eb << std::endl;
            // }
            T diff = data - pred;
            int quant_index = (int) (fabs(diff) / eb) + 1;
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
                T decompressed_data = pred + quant_index * eb;
                // if(fabs(data + 14.3927) < 0.0001){
                //     std::cout << decompressed_data << ", err = " << decompressed_data - data << std::endl;
                // }
                if (fabs(decompressed_data - data) > eb) {
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

        int quantize_and_overwrite(T ori, T pred, T_eb eb, T &dest) {
            if(eb == 0){
                unpred.push_back(ori);
                dest = ori;
                return 0;
            }
            T diff = ori - pred;
            int quant_index = (int) (fabs(diff) / eb) + 1;
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
                T decompressed_data = pred + quant_index * eb;
                if (fabs(decompressed_data - ori) > eb) {
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
        T recover(T pred, int quant_index, T_eb eb) {
            if (quant_index) {
                return recover_pred(pred, quant_index, eb);
            } else {
                return recover_unpred();
            }
        }

        T recover_pred(T pred, int quant_index, T_eb eb) {
            return pred + 2 * (quant_index - this->radius) * eb;
        }

        T recover_unpred() {
            return unpred[index++];
        }

        // required function in Quantizer interface
        int quantize(T data, T pred) {
            return 0;
        }

        // required function in Quantizer interface
        int quantize_and_overwrite(T &data, T pred) {
            return 0;
        }

        // required function in Quantizer interface
        T recover(T pred, int quant_index){
            return 0;
        }

        size_t size_est() {
            return unpred.size() * sizeof(T);
        }

        void save(unsigned char *&c) const {
            // std::string serialized(sizeof(uint8_t) + sizeof(T) + sizeof(int),0);
            c[0] = 0b00000010;
            c += 1;
            // std::cout << "saving eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
            *reinterpret_cast<int *>(c) = this->radius;
            c += sizeof(int);
            *reinterpret_cast<size_t *>(c) = unpred.size();
            c += sizeof(size_t);
            memcpy(c, unpred.data(), unpred.size() * sizeof(T));
            c += unpred.size() * sizeof(T);
            // std::cout << "unpred size = " << unpred.size() << std::endl;
        };

        void load(const unsigned char *&c, size_t &remaining_length) {
            assert(remaining_length > (sizeof(uint8_t) + sizeof(T) + sizeof(int)));
            c += sizeof(uint8_t);
            this->radius = *reinterpret_cast<const int *>(c);
            c += sizeof(int);
            size_t unpred_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            this->unpred = std::vector<T>(reinterpret_cast<const T *>(c), reinterpret_cast<const T *>(c) + unpred_size);
            c += unpred_size * sizeof(T);
            // reset index
            index = 0;
        }

        void clear() {
            // std::cout << "unpred size = " << unpred.size() << std::endl;
            unpred.clear();
            index = 0;
        }

        virtual void postcompress_data() {}

        virtual void postdecompress_data() {}

        virtual void precompress_data() {}

        virtual void predecompress_data() {}


    private:
        std::vector<T> unpred;
        size_t index = 0; // used in decompression only
        int radius; // quantization interval radius
    };

    template<class T>
    class EBLogQuantizer : public concepts::QuantizerInterface<T> {
    public:
        EBLogQuantizer(T eb_base = std::numeric_limits<T>::epsilon(), T log_base = 2, int r = 32) : 
                eb_base(eb_base), log_base(log_base), radius(r) {
                    // TODO: adjust for int data
                    set_reciprocal();
                    //printf("eb_base = %.4f, log_base = %.4f\n", (double) eb_base, (double) log_base);
                }

        int get_radius() const { return radius; }

        // quantize the data with a prediction value, and returns the quantization index
        int quantize(T eb) {
            if(eb <= eb_base){
                eb = 0;
                return 0;
            }
            int id = log2(eb * eb_base_reciprocal) * log_of_base_reciprocal;
            return id;
        }

        // quantize the error bound, and returns the quantization index and the decompressed data
        int quantize_and_overwrite(T &eb) {
            // std::cout << eb << " ";
            if(eb <= eb_base){
                eb = 0;
                return 0;
            }
            int id = log2(eb * eb_base_reciprocal) * log_of_base_reciprocal;
            // need to check if id = 0
            if(id == 0){
                eb = 0;
                return 0;
            }
            id = std::min(id, radius);
            eb = pow(log_base, id) * eb_base;
            return id;            
        }

        // recover the error bound using the quantization index
        T recover(int quant_index) {
            if (quant_index == 0) {
                return 0;
            } else {
                return pow(log_base, quant_index) * eb_base;
            }
        }

        // required function in Quantizer interface
        T recover(T pred, int quant_index) {
            return recover(quant_index);
        }

        // required function in Quantizer interface
        int quantize(T data, T pred) {
            return quantize(data);
        }

        // required function in Quantizer interface
        int quantize_and_overwrite(T &data, T pred) {
            quantize_and_overwrite(data);
            return 0;
        }

        void save(unsigned char *&c) const {
            // std::string serialized(sizeof(uint8_t) + sizeof(T) + sizeof(int),0);
            c[0] = 0b00000010;
            c += 1;
            // std::cout << "saving eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
            *reinterpret_cast<T *>(c) = eb_base;
            c += sizeof(T);
            *reinterpret_cast<T *>(c) = log_base;
            c += sizeof(T);
            *reinterpret_cast<int *>(c) = this->radius;
            c += sizeof(int);
        };

        void load(const unsigned char *&c, size_t &remaining_length) {
            assert(remaining_length > (sizeof(uint8_t) + sizeof(T) + sizeof(int)));
            c += sizeof(uint8_t);
            eb_base = *reinterpret_cast<const T *>(c);            
            c += sizeof(T);
            log_base = *reinterpret_cast<const T *>(c);
            c += sizeof(T);
            set_reciprocal();
            this->radius = *reinterpret_cast<const int *>(c);
            c += sizeof(int);
        }

        void clear() {}

        virtual void postcompress_data() {}

        virtual void postdecompress_data() {}

        virtual void precompress_data() {}

        virtual void predecompress_data() {}

    private:
        void set_reciprocal(){
            eb_base_reciprocal = ((T) 1.0) / eb_base;
            log_of_base_reciprocal = ((T) 1.0) / log2(log_base);
        }

        int radius; // quantization interval radius
        T eb_base;
        T log_base;
        T eb_base_reciprocal;
        T log_of_base_reciprocal;
    };

}
#endif
