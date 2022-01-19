//
// Created by Kai Zhao on 4/28/20.
//

#ifndef SZ_Config_HPP
#define SZ_Config_HPP

#include <vector>
#include <numeric>
#include "SZ3/def.hpp"
#include "MemoryUtil.hpp"

namespace SZ {


    class Config {
    public:
        template<class ... Dims>
        Config(Dims ... args) {
            dims = std::vector<size_t>{static_cast<size_t>(std::forward<Dims>(args))...};
            N = dims.size();
            num = std::accumulate(dims.begin(), dims.end(), (size_t) 1, std::multiplies<>());
            block_size = (N == 1 ? 128 : (N == 2 ? 16 : 6));
            pred_dim = N;
            stride = block_size;
        }

        template<class Iter>
        size_t setDims(Iter begin, Iter end) {
            dims = std::vector<size_t>(begin, end);
            N = dims.size();
            num = std::accumulate(dims.begin(), dims.end(), (size_t) 1, std::multiplies<>());
            return num;
        }


        void save(unsigned char *&c) {
            write(N, c);
            write(dims.data(), dims.size(), c);
            write(num, c);
            write(cmprMethod, c);
            write(errorBoundMode, c);
            write(absErrorBound, c);
            write(relErrorBound, c);
            write(enable_lorenzo, c);
            write(enable_2ndlorenzo, c);
            write(enable_regression, c);
            write(enable_2ndregression, c);
            write(interp_op, c);
            write(interp_direction_op, c);
            write(interp_block_size, c);
            write(lossless_op, c);
            write(encoder_op, c);
            write(quant_state_num, c);
            write(block_size, c);
            write(stride, c);
            write(pred_dim, c);
        };

        void load(const unsigned char *&c) {
            read(N, c);
            dims.resize(N);
            read(dims.data(), N, c);
            read(num, c);
            read(cmprMethod, c);
            read(errorBoundMode, c);
            read(absErrorBound, c);
            read(relErrorBound, c);
            read(enable_lorenzo, c);
            read(enable_2ndlorenzo, c);
            read(enable_regression, c);
            read(enable_2ndregression, c);
            read(interp_op, c);
            read(interp_direction_op, c);
            read(interp_block_size, c);
            read(lossless_op, c);
            read(encoder_op, c);
            read(quant_state_num, c);
            read(block_size, c);
            read(stride, c);
            read(pred_dim, c);
        }

        int N;
        std::vector<size_t> dims;
        size_t num;
        bool enable_lorenzo = true;
        bool enable_2ndlorenzo = false;
        bool enable_regression = true;
        bool enable_2ndregression = false;
        int interp_op = 1;
        int interp_direction_op = 0;
        int interp_block_size = 32;
        int lossless_op = 1; // 0-> skip lossless(use lossless_bypass); 1-> zstd
        int encoder_op = 1;// 0-> skip encoder; 1->HuffmanEncoder; 2->ArithmeticEncoder
        size_t quant_state_num = 65536;
        int block_size, stride, pred_dim;
        double absErrorBound;
        double relErrorBound;
        int errorBoundMode = ABS;
        char cmprMethod = METHOD_LORENZO_REG;

    };


}

#endif //SZ_CONFIG_HPP
