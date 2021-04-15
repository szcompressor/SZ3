//
// Created by Kai Zhao on 4/28/20.
//

#ifndef SZ_CONFIG_HPP
#define SZ_CONFIG_HPP

#include <array>

namespace SZ {
    template<class T, uint N>
    class Config {
    public:

        Config(std::array<size_t, N> _dims) {
            init(_dims);
        }

        Config(T _eb, std::array<size_t, N> _dims) : eb(_eb) {
            init(_dims);
        }

        std::array<size_t, N> dims = {0};
        size_t num;
        bool enable_lorenzo = true;
        bool enable_2ndlorenzo = false;
        bool enable_regression = true;
        int lossless_op = 1; // 0-> skip lossless(use lossless_bypass); 1-> zstd
        int encoder_op = 1;// 0-> skip encoder(use PQLCompressor); 1->HuffmanEncoder; 2->ArithmeticEncoder
        size_t quant_state_num = 65536;
        uint block_size, stride, pred_dim = 0;
        T eb = 0;
        T relative_eb = 0;
        int eb_mode = 0;
        std::string src_file_name;
        size_t timestep_batch = 0;

    private:

        void init(std::array<size_t, N> _dims) {
            dims = _dims;
            switch (N) {
                case 1:
                    block_size = 128;
                    break;
                case 2:
                    block_size = 16;
                    break;
                default:
                    // >= 3D
                    block_size = 6;
                    break;
            }
            stride = block_size;
            num = 1;
            for (const auto &d:_dims) {
                num *= d;
            }
        }
    };
}

#endif //SZ_CONFIG_HPP
