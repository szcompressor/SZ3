//
// Created by Kai Zhao on 4/28/20.
//

#ifndef SZ_CONFIG_HPP
#define SZ_CONFIG_HPP

namespace SZ {
    template<class T, uint N>
    class Config {
    public:
        Config(T _eb, std::array<size_t, N> _dims) : eb(_eb), dims(_dims) {
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

        std::array<size_t, N> dims = {0};
        size_t num;
        bool enable_lorenzo = true;
        bool enable_2ndlorenzo = false;
        bool enable_regression = true;
        bool enable_lossless = true;
        size_t quant_bin = 32768;
        uint block_size, stride, pred_dim = 0;
        T eb;
    };
}

#endif //SZ_CONFIG_HPP
