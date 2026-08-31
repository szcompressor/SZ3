/**
 * @file PreFilter.hpp
 * @ingroup Preprocessor
 */

#ifndef SZ3_PREFILTER_HPP
#define SZ3_PREFILTER_HPP

#include "SZ3/preprocessor/PreProcessor.hpp"

namespace SZ3 {

    template<class T, uint N>

    class PreFilter : public concepts::PreprocessorInterface<T, N> {

        void preprocess(T *data, std::array<size_t, N> dims, std::pair<T, T> range, T defaultValue) {
            size_t num = 1;
            for (size_t i = 0; i < N; i++) {
                num *= dims[i];
            }
            for (size_t i = 0; i < num; i++) {
                if (data[i] > range.second || data[i] < range.first) {
                    data[i] = defaultValue;
                }
            }
        }
    };
}
#endif //SZ3_PRETRANSPOSE_H
