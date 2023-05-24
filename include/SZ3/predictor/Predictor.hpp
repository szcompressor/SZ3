#ifndef _SZ_PREDICTOR_HPP
#define _SZ_PREDICTOR_HPP

#include "SZ3/utils/MDdata.hpp"
#include "SZ3/def.hpp"

namespace SZ {


    namespace concepts {

        template<class T, uint N>
        class PredictorInterface {
        public:
            using block_iter = typename multi_dimensional_data<T, N>::block_iterator;

            virtual ~PredictorInterface() = default;

            virtual inline bool precompress(const block_iter &) = 0;

            virtual inline void compress(const block_iter &, std::vector<int> &) = 0;

            virtual inline bool predecompress(const block_iter &) = 0;

            virtual inline void decompress(const block_iter &, int *&) = 0;

            virtual inline void save(uchar *&c) const = 0;

            virtual inline void load(const uchar *&c, size_t &remaining_length)  = 0;

            virtual inline T est_error(const block_iter &) = 0;

            virtual inline size_t size_est() = 0;

            virtual inline void print() = 0;

            virtual inline void clear() = 0;
        };

    }

}

#endif
