#ifndef SZ3_QOI_INTERFACE
#define SZ3_QOI_INTERFACE
/**
 * Interface for some specific quantities of interest (QoIs)
 * Created by Xin Liang on 12/06/2021.
 */

#include "SZ3/utils/Iterator.hpp"

namespace SZ {


    namespace concepts {

        template<class T, uint N>
        class QoIInterface {
        public:

            using Range = multi_dimensional_range<T, N>;
            using iterator = typename multi_dimensional_range<T, N>::iterator;

            virtual ~QoIInterface() = default;

            virtual T interpret_eb(T data) const = 0;

            // interpret eb with iterator (Lorenzo)
            virtual T interpret_eb(const iterator &iter) const = 0;

            // interpret eb with data pointer (Interpolation)
            virtual T interpret_eb(const T * data, ptrdiff_t offset) = 0;

            virtual void update_tolerance(T data, T dec_data) = 0;

            virtual bool check_compliance(T data, T dec_data, bool verbose=false) const = 0;

            virtual void precompress_block(const std::shared_ptr<Range> &range) = 0;

            virtual void postcompress_block() = 0;

            virtual void print() = 0;

            // for interpolation compressors
            virtual T get_global_eb() const = 0;

            virtual void set_global_eb(T eb) = 0;

            virtual void init() = 0;

            virtual void set_dims(const std::vector<size_t>& new_dims) = 0;

            int id = 0;
        };

    }

}

#endif
