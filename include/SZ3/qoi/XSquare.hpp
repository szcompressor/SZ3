//
// Created by Xin Liang on 12/06/2021.
//

#ifndef SZ_QOI_X_SQUARE_HPP
#define SZ_QOI_X_SQUARE_HPP

#include <algorithm>
#include "SZ3/def.hpp"
#include "SZ3/qoi/QoI.hpp"
#include "SZ3/utils/Iterator.hpp"

namespace SZ {
    template<class T, uint N>
    class QoI_X_Square : public concepts::QoIInterface<T, N> {

    public:
        QoI_X_Square(T tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 1;
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            // Implementation 1
            // e^2 + 2ex - t < 0
            // T x1 = - data + sqrt(data * data + tolerance);
            // T x2 = data + sqrt(data * data + tolerance);
            // T eb = std::min(x1, x2);
            // e^2 + 2ex + t > 0
            // if ((data < 0) && (data * data - tolerance >= 0)){
            //     if (eb > - data - sqrt(data * data - tolerance))
            //         eb = - data - sqrt(data * data - tolerance);
            // }
            // Implementation 2
            // T eb;
            // if(data >= 0){
            //     eb = - data + sqrt(data * data + tolerance);
            // }
            // else{
            //     eb = data + sqrt(data * data + tolerance);
            //     if(data * data - tolerance >= 0){
            //         T eb2 = - data - sqrt(data * data - tolerance);
            //         eb = std::min(eb, eb2);
            //     }
            // }
            // Implementation 3
            T eb = - fabs(data) + sqrt(data * data + tolerance);
            return std::min(eb, global_eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return (fabs(data*data - dec_data*dec_data) < tolerance);
        }

        void update_tolerance(T data, T dec_data){}

        void precompress_block(const std::shared_ptr<Range> &range){}

        void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

    private:
        T tolerance;
        T global_eb;
    };
}
#endif 
