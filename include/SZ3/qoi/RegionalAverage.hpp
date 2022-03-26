//
// Created by Xin Liang on 03/09/2021.
//

#ifndef SZ_QOI_REGIONAL_AVERAGE_HPP
#define SZ_QOI_REGIONAL_AVERAGE_HPP

#include <algorithm>
#include "SZ3/def.hpp"
#include "SZ3/qoi/QoI.hpp"
#include "SZ3/utils/Iterator.hpp"

namespace SZ {
    template<class T, uint N>
    class QoI_RegionalAverage : public concepts::QoIInterface<T, N> {

    public:
        QoI_RegionalAverage(T tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            printf("tolerance = %.4f\n", (double) tolerance);
            printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 3;
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            // T eb = tolerance - fabs(error);
            // T eb = (aggregated_tolerance - fabs(error)) / rest_elements;
            // std::cout << eb << std::endl;
            T eb;
            if(rest_elements > 0.5 * block_elements){
                eb = (aggregated_tolerance - fabs(error)) * 2 / rest_elements;
            }
            else{
                eb = (aggregated_tolerance - fabs(error)) / rest_elements;
            }
            return std::min(eb, global_eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) const {
            return interpret_eb(*data);
        }

        void update_tolerance(T data, T dec_data){
            error += data - dec_data;
            rest_elements --;
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return true;
        }

        void precompress_block(const std::shared_ptr<Range> &range){
            // compute number of elements
            auto dims = range->get_dimensions();
            size_t num_elements = 1;
            for (const auto &dim: dims) {
                num_elements *= dim;
            }
            // assignment
            rest_elements = num_elements;
            block_elements = num_elements;
            aggregated_tolerance = tolerance * num_elements;
        }

        void postcompress_block(){
            error = 0;
        }

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

    private:
        T tolerance;
        T global_eb;
        double error = 0;
        int rest_elements;
        int block_elements;
        double aggregated_tolerance;
    };

    template<class T, uint N>
    class QoI_RegionalAverageInterp : public concepts::QoIInterface<T, N> {

    public:
        QoI_RegionalAverageInterp(T tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            printf("tolerance = %.4f\n", (double) tolerance);
            printf("global_eb = %.4f\n", (double) global_eb);
            tolerance = std::min(tolerance, global_eb);
            concepts::QoIInterface<T, N>::id = 3;
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            return tolerance;
        }

        T interpret_eb(const iterator &iter) const {
            return tolerance;
        }

        T interpret_eb(const T * data, ptrdiff_t offset) const {
            return tolerance;
        }

        void update_tolerance(T data, T dec_data){}

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return true;
        }

        void precompress_block(const std::shared_ptr<Range> &range){}

        void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

    private:
        T tolerance;
        T global_eb;
    };

}
#endif 