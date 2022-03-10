//
// Created by Xin Liang on 03/09/2021.
//

#ifndef SZ_QOI_REGIONAL_AVERAGE_HPP
#define SZ_QOI_REGIONAL_AVERAGE_HPP

#include <algorithm>
#include "SZ3/def.hpp"
#include "SZ3/qoi/QoI.hpp"

namespace SZ {
    template<class T>
    class QoI_RegionalAverage : public concepts::QoIInterface<T, 1> {

    public:
        QoI_RegionalAverage(T tolerance, T global_eb, int block_elements=1) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            printf("tolerance = %.4f\n", (double) tolerance);
            printf("global_eb = %.4f\n", (double) global_eb);
        }

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

        void update_tolerance(T data, T dec_data){
            error += data - dec_data;
            rest_elements --;
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return true;
        }

        void precompress_block(int num_elements){
            rest_elements = num_elements;
            block_elements = num_elements;
            aggregated_tolerance = tolerance * num_elements;
        }

        void postcompress_block(){
            error = 0;
        }

        void print(){}

    private:
        T tolerance;
        T global_eb;
        double error = 0;
        int rest_elements;
        int block_elements;
        double aggregated_tolerance;
    };
}
#endif 