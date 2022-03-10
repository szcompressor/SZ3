//
// Created by Xin Liang on 03/09/2022.
//

#ifndef SZ_QOI_LOG_X_HPP
#define SZ_QOI_LOG_X_HPP

#include <cmath>
#include <algorithm>
#include "SZ3/def.hpp"
#include "SZ3/qoi/QoI.hpp"

namespace SZ {
    template<class T>
    class QoI_Log_X : public concepts::QoIInterface<T, 1> {

    public:
        QoI_Log_X(T tolerance, T global_eb, T base=2) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            printf("global_eb = %.4f\n", (double) global_eb);
            printf("tolerance = %.4f\n", (double) tolerance);
            // assuming b > 1
            T coeff1 = fabs(1 - pow(base, - tolerance));
            T coeff2 = fabs(pow(base, tolerance) - 1);
            coeff = std::min(coeff1, coeff2);
            printf("coeff1 = %.4f, coeff2 = %.4f\n", (double) coeff1, (double) coeff2);
            printf("coeff = %.4f\n", (double) coeff);
            log_b = log(base);
            printf("log base = %.4f\n", log_b);
        }

        T interpret_eb(T data) const {
            // if b > 1
            // e = min{(1 - b^{-t})|x|, (b^{t} - 1)|x|}
            // return 0;
            if(data == 0) return 0;
            T eb = coeff * fabs(data);
            return std::min(eb, global_eb);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            if(data == 0) return (dec_data == 0);
            if(verbose){
                std::cout << data << " " << dec_data << std::endl;
                std::cout << log_b_a(fabs(data)) << " " << log_b_a(fabs(dec_data)) << std::endl;
            }
            return (fabs(log_b_a(fabs(data)) - log_b_a(fabs(dec_data))) < tolerance);
        }

        void update_tolerance(T data, T dec_data){}

        void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

    private:
        inline T log_b_a(T a) const {
            return log(a) / log_b;
        }
        T tolerance;
        T global_eb;
        double coeff;
        double log_b;
    };
}
#endif 