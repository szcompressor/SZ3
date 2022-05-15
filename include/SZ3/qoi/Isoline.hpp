//
// Created by Xin Liang on 03/23/2022.
//

#ifndef SZ_QOI_ISOLINE_HPP
#define SZ_QOI_ISOLINE_HPP

#include <algorithm>
#include "SZ3/def.hpp"
#include "SZ3/qoi/QoI.hpp"
#include "SZ3/utils/Iterator.hpp"

namespace SZ {
    template<class T, uint N>
    class QoI_Isoline : public concepts::QoIInterface<T, N> {

    public:
        QoI_Isoline(std::vector<size_t> dimensions, std::vector<T> values, T global_eb) : 
                dims(dimensions),
                isovalues(values),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 4;
            std::sort(isovalues.begin(), isovalues.end());
            std::cout << "isovalues: ";
            for(int i=0; i<isovalues.size(); i++){
                // std::cout << isovalues[i] << " ";
                printf("%.20f ", (float) isovalues[i]);
            }
            std::cout << "\n";            
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            T eb = global_eb;
            if(isovalues.size() > 5){
                // binary search
                auto iter = std::lower_bound(isovalues.begin(), isovalues.end(), data);
                eb = min(eb, fabs(data - *iter));
                if(iter != isovalues.begin()) eb = min(eb, fabs(data - *(iter - 1)));
            }
            else{
                for(int i=0; i<isovalues.size(); i++){
                    eb = min(eb, fabs(data - isovalues[i]));
                }                
            }
            return eb;
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            if(isovalues.size() > 5){
                // binary search
                auto iter = std::lower_bound(isovalues.begin(), isovalues.end(), data);
                if((data - *iter)*(dec_data - *iter) < 0) return false;
                if(iter != isovalues.begin()){
                    if((data - *(iter - 1))*(dec_data - *(iter - 1)) < 0) return false;
                }
            }
            else{            
                for(int i=0; i<isovalues.size(); i++){
                    if((data - isovalues[i])*(dec_data - isovalues[i]) < 0) return false;
                }
            }
            return true;
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
        inline float min(float a, float b) const noexcept{
            return std::min(a, b);
        }

        std::vector<size_t> dims;
        std::vector<T> isovalues;
        T global_eb;
    };
}
#endif 