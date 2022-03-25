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
            // compute offsets
            // offsets = std::vector<size_t>(dimensions.size());
            // size_t offset = 1;
            // for(int i=dimensions.size()-1; i>=0; i--){
            //     offsets[i] = offset;
            //     offset *= dimensions[i];
            // }
            concepts::QoIInterface<T, N>::id = 4;
            std::sort(isovalues.begin(), isovalues.end());
            std::cout << "isovalues: ";
            for(int i=0; i<isovalues.size(); i++){
                std::cout << isovalues[i] << " ";
            }
            std::cout << "\n";            
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            T eb = global_eb;
            if(isovalues.size() > 3){
                // binary search
                auto iter = std::lower_bound(isovalues.begin(), isovalues.end(), data);
                eb = min(eb, data - *iter);
                if(iter != isovalues.begin()) eb = min(eb, data - *(iter - 1));
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

        T interpret_eb(const T * data, ptrdiff_t offset) const {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return true;
        }

        void update_tolerance(T data, T dec_data){}

        void precompress_block(const std::shared_ptr<Range> &range){}

        void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

    private:
        inline float min(float a, float b) const noexcept{
            return std::min(a, b);
        }

        // template<uint NN = N>
        // inline typename std::enable_if<NN == 1, T>::type do_interpret_eb(const iterator &iter) const noexcept {
        //     std::cerr << "Isoline for 1D is not implemented\n";
        //     exit(-1);
        //     return 0;
        // }

        // template<uint NN = N>
        // inline typename std::enable_if<NN == 2, T>::type do_interpret_eb(const iterator &iter) const noexcept {
        //     // check with isovalue
        //     T x0 = *iter;
        //     T eb = global_eb;
        //     for(int i=0; i<isovalues.size(); i++){
        //         eb = min(eb, fabs(x0 - isovalues[i]));
        //     }
        //     auto index = iter.get_global_index();
        //     // four directions
        //     if(index[1] > 0) eb = min(eb, fabs(x0 - iter.prev(0, 1)));
        //     if(index[0] > 0) eb = min(eb, fabs(x0 - iter.prev(1, 0)));
        //     if(index[1] < dims[1]) eb = min(eb, fabs(x0 - iter.prev(0, -1)));
        //     if(index[0] < dims[0]) eb = min(eb, fabs(x0 - iter.prev(-1, 0)));
        //     return eb;
        // }

        // template<uint NN = N>
        // inline typename std::enable_if<NN == 3, T>::type do_interpret_eb(const iterator &iter) const noexcept {
        //     // check with isovalue
        //     T x0 = *iter;
        //     T eb = global_eb;
        //     for(int i=0; i<isovalues.size(); i++){
        //         eb = min(eb, fabs(x0 - isovalues[i]));
        //     }
        //     auto index = iter.get_global_index();
        //     // six directions
        //     if(index[2] > 0) eb = min(eb, fabs(x0 - iter.prev(0, 0, 1)));
        //     if(index[1] > 0) eb = min(eb, fabs(x0 - iter.prev(0, 1, 0)));
        //     if(index[0] > 0) eb = min(eb, fabs(x0 - iter.prev(1, 0, 0)));
        //     if(index[2] < dims[2]) eb = min(eb, fabs(x0 - iter.prev(0, 0, -1)));
        //     if(index[1] < dims[1]) eb = min(eb, fabs(x0 - iter.prev(0, -1, 0)));
        //     if(index[0] < dims[0]) eb = min(eb, fabs(x0 - iter.prev(-1, 0, 0)));
        //     return eb;
        // }

        // template<uint NN = N>
        // inline typename std::enable_if<NN == 4, T>::type do_interpret_eb(const iterator &iter) const noexcept {
        //     std::cerr << "Isoline for 4D is not implemented\n";
        //     exit(-1);
        //     return 0;
        // }

        // template<uint NN = N>
        // inline typename std::enable_if<NN == 1, T>::type do_interpret_eb_with_ptr(const T * data, ptrdiff_t offset) const noexcept {
        //     std::cerr << "Isoline for 1D is not implemented\n";
        //     exit(-1);
        //     return 0;
        // }

        // template<uint NN = N>
        // inline typename std::enable_if<NN == 2, T>::type do_interpret_eb_with_ptr(const T * data, ptrdiff_t offset) const noexcept {
        //     // check with isovalue
        //     T x0 = *data;
        //     T eb = global_eb;
        //     for(int i=0; i<isovalues.size(); i++){
        //         eb = min(eb, fabs(x0 - isovalues[i]));
        //     }
        //     size_t x = offset / offsets[0];
        //     size_t y = offset % offsets[0];
        //     // four directions
        //     if(y > 0) eb = min(eb, fabs(x0 - *(data - 1)));
        //     if(x > 0) eb = min(eb, fabs(x0 - *(data - offsets[0])));
        //     // check borders
        //     if(y < dims[1]) eb = min(eb, fabs(x0 - *(data + 1)));
        //     if(x < dims[0]) eb = min(eb, fabs(x0 - *(data + offsets[0])));
        //     return eb;
        // }

        // template<uint NN = N>
        // inline typename std::enable_if<NN == 3, T>::type do_interpret_eb_with_ptr(const T * data, ptrdiff_t offset) const noexcept {
        //     // check with isovalue
        //     T x0 = *data;
        //     T eb = global_eb;
        //     for(int i=0; i<isovalues.size(); i++){
        //         eb = min(eb, fabs(x0 - isovalues[i]));
        //     }
        //     size_t x = offset / offsets[0];
        //     offset = offset % offsets[0];
        //     size_t y = offset / offsets[1];
        //     size_t z = offset % offsets[1];
        //     // six directions
        //     if(z > 0) eb = min(eb, fabs(x0 - *(data - 1)));
        //     if(y > 0) eb = min(eb, fabs(x0 - *(data - offsets[1])));
        //     if(x > 0) eb = min(eb, fabs(x0 - *(data - offsets[0])));
        //     // check borders
        //     if(z < dims[2]) eb = min(eb, fabs(x0 - *(data + 1)));
        //     if(y < dims[1]) eb = min(eb, fabs(x0 - *(data + offsets[1])));
        //     if(x < dims[0]) eb = min(eb, fabs(x0 - *(data + offsets[0])));
        //     return eb;
        // }

        template<uint NN = N>
        inline typename std::enable_if<NN == 4, T>::type do_interpret_eb_with_ptr(const T * data, ptrdiff_t offset) const noexcept {
            std::cerr << "Isoline for 4D is not implemented\n";
            exit(-1);
            return 0;
        }

        std::vector<size_t> dims;
        // std::vector<size_t> offsets;
        std::vector<T> isovalues;
        T global_eb;
    };
}
#endif 