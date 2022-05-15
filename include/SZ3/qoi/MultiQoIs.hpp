//
// Created by Xin Liang on 03/26/2022.
//

#ifndef SZ_QOI_MULTI_HPP
#define SZ_QOI_MULTI_HPP

#include <algorithm>
#include "SZ3/def.hpp"
#include "SZ3/qoi/QoI.hpp"
#include "SZ3/utils/Iterator.hpp"

namespace SZ {
    template<class T, uint N>
    class QoI_MultiQoIs : public concepts::QoIInterface<T, N> {

    public:
        QoI_MultiQoIs(std::vector<std::shared_ptr<concepts::QoIInterface< T, N>>> qois){
            // TODO: adjust type for int data
            this->qois = qois;
            concepts::QoIInterface<T, N>::id = 999;
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            T eb = qois[0]->interpret_eb(data);
            for(int i=1; i<qois.size(); i++){
                eb = std::min(eb, qois[i]->interpret_eb(data));
            }
            return eb;
        }

        T interpret_eb(const iterator &iter) const {
            T eb = qois[0]->interpret_eb(iter);
            for(int i=1; i<qois.size(); i++){
                eb = std::min(eb, qois[i]->interpret_eb(iter));
            }
            return eb;
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            T eb = qois[0]->interpret_eb(data, offset);
            for(int i=1; i<qois.size(); i++){
                eb = std::min(eb, qois[i]->interpret_eb(data, offset));
            }
            return eb;
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            for(int i=0; i<qois.size(); i++){
                if(!qois[i]->check_compliance(data, dec_data)) return false;
            }
            return true;
        }

        void update_tolerance(T data, T dec_data){
            for(int i=0; i<qois.size(); i++){
                qois[i]->update_tolerance(data, dec_data);
            }
        }

        void precompress_block(const std::shared_ptr<Range> &range){
            for(int i=0; i<qois.size(); i++){
                qois[i]->precompress_block(range);
            }            
        }

        void postcompress_block(){
            for(int i=0; i<qois.size(); i++){
                qois[i]->postcompress_block();
            }                        
        }

        void print(){}

        T get_global_eb() const {
            return qois[0]->get_global_eb();                         
        }

        void set_global_eb(T eb) {
            for(int i=0; i<qois.size(); i++){
                qois[i]->set_global_eb(eb);
            }                                    
        }

        void init(){
            for(int i=0; i<qois.size(); i++){
                qois[i]->init();
            }                                                
        }

        void set_dims(const std::vector<size_t>& new_dims){
            for(int i=0; i<qois.size(); i++){
                qois[i]->set_dims(new_dims);
            }                                                            
        }

    private:
        std::vector<std::shared_ptr<concepts::QoIInterface< T, N>>> qois;
    };
}
#endif 