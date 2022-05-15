//
// Created by Xin Liang on 03/09/2021.
//

#ifndef SZ_QOI_REGIONAL_AVERAGE_OF_SQUARE_HPP
#define SZ_QOI_REGIONAL_AVERAGE_OF_SQUARE_HPP

#include <algorithm>
#include "SZ3/def.hpp"
#include "SZ3/qoi/QoI.hpp"
#include "SZ3/utils/Iterator.hpp"

namespace SZ {
    template<class T, uint N>
    class QoI_RegionalAverageOfSquare : public concepts::QoIInterface<T, N> {

    public:
        QoI_RegionalAverageOfSquare(T tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            printf("QoI_RegionalAverageOfSquare\n");
            printf("tolerance = %.4e\n", (double) tolerance);
            printf("global_eb = %.4e\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 3;
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            // eb for data^2
            double eb_x2 = (aggregated_tolerance - fabs(error)) / rest_elements;
            // compute eb based on formula of x^2
            T eb = - fabs(data) + sqrt(data * data + eb_x2);
            return std::min(eb, global_eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        void update_tolerance(T data, T dec_data){
            error += data*data - dec_data*dec_data;
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

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

    private:
        T tolerance;
        T global_eb;
        double error = 0;
        int rest_elements;
        int block_elements;
        double aggregated_tolerance;
    };

    template<class T, uint N>
    class QoI_RegionalAverageOfSquareInterp : public concepts::QoIInterface<T, N> {

    public:
        QoI_RegionalAverageOfSquareInterp(T tolerance, T global_eb, int block_size, std::vector<size_t> dims) : 
                tolerance(tolerance),
                global_eb(global_eb),
                dims(dims),
                block_size(block_size) {
            printf("tolerance = %.4e\n", (double) tolerance);
            printf("global_eb = %.4e\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 3;
            init();
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            std::cerr << "Not implemented\n";
            exit(-1);
            return 0;
        }

        T interpret_eb(const iterator &iter) const {
            std::cerr << "Not implemented\n";
            exit(-1);
            return 0;
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            block_id = compute_block_id(offset);
            double eb_x2 = (aggregated_tolerance[block_id] - fabs(accumulated_error[block_id])) / rest_elements[block_id];
            T data_val = *data;
            T eb = - fabs(data_val) + sqrt(data_val * data_val + eb_x2);
            return std::min(eb, global_eb);
        }

        void update_tolerance(T data, T dec_data){
            accumulated_error[block_id] += data * data - dec_data * dec_data;
            rest_elements[block_id] --;
            if(accumulated_error[block_id] > aggregated_tolerance[block_id]){
                printf("%d: %.4e / %.4e\n", block_id, accumulated_error[block_id], aggregated_tolerance[block_id]);
                printf("%d / %d\n", rest_elements[block_id], block_elements[block_id]);
                exit(-1);
            }
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return true;
        }

        void precompress_block(const std::shared_ptr<Range> &range){}

        void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){
            block_dims = std::vector<size_t>(dims.size());
            size_t num_blocks = 1;
            for(int i=0; i<dims.size(); i++){
                block_dims[i] = (dims[i] - 1) / block_size + 1;
                num_blocks *= block_dims[i];
                std::cout << block_dims[i] << " ";
            }
            std::cout << std::endl;
            aggregated_tolerance = std::vector<double>(num_blocks);
            block_elements = std::vector<int>(num_blocks, 0);
            rest_elements = std::vector<int>(num_blocks, 0);
            if(dims.size() == 2){
                for(int i=0; i<block_dims[0]; i++){
                    int size_x = (i < block_dims[0] - 1) ? block_size : dims[0] - i * block_size;
                    for(int j=0; j<block_dims[1]; j++){
                        int size_y = (j < block_dims[1] - 1) ? block_size : dims[1] - j * block_size;
                        int num_block_elements = size_x * size_y;
                        aggregated_tolerance[i * block_dims[1] + j] = num_block_elements * tolerance;
                        block_elements[i * block_dims[1] + j] = num_block_elements;
                        rest_elements[i * block_dims[1] + j] = num_block_elements;
                    }
                }
            }
            else if(dims.size() == 3){
                for(int i=0; i<block_dims[0]; i++){
                    int size_x = (i < block_dims[0] - 1) ? block_size : dims[0] - i * block_size;
                    for(int j=0; j<block_dims[1]; j++){
                        int size_y = (j < block_dims[1] - 1) ? block_size : dims[1] - j * block_size;
                        for(int k=0; k<block_dims[2]; k++){
                            int size_z = (k < block_dims[2] - 1) ? block_size : dims[2] - k * block_size;
                            int num_block_elements = size_x * size_y * size_z;
                            // printf("%d, %d, %d: %d * %d * %d = %d\n", i, j, k, size_x, size_y, size_z, num_block_elements);
                            aggregated_tolerance[i * block_dims[1] * block_dims[2] + j * block_dims[2] + k] = num_block_elements * tolerance;
                            block_elements[i * block_dims[1] * block_dims[2] + j * block_dims[2] + k] = num_block_elements;
                            rest_elements[i * block_dims[1] * block_dims[2] + j * block_dims[2] + k] = num_block_elements;
                        }
                    }
                }
            }
            else{
                std::cerr << "dims other than 2 or 3 are not implemented" << std::endl;
                exit(-1);
            }

            accumulated_error = std::vector<double>(num_blocks, 0);
            std::cout << "end of init\n";            
        }

        void set_dims(const std::vector<size_t>& new_dims){
            dims = new_dims;
        }

    private:
        template<uint NN = N>
        inline typename std::enable_if<NN == 1, int>::type compute_block_id(ptrdiff_t offset) const noexcept {
            // 1D data
            return offset / block_size;
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 2, int>::type compute_block_id(ptrdiff_t offset) const noexcept {
            // 3D data
            int i = offset / dims[1];
            int j = offset % dims[1];
            return (i / block_size) * block_dims[1] + (j / block_size);
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 3, int>::type compute_block_id(ptrdiff_t offset) const noexcept {
            // 3D data
            int i = offset / (dims[1] * dims[2]);
            offset = offset % (dims[1] * dims[2]);
            int j = offset / dims[2];
            int k = offset % dims[2];
            return (i / block_size) * block_dims[1] * block_dims[2] + (j / block_size) * block_dims[2] + (k / block_size);
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 4, int>::type compute_block_id(ptrdiff_t offset) const noexcept {
            // 4D data
            std::cerr << "Not implemented!\n";
            exit(-1);
            return 0;
        }

        T tolerance;
        T global_eb;
        int block_id;
        int block_size;
        std::vector<double> aggregated_tolerance;
        std::vector<double> accumulated_error;
        std::vector<int> block_elements;
        std::vector<int> rest_elements;
        std::vector<size_t> dims;
        std::vector<size_t> block_dims;
    };


}
#endif 