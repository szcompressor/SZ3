#ifndef _meta_encode_h
#define _meta_encode_h

#include "meta_huffman.hpp"
#include "meta_utils.hpp"

namespace SZMETA {

    int *
    Huffman_decode_tree_and_data(size_t state_num, size_t num_elements, const unsigned char *&compressed_pos) {
        HuffmanTree *huffman = createHuffmanTree(state_num);
        size_t node_count = 0;
        read_variable_from_src(compressed_pos, node_count);
        unsigned int tree_size = 0;
        read_variable_from_src(compressed_pos, tree_size);
        node root = reconstruct_HuffTree_from_bytes_anyStates(huffman, compressed_pos, node_count);
        compressed_pos += tree_size;
        size_t type_array_size = 0;
        read_variable_from_src(compressed_pos, type_array_size);
        int *type = (int *) malloc(num_elements * sizeof(int));
        decode(compressed_pos, num_elements, root, type);
        compressed_pos += type_array_size;
        META_ReleaseHuffman(huffman);
        return type;
    }

    template<typename T>
    float *
    decode_regression_coefficients(const unsigned char *&compressed_pos, size_t reg_count, int block_size, T precision,
                                   const meta_params &params) {
        size_t reg_unpredictable_count = 0;
        read_variable_from_src(compressed_pos, reg_unpredictable_count);
        const float *reg_unpredictable_data_pos = (const float *) compressed_pos;
        compressed_pos += reg_unpredictable_count * sizeof(float);
        int *reg_type = Huffman_decode_tree_and_data(2 * RegCoeffCapacity, RegCoeffNum3d * reg_count, compressed_pos);
        float *reg_params = (float *) malloc(RegCoeffNum3d * (reg_count + 1) * sizeof(float));
        for (int i = 0; i < RegCoeffNum3d; i++)
            reg_params[i] = 0;
        T reg_precisions[RegCoeffNum3d];
        for (int i = 0; i < RegCoeffNum3d - 1; i++) {
            reg_precisions[i] = params.regression_param_eb_linear;
        }
        reg_precisions[RegCoeffNum3d - 1] = params.regression_param_eb_independent;
        float *prev_reg_params = reg_params;
        float *reg_params_pos = reg_params + RegCoeffNum3d;
        const int *type_pos = (const int *) reg_type;
        for (int i = 0; i < reg_count; i++) {
            for (int j = 0; j < RegCoeffNum3d; j++) {
                *reg_params_pos = recover(*prev_reg_params, reg_precisions[j], *(type_pos++), RegCoeffRadius,
                                          reg_unpredictable_data_pos);
                prev_reg_params++, reg_params_pos++;
            }
        }
        free(reg_type);
        return reg_params;
    }

    HuffmanTree *
    build_Huffman_tree(size_t state_num, const int *type, size_t num_elements) {
        HuffmanTree *huffman = createHuffmanTree(state_num);
        init(huffman, type, num_elements);
        return huffman;
    }

    void
    Huffman_encode_tree_and_data(size_t state_num, const int *type, size_t num_elements, unsigned char *&compressed_pos) {
        HuffmanTree *huffman = build_Huffman_tree(state_num, type, num_elements);
        size_t node_count = 0;
        size_t i = 0;
        for (i = 0; i < state_num; i++)
            if (huffman->code[i]) node_count++;
        node_count = node_count * 2 - 1;
        unsigned char *tree_structure = NULL;
        unsigned int tree_size = convert_HuffTree_to_bytes_anyStates(huffman, node_count, &tree_structure);
        write_variable_to_dst(compressed_pos, node_count);
        write_variable_to_dst(compressed_pos, tree_size);
        write_array_to_dst(compressed_pos, tree_structure, tree_size);
        unsigned char *type_array_size_pos = compressed_pos;
        compressed_pos += sizeof(size_t);
        size_t type_array_size = 0;
        encode(huffman, type, num_elements, compressed_pos, &type_array_size);
        write_variable_to_dst(type_array_size_pos, type_array_size);
        compressed_pos += type_array_size;
        free(tree_structure);
        META_ReleaseHuffman(huffman);
    }


    void
    encode_regression_coefficients(const int *reg_params_type, const float *reg_unpredictable_data, size_t reg_count,
                                   size_t reg_unpredictable_count, unsigned char *&compressed_pos) {
        write_variable_to_dst(compressed_pos, reg_unpredictable_count);
        write_array_to_dst(compressed_pos, reg_unpredictable_data, reg_unpredictable_count);
        Huffman_encode_tree_and_data(2 * RegCoeffCapacity, reg_params_type, reg_count, compressed_pos);
    }

}

#endif