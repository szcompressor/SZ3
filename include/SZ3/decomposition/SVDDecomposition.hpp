#ifndef SZ3_SVD_DECOMPOSITION_HPP
#define SZ3_SVD_DECOMPOSITION_HPP

#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <unsupported/Eigen/CXX11/Tensor>
#include <iostream>

namespace SZ3 {

template<class T, uint N, class Quantizer>
class SVDDecomposition : public concepts::DecompositionInterface<T, int, N> {
public:
    SVDDecomposition(const Config& conf, Quantizer svd_quantizer, Quantizer res_quantizer) : config(conf), svd_quantizer(svd_quantizer), res_quantizer(res_quantizer) {}

    ~SVDDecomposition() override = default;

    std::vector<int> compress(const Config& conf, T* data) override {
        // 1. SVD Decomposition
        Eigen::array<Eigen::Index, N> tensor_dims;
        for (uint i = 0; i < N; ++i) {
            tensor_dims[i] = conf.dims[i];
        }
        Eigen::Tensor<T, N> tensor(tensor_dims);
        for (size_t i = 0; i < conf.num; ++i) {
            tensor.data()[i] = data[i];
        }

        Eigen::Tensor<T, N> temp_core_tensor = tensor;
        std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> temp_factor_matrices;
        for (uint n = 0; n < N; ++n) {
            auto unfolded_matrix = unfold(temp_core_tensor, n);
            int rank = config.svd_target_rank > 0 ? config.svd_target_rank : unfolded_matrix.cols() / 4;
            auto U = perform_rsvd(unfolded_matrix, rank);
            temp_factor_matrices.push_back(U);
            Eigen::Tensor<T, N> contracted_tensor = contract_tensor_matrix(temp_core_tensor, U.transpose(), n);
            temp_core_tensor = contracted_tensor;
        }

        // 2. Quantize SVD Tensors
        quantized_core_tensor_data.resize(temp_core_tensor.size());
        for (int i = 0; i < temp_core_tensor.size(); ++i) {
            T val = temp_core_tensor.data()[i];
            quantized_core_tensor_data[i] = svd_quantizer.quantize_and_overwrite(val, 0);
        }
        quantized_factor_matrices_data.resize(temp_factor_matrices.size());
        for (size_t i = 0; i < temp_factor_matrices.size(); ++i) {
            quantized_factor_matrices_data[i].resize(temp_factor_matrices[i].size());
            for (int j = 0; j < temp_factor_matrices[i].size(); ++j) {
                T val = temp_factor_matrices[i].data()[j];
                quantized_factor_matrices_data[i][j] = svd_quantizer.quantize_and_overwrite(val, 0);
            }
        }
        
        // Store dimensions for decompression
        for(uint i=0; i<N; ++i) core_tensor_dims[i] = temp_core_tensor.dimension(i);
        for(const auto& m : temp_factor_matrices) {
            factor_matrices_dims.push_back({(long)m.rows(), (long)m.cols()});
        }


        // 3. Reverse SVD (Reconstruction)
        Eigen::Tensor<T, N> reconstructed_tensor = dequantize_and_reconstruct(svd_quantizer);

        // 4. Calculate Difference
        std::vector<T> diff(conf.num);
        for (size_t i = 0; i < conf.num; ++i) {
            diff[i] = data[i] - reconstructed_tensor.data()[i];
        }

        // 5. Quantize Difference
        std::vector<int> quantized_diff(conf.num);
        for (size_t i = 0; i < conf.num; ++i) {
            quantized_diff[i] = res_quantizer.quantize_and_overwrite(diff[i], 0);
        }

        return quantized_diff;
    }

    T* decompress(const Config& conf, std::vector<int>& quant_inds, T* dec_data) override {

        // 1. Dequantize Difference
        std::vector<T> diff(conf.num);
        for (size_t i = 0; i < conf.num; ++i) {
            diff[i] = res_quantizer.recover(0, quant_inds[i]);
        }

        // 2. Reconstruct from SVD
        Eigen::Tensor<T, N> reconstructed_tensor = dequantize_and_reconstruct(svd_quantizer);
        // 3. Add Difference
        for (size_t i = 0; i < conf.num; ++i) {
            dec_data[i] = reconstructed_tensor.data()[i] + diff[i];
        }

        return dec_data;
    }

    void save(uchar*& c) override {
        write(static_cast<uint>(quantized_factor_matrices_data.size()), c);
        for (size_t i = 0; i < quantized_factor_matrices_data.size(); ++i) {
            write(factor_matrices_dims[i][0], c);
            write(factor_matrices_dims[i][1], c);
            write(quantized_factor_matrices_data[i].data(), quantized_factor_matrices_data[i].size(), c);
        }
        write(static_cast<uint>(core_tensor_dims.size()), c);
        for (int i = 0; i < core_tensor_dims.size(); ++i) {
            write(core_tensor_dims[i], c);
        }
        write(quantized_core_tensor_data.data(), quantized_core_tensor_data.size(), c);
        svd_quantizer.save(c);
        res_quantizer.save(c);
    }

    void load(const uchar*& c, size_t& remaining_length) override {
        uint num_factor_matrices;
        read(num_factor_matrices, c);
        quantized_factor_matrices_data.resize(num_factor_matrices);
        factor_matrices_dims.resize(num_factor_matrices);
        for (uint i = 0; i < num_factor_matrices; ++i) {
            read(factor_matrices_dims[i][0], c);
            read(factor_matrices_dims[i][1], c);
            quantized_factor_matrices_data[i].resize(factor_matrices_dims[i][0] * factor_matrices_dims[i][1]);
            read(quantized_factor_matrices_data[i].data(), quantized_factor_matrices_data[i].size(), c);
        }
        uint rank;
        read(rank, c);
        for (uint i = 0; i < rank; ++i) {
            read(core_tensor_dims[i], c);
        }
        size_t core_size = 1;
        for(long dim : core_tensor_dims) core_size *= dim;
        quantized_core_tensor_data.resize(core_size);
        read(quantized_core_tensor_data.data(), quantized_core_tensor_data.size(), c);
        svd_quantizer.load(c, remaining_length);
        res_quantizer.load(c, remaining_length);
    }

    std::pair<int, int> get_out_range() override { return res_quantizer.get_out_range(); }

    void print() override { std::cout << "SVDDecomposition" << std::endl; }

private:
    Config config;
    std::vector<std::vector<int>> quantized_factor_matrices_data;
    std::vector<int> quantized_core_tensor_data;
    Eigen::array<long, N> core_tensor_dims;
    std::vector<std::array<long, 2>> factor_matrices_dims;
    Quantizer svd_quantizer;
    Quantizer res_quantizer;


    Eigen::Tensor<T, 2> matrix_to_tensor(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix) {
        Eigen::Tensor<T, 2> tensor(matrix.rows(), matrix.cols());
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                tensor(i, j) = matrix(i, j);
            }
        }
        return tensor;
    }

    Eigen::Tensor<T, N> contract_tensor_matrix(const Eigen::Tensor<T, N>& tensor, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix, int mode) {
        Eigen::Tensor<T, 2> matrix_tensor = matrix_to_tensor(matrix);
        Eigen::array<Eigen::IndexPair<int>, 1> product_dims = { Eigen::IndexPair<int>(mode, 1) };
        return tensor.contract(matrix_tensor, product_dims);
    }

    Eigen::Tensor<T, N> dequantize_and_reconstruct(Quantizer quantizer) {
        Eigen::Tensor<T, N> core_tensor(core_tensor_dims);
        for(size_t i=0; i<quantized_core_tensor_data.size(); ++i) {
            core_tensor.data()[i] = quantizer.recover(0, quantized_core_tensor_data[i]);
        }

        std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> factor_matrices(quantized_factor_matrices_data.size());
        for(size_t i=0; i<quantized_factor_matrices_data.size(); ++i) {
            factor_matrices[i].resize(factor_matrices_dims[i][0], factor_matrices_dims[i][1]);
            for(size_t j=0; j<quantized_factor_matrices_data[i].size(); ++j) {
                factor_matrices[i].data()[j] = quantizer.recover(0, quantized_factor_matrices_data[i][j]);
            }
        }

        Eigen::Tensor<T, N> reconstructed_tensor = core_tensor;
        for (int n = N - 1; n >= 0; --n) {
            Eigen::Tensor<T, N> contracted_tensor = contract_tensor_matrix(reconstructed_tensor, factor_matrices[n], n);
            reconstructed_tensor = contracted_tensor;
        }
        return reconstructed_tensor;
    }

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> unfold(const Eigen::Tensor<T, N>& tensor, int mode) {
        Eigen::array<int, N> shuffle;
        std::iota(shuffle.begin(), shuffle.end(), 0);
        std::swap(shuffle[0], shuffle[mode]);
        auto shuffled_tensor = tensor.shuffle(shuffle);

        long rows = tensor.dimension(mode);
        long cols = tensor.size() / rows;
        
        Eigen::array<Eigen::Index, 2> unfolded_dims = {rows, cols};
        Eigen::Tensor<T, 2> unfolded_tensor = shuffled_tensor.reshape(unfolded_dims);

        return Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(unfolded_tensor.data(), rows, cols);
    }

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> perform_rsvd(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, int rank) {
        if (rank >= A.cols()) {
            Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(A, Eigen::ComputeThinU);
            return svd.matrixU();
        }

        int l = rank + config.svd_oversampling_param;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omega = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Random(A.cols(), l);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Y = A * Omega;
        Eigen::HouseholderQR<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> qr(Y);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Q = qr.householderQ() * Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(Y.rows(), l);

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B = Q.transpose() * A;
        Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd_B(B, Eigen::ComputeThinU);
        
        return Q * svd_B.matrixU();
    }
};

} // namespace SZ3

#endif
