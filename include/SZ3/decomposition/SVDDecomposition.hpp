/**
 * @file SVDDecomposition.hpp
 * @ingroup Decomposition
 */

#ifndef SZ3_SVD_DECOMPOSITION_HPP
#define SZ3_SVD_DECOMPOSITION_HPP

#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/stdafx.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>

namespace SZ3 {


template<class T, uint N, class Quantizer>
class SVDDecomposition : public concepts::DecompositionInterface<T, int, N> {
public:
    SVDDecomposition(const Config& conf, Quantizer svd_quantizer, Quantizer res_quantizer) : config(conf), svd_quantizer(svd_quantizer), res_quantizer(res_quantizer) {}

    ~SVDDecomposition() override = default;

    std::vector<int> compress(const Config& conf, T* data) override {
        // 1. SVD Decomposition (ST-HOSVD)
        Eigen::array<Eigen::Index, N> tensor_dims;
        for (uint i = 0; i < N; ++i) {
            tensor_dims[i] = conf.dims[i];
        }
        Eigen::Tensor<T, N> temp_core_tensor(tensor_dims);
        memcpy(temp_core_tensor.data(), data, conf.num * sizeof(T));

        std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> temp_factor_matrices;
        
        T tau = config.svd_energy_threshold;

        for (uint n = 0; n < N; ++n) {
            auto unfolded_matrix = unfold(temp_core_tensor, n);
            
            auto rsvd_res = adaptive_rsvd(unfolded_matrix, tau);
            auto U = rsvd_res.first;
            
            temp_factor_matrices.push_back(U);

            temp_core_tensor = contract_and_shuffle(temp_core_tensor, U, n, true);
        }

        // 2. Quantize SVD components and store them
        quantized_core.resize(temp_core_tensor.size());
        core_dims.resize(N);
        for(int i=0; i<N; ++i) core_dims[i] = temp_core_tensor.dimension(i);

        for (int i = 0; i < temp_core_tensor.size(); ++i) {
            T val = temp_core_tensor.data()[i];
            quantized_core[i] = svd_quantizer.quantize_and_overwrite(val, 0);
        }

        quantized_factors.resize(temp_factor_matrices.size());
        factor_dims.resize(temp_factor_matrices.size());
        for (size_t i = 0; i < temp_factor_matrices.size(); ++i) {
            factor_dims[i] = {temp_factor_matrices[i].rows(), temp_factor_matrices[i].cols()};
            quantized_factors[i].resize(temp_factor_matrices[i].size());
            for (int j = 0; j < temp_factor_matrices[i].size(); ++j) {
                T val = temp_factor_matrices[i].data()[j];
                quantized_factors[i][j] = svd_quantizer.quantize_and_overwrite(val, 0);
            }
        }
        
        // 3. Dequantize and reconstruct
        Eigen::Tensor<T, N> reconstructed_tensor = dequantize_and_reconstruct();

        // 4. Calculate residual and quantize it
        std::vector<int> quantized_residual(conf.num);
        for(size_t i=0; i<conf.num; ++i){
            T residual = data[i] - reconstructed_tensor.data()[i];
            quantized_residual[i] = res_quantizer.quantize_and_overwrite(residual, 0);
        }

        return quantized_residual;
    }

    T* decompress(const Config& conf, std::vector<int>& quant_inds, T* dec_data) override {
        // 1. Dequantize residual
        std::vector<T> residual(conf.num);
        for(size_t i=0; i<conf.num; ++i){
            residual[i] = res_quantizer.recover(0, quant_inds[i]);
        }

        // 2. Dequantize SVD components and reconstruct
        Eigen::Tensor<T, N> reconstructed_tensor = dequantize_and_reconstruct();

        // 3. Add residual
        for(size_t i=0; i<conf.num; ++i){
            dec_data[i] = reconstructed_tensor.data()[i] + residual[i];
        }
        return dec_data;
    }

    void save(uchar*& c) override {
        write(core_dims.size(), c);
        write(core_dims.data(), core_dims.size(), c);
        write(quantized_core.size(), c);
        write(quantized_core.data(), quantized_core.size(), c);

        write(factor_dims.size(), c);
        for(const auto& p : factor_dims){
            write(p.first, c);
            write(p.second, c);
        }
        for(const auto& v : quantized_factors){
            write(v.size(), c);
            write(v.data(), v.size(), c);
        }
        svd_quantizer.save(c);
        res_quantizer.save(c);
    }

    void load(const uchar*& c, size_t& remaining_length) override {
        size_t core_dims_size;
        read(core_dims_size, c, remaining_length);
        core_dims.resize(core_dims_size);
        read(core_dims.data(), core_dims_size, c, remaining_length);

        size_t quantized_core_size;
        read(quantized_core_size, c, remaining_length);
        quantized_core.resize(quantized_core_size);
        read(quantized_core.data(), quantized_core_size, c, remaining_length);

        size_t factor_dims_size;
        read(factor_dims_size, c, remaining_length);
        factor_dims.resize(factor_dims_size);
        for(size_t i=0; i<factor_dims_size; ++i){
            read(factor_dims[i].first, c, remaining_length);
            read(factor_dims[i].second, c, remaining_length);
        }

        quantized_factors.resize(factor_dims_size);
        for(size_t i=0; i<factor_dims_size; ++i){
            size_t v_size;
            read(v_size, c, remaining_length);
            quantized_factors[i].resize(v_size);
            read(quantized_factors[i].data(), v_size, c, remaining_length);
        }
        svd_quantizer.load(c, remaining_length);
        res_quantizer.load(c, remaining_length);
    }

    std::pair<int, int> get_out_range() override { return res_quantizer.get_out_range(); }

    void print() override { std::cout << "SVDDecomposition" << std::endl; }

private:
    Config config;
    std::vector<std::vector<int>> quantized_factors;
    std::vector<int> quantized_core;
    std::vector<size_t> core_dims;
    std::vector<std::pair<size_t, size_t>> factor_dims;
    Quantizer svd_quantizer;
    Quantizer res_quantizer;


    // From paper Algorithm 2: Adaptive rSVD (with Rank Finding)
    std::pair<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::VectorX<T>>
    adaptive_rsvd(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, T tau) {
        T energy_A = A.squaredNorm();
        if (energy_A == 0) {
            return {Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(A.rows(), 1), Eigen::VectorX<T>::Zero(1)};
        }
        long min_dim = std::min(A.rows(), A.cols());

        // Predefined sequence of ranks to test
        std::vector<int> ranks_to_test;
        for (int r = 8; r < min_dim; r *= 2) {
            ranks_to_test.push_back(r);
        }
        ranks_to_test.push_back(min_dim);

        for (int rank : ranks_to_test) {
            auto rsvd_res = perform_rsvd(A, rank);
            auto sigma = rsvd_res.second;
            T captured_energy = sigma.squaredNorm();
            if (captured_energy / energy_A >= tau) {
                return rsvd_res;
            }
        }
        // Fallback: return result of largest rank tested
        return perform_rsvd(A, min_dim);
    }

    std::pair<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::VectorX<T>>
    perform_rsvd(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, int rank) {
        // unfolded matrix is short and wide, so A.rows() < A.cols()
        rank = std::min({rank, (int)A.rows(), (int)A.cols()});

        if (rank >= A.rows()) {
            Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
            return {svd.matrixU(), svd.singularValues()};
        }

        int l = rank + config.svd_oversampling_param;
        l = std::min({l, (int)A.rows(), (int)A.cols()});

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Omega = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Random(A.cols(), l);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Y = A * Omega;
        Eigen::HouseholderQR<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> qr(Y);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Q = qr.householderQ() * Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(Y.rows(), l);

        // Gram-based rSVD from paper's Algorithm 1
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B = (Q.transpose() * A) * (Q.transpose() * A).transpose();

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> eigen_solver(B);
        auto eigenvalues = eigen_solver.eigenvalues();
        auto eigenvectors = eigen_solver.eigenvectors();

        // Sort eigenvalues and eigenvectors in descending order
        std::vector<std::pair<T, Eigen::VectorX<T>>> eigen_pairs;
        for (int i = 0; i < eigenvalues.size(); ++i) {
            eigen_pairs.push_back({eigenvalues(i), eigenvectors.col(i)});
        }
        std::sort(eigen_pairs.begin(), eigen_pairs.end(), [](const auto& a, const auto& b) {
            return a.first > b.first;
        });

        Eigen::Vector<T, Eigen::Dynamic> sorted_eigenvalues(l);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> U0(l, l);

        for (int i = 0; i < l; ++i) {
            sorted_eigenvalues(i) = eigen_pairs[i].first > 0 ? eigen_pairs[i].first : 0;
            U0.col(i) = eigen_pairs[i].second;
        }

        Eigen::Vector<T, Eigen::Dynamic> sigma = sorted_eigenvalues.cwiseSqrt();
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> U = Q * U0;

        return {U.leftCols(rank), sigma.head(rank)};
    }

    Eigen::Tensor<T, 2> matrix_to_tensor(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix) {
        Eigen::Tensor<T, 2> tensor(matrix.rows(), matrix.cols());
        memcpy(tensor.data(), matrix.data(), matrix.size() * sizeof(T));
        return tensor;
    }

    Eigen::Tensor<T, N> contract_and_shuffle(const Eigen::Tensor<T, N>& tensor, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix, int mode, bool transpose_matrix) {
        Eigen::Tensor<T, 2> matrix_tensor = transpose_matrix ? matrix_to_tensor(matrix.transpose()) : matrix_to_tensor(matrix);

        Eigen::array<Eigen::IndexPair<int>, 1> product_dims = { Eigen::IndexPair<int>(mode, 1) };
        
        auto contracted_tensor = tensor.contract(matrix_tensor, product_dims);

        std::vector<int> p(N);
        std::iota(p.begin(), p.end(), 0);
        int last = p.back();
        p.pop_back();
        p.insert(p.begin() + mode, last);

        Eigen::array<int, N> shuffle_perm;
        for(int i=0; i<N; ++i) shuffle_perm[i] = p[i];

        return contracted_tensor.shuffle(shuffle_perm);
    }


    Eigen::Tensor<T, N> dequantize_and_reconstruct() {
        Eigen::array<Eigen::Index, N> eigen_core_dims;
        for(int i=0; i<N; ++i) eigen_core_dims[i] = core_dims[i];
        Eigen::Tensor<T, N> core_tensor(eigen_core_dims);
        for(size_t i=0; i<quantized_core.size(); ++i) {
            core_tensor.data()[i] = svd_quantizer.recover(0, quantized_core[i]);
        }

        std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> factor_matrices(quantized_factors.size());
        for(size_t i=0; i<quantized_factors.size(); ++i) {
            factor_matrices[i].resize(factor_dims[i].first, factor_dims[i].second);
            for(size_t j=0; j<quantized_factors[i].size(); ++j) {
                factor_matrices[i].data()[j] = svd_quantizer.recover(0, quantized_factors[i][j]);
            }
        }

        Eigen::Tensor<T, N> reconstructed_tensor = core_tensor;
        for (int mode = N - 1; mode >= 0; --mode) {
            reconstructed_tensor = contract_and_shuffle(reconstructed_tensor, factor_matrices[mode], mode, false);
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

        return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(unfolded_tensor.data(), rows, cols));
    }
};

} // namespace SZ3

#endif
