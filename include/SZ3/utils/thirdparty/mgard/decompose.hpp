#ifndef SZ3_THIRDPARTY_MGARD_DECOMPOSE_HPP
#define SZ3_THIRDPARTY_MGARD_DECOMPOSE_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "SZ3/utils/thirdparty/mgard/correction.hpp"
#include "SZ3/utils/thirdparty/mgard/reorder.hpp"
#include "SZ3/utils/thirdparty/mgard/utils.hpp"

namespace SZ3 {
namespace MGARD {

template <class T>
class Decomposer {
   public:
    Decomposer() = default;

    ~Decomposer() {
        if (data_buffer) std::free(data_buffer);
        if (correction_buffer) std::free(correction_buffer);
        if (load_v_buffer) std::free(load_v_buffer);
    }

    Decomposer(const Decomposer&) = delete;
    Decomposer& operator=(const Decomposer&) = delete;

    /**
     * @brief Apply the MGARD multigrid forward transform in-place.
     *
     * Walks `target_level` levels from finest to coarsest. After each level the
     * coarse nodal block is overwritten with smoothed values; the complementary
     * coefficient block is overwritten with the L2-orthogonal multilevel detail.
     */
    int decompose(T* data_, const std::vector<size_t>& dims, size_t target_level, bool hierarchical = false,
                  std::vector<size_t> strides = std::vector<size_t>()) {
        data = data_;
        size_t num_elements = 1;
        for (const auto& d : dims) {
            num_elements *= d;
        }
        if (strides.empty()) {
            strides = std::vector<size_t>(dims.size());
            size_t stride = 1;
            for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(dims.size()) - 1; i >= 0; i--) {
                strides[i] = stride;
                stride *= dims[i];
            }
        }
        data_buffer_size = num_elements * sizeof(T);
        size_t min_dim = *std::min_element(dims.begin(), dims.end());
        size_t max_level = (min_dim > 1) ? static_cast<size_t>(std::log2(static_cast<double>(min_dim))) : 0;
        if (target_level > max_level) target_level = max_level;
        init(dims);
        if (dims.size() == 1) {
            size_t h = 1;
            size_t n = dims[0];
            for (size_t i = 0; i < target_level; i++) {
                hierarchical ? decompose_level_1D_with_hierarchical_basis(data, n, T(h))
                             : decompose_level_1D(data, n, T(h));
                n = (n >> 1) + 1;
                h <<= 1;
            }
        } else if (dims.size() == 2) {
            size_t h = 1;
            size_t n1 = dims[0];
            size_t n2 = dims[1];
            for (size_t i = 0; i < target_level; i++) {
                hierarchical ? decompose_level_2D_with_hierarchical_basis(data, n1, n2, T(h), strides[0])
                             : decompose_level_2D(data, n1, n2, T(h), strides[0]);
                n1 = (n1 >> 1) + 1;
                n2 = (n2 >> 1) + 1;
                h <<= 1;
            }
        } else if (dims.size() == 3) {
            size_t h = 1;
            size_t n1 = dims[0];
            size_t n2 = dims[1];
            size_t n3 = dims[2];
            for (size_t i = 0; i < target_level; i++) {
                hierarchical ? decompose_level_3D_with_hierarchical_basis(data, n1, n2, n3, T(h), strides[0],
                                                                          strides[1])
                             : decompose_level_3D(data, n1, n2, n3, T(h), strides[0], strides[1]);
                n1 = (n1 >> 1) + 1;
                n2 = (n2 >> 1) + 1;
                n3 = (n3 >> 1) + 1;
                h <<= 1;
            }
        }
        return static_cast<int>(target_level);
    }

   private:
    static constexpr unsigned int default_batch_size = 32;
    size_t data_buffer_size = 0;
    T* data = nullptr;
    T* data_buffer = nullptr;
    T* load_v_buffer = nullptr;
    T* correction_buffer = nullptr;

    void init(const std::vector<size_t>& dims) {
        size_t buffer_size = default_batch_size * (*std::max_element(dims.begin(), dims.end())) * sizeof(T);
        if (data_buffer) std::free(data_buffer);
        if (correction_buffer) std::free(correction_buffer);
        if (load_v_buffer) std::free(load_v_buffer);
        data_buffer = static_cast<T*>(std::malloc(data_buffer_size));
        correction_buffer = static_cast<T*>(std::malloc(buffer_size));
        load_v_buffer = static_cast<T*>(std::malloc(buffer_size));
    }

    void compute_interpolant_difference_1D(size_t n_coeff, const T* nodal_buffer, T* coeff_buffer) {
        for (size_t i = 0; i < n_coeff; i++) {
            coeff_buffer[i] -= (nodal_buffer[i] + nodal_buffer[i + 1]) / 2;
        }
    }

    void add_correction(size_t n_nodal, T* nodal_buffer) {
        for (size_t i = 0; i < n_nodal; i++) {
            nodal_buffer[i] += correction_buffer[i];
        }
    }

    void decompose_level_1D(T* data_pos, size_t n, T h, bool nodal_row = true) {
        size_t n_nodal = (n >> 1) + 1;
        size_t n_coeff = n - n_nodal;
        T* nodal_buffer = data_buffer;
        T* coeff_buffer = data_buffer + n_nodal;
        data_reorder_1D(data_pos, n_nodal, n_coeff, nodal_buffer, coeff_buffer);
        compute_interpolant_difference_1D(n_coeff, nodal_buffer, coeff_buffer);
        if (nodal_row)
            compute_load_vector_nodal_row(load_v_buffer, n_nodal, n_coeff, h, coeff_buffer);
        else
            compute_load_vector_coeff_row(load_v_buffer, n_nodal, n_coeff, h, nodal_buffer, coeff_buffer);
        compute_correction(correction_buffer, n_nodal, h, load_v_buffer);
        add_correction(n_nodal, nodal_buffer);
        std::memcpy(data_pos, data_buffer, n * sizeof(T));
    }

    void decompose_level_1D_with_hierarchical_basis(T* data_pos, size_t n, T /*h*/, bool /*nodal_row*/ = true) {
        size_t n_nodal = (n >> 1) + 1;
        size_t n_coeff = n - n_nodal;
        T* nodal_buffer = data_buffer;
        T* coeff_buffer = data_buffer + n_nodal;
        data_reorder_1D(data_pos, n_nodal, n_coeff, nodal_buffer, coeff_buffer);
        compute_interpolant_difference_1D(n_coeff, nodal_buffer, coeff_buffer);
        std::memcpy(data_pos, data_buffer, n * sizeof(T));
    }

    void compute_interpolant_difference_2D_vertical(T* data_pos, size_t n1, size_t n2, size_t stride) {
        size_t n1_nodal = (n1 >> 1) + 1;
        size_t n1_coeff = n1 - n1_nodal;
        size_t n2_nodal = (n2 >> 1) + 1;
        size_t n2_coeff = n2 - n2_nodal;
        bool even_n2 = !(n2 & 1);
        T* n1_nodal_data = data_pos;
        T* n1_coeff_data = data_pos + n1_nodal * stride;
        for (size_t i = 0; i < n1_coeff; i++) {
            const T* nodal_pos = n1_nodal_data + i * stride;
            T* coeff_pos = n1_coeff_data + i * stride;
            T* nodal_coeff_pos = coeff_pos;
            T* coeff_coeff_pos = coeff_pos + n2_nodal;
            for (size_t j = 0; j < n2_coeff; j++) {
                *(nodal_coeff_pos++) -= (nodal_pos[j] + nodal_pos[stride + j]) / 2;
                *(coeff_coeff_pos++) -= (nodal_pos[j] + nodal_pos[j + 1] + nodal_pos[stride + j] +
                                         nodal_pos[stride + j + 1]) / 4;
            }
            *(nodal_coeff_pos++) -= (nodal_pos[n2_coeff] + nodal_pos[stride + n2_coeff]) / 2;
            if (even_n2) {
                *(nodal_coeff_pos++) -= (nodal_pos[n2_coeff + 1] + nodal_pos[stride + n2_coeff + 1]) / 2;
            }
        }
    }

    void compute_interpolant_difference_2D(T* data_pos, size_t n1, size_t n2, size_t stride) {
        size_t n1_nodal = (n1 >> 1) + 1;
        size_t n2_nodal = (n2 >> 1) + 1;
        size_t n2_coeff = n2 - n2_nodal;
        const T* nodal_pos = data_pos;
        T* coeff_pos = data_pos + n2_nodal;
        for (size_t i = 0; i < n1_nodal; i++) {
            compute_interpolant_difference_1D(n2_coeff, nodal_pos, coeff_pos);
            nodal_pos += stride;
            coeff_pos += stride;
        }
        compute_interpolant_difference_2D_vertical(data_pos, n1, n2, stride);
    }

    void decompose_level_2D(T* data_pos, size_t n1, size_t n2, T h, size_t stride) {
        size_t n1_nodal = (n1 >> 1) + 1;
        size_t n2_nodal = (n2 >> 1) + 1;
        data_reorder_2D(data_pos, data_buffer, n1, n2, stride);
        compute_interpolant_difference_2D(data_pos, n1, n2, stride);
        std::vector<T> w1(n1_nodal);
        std::vector<T> b1(n1_nodal);
        std::vector<T> w2(n2_nodal);
        std::vector<T> b2(n2_nodal);
        precompute_w_and_b(w1.data(), b1.data(), n1_nodal);
        precompute_w_and_b(w2.data(), b2.data(), n2_nodal);
        compute_correction_2D(data_pos, data_buffer, load_v_buffer, n1, n2, n1_nodal, h, stride, w1.data(), b1.data(),
                              w2.data(), b2.data(), default_batch_size);
        apply_correction_batched(data_pos, data_buffer, static_cast<int>(n1_nodal), static_cast<int>(stride),
                                 static_cast<int>(n2_nodal), true);
    }

    void decompose_level_2D_with_hierarchical_basis(T* data_pos, size_t n1, size_t n2, T /*h*/, size_t stride) {
        data_reorder_2D(data_pos, data_buffer, n1, n2, stride);
        compute_interpolant_difference_2D(data_pos, n1, n2, stride);
    }

    void compute_interpolant_difference_3D(T* data_pos, size_t n1, size_t n2, size_t n3, size_t dim0_stride,
                                           size_t dim1_stride) {
        size_t n1_nodal = (n1 >> 1) + 1;
        size_t n1_coeff = n1 - n1_nodal;
        size_t n2_nodal = (n2 >> 1) + 1;
        size_t n2_coeff = n2 - n2_nodal;
        size_t n3_nodal = (n3 >> 1) + 1;
        size_t n3_coeff = n3 - n3_nodal;
        bool even_n2 = (!(n2 & 1));
        bool even_n3 = (!(n3 & 1));
        T* cur_data_pos = data_pos;
        for (size_t i = 0; i < n1_nodal; i++) {
            compute_interpolant_difference_2D(cur_data_pos, n2, n3, dim1_stride);
            cur_data_pos += dim0_stride;
        }
        const T* nodal_pos = data_pos;
        T* coeff_pos = data_pos + n1_nodal * dim0_stride;
        for (size_t i = 0; i < n1_coeff; i++) {
            const T* nodal_nodal_nodal_pos = nodal_pos;
            T* coeff_nodal_nodal_pos = coeff_pos;
            T* coeff_nodal_coeff_pos = coeff_pos + n3_nodal;
            T* coeff_coeff_nodal_pos = coeff_pos + n2_nodal * dim1_stride;
            T* coeff_coeff_coeff_pos = coeff_coeff_nodal_pos + n3_nodal;
            for (size_t j = 0; j < n2_coeff; j++) {
                for (size_t k = 0; k < n3_coeff; k++) {
                    coeff_nodal_nodal_pos[k] -=
                        (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k]) / 2;
                    coeff_nodal_coeff_pos[k] -=
                        (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k] +
                         nodal_nodal_nodal_pos[k + 1] + nodal_nodal_nodal_pos[dim0_stride + k + 1]) / 4;
                    coeff_coeff_nodal_pos[k] -=
                        (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k] +
                         nodal_nodal_nodal_pos[k + dim1_stride] +
                         nodal_nodal_nodal_pos[dim0_stride + k + dim1_stride]) / 4;
                    coeff_coeff_coeff_pos[k] -=
                        (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k] +
                         nodal_nodal_nodal_pos[k + 1] + nodal_nodal_nodal_pos[dim0_stride + k + 1] +
                         nodal_nodal_nodal_pos[k + dim1_stride] +
                         nodal_nodal_nodal_pos[dim0_stride + k + dim1_stride] +
                         nodal_nodal_nodal_pos[k + dim1_stride + 1] +
                         nodal_nodal_nodal_pos[dim0_stride + k + dim1_stride + 1]) / 8;
                }
                coeff_nodal_nodal_pos[n3_coeff] -=
                    (nodal_nodal_nodal_pos[n3_coeff] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff]) / 2;
                coeff_coeff_nodal_pos[n3_coeff] -=
                    (nodal_nodal_nodal_pos[n3_coeff] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff] +
                     nodal_nodal_nodal_pos[n3_coeff + dim1_stride] +
                     nodal_nodal_nodal_pos[dim0_stride + n3_coeff + dim1_stride]) / 4;
                if (even_n3) {
                    coeff_nodal_nodal_pos[n3_coeff + 1] -=
                        (nodal_nodal_nodal_pos[n3_coeff + 1] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]) / 2;
                    coeff_coeff_nodal_pos[n3_coeff + 1] -=
                        (nodal_nodal_nodal_pos[n3_coeff + 1] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1] +
                         nodal_nodal_nodal_pos[n3_coeff + 1 + dim1_stride] +
                         nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1 + dim1_stride]) / 4;
                }
                coeff_nodal_nodal_pos += dim1_stride;
                coeff_nodal_coeff_pos += dim1_stride;
                coeff_coeff_nodal_pos += dim1_stride;
                coeff_coeff_coeff_pos += dim1_stride;
                nodal_nodal_nodal_pos += dim1_stride;
            }
            {
                for (size_t k = 0; k < n3_coeff; k++) {
                    coeff_nodal_nodal_pos[k] -=
                        (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k]) / 2;
                    coeff_nodal_coeff_pos[k] -=
                        (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k] +
                         nodal_nodal_nodal_pos[k + 1] + nodal_nodal_nodal_pos[dim0_stride + k + 1]) / 4;
                }
                coeff_nodal_nodal_pos[n3_coeff] -=
                    (nodal_nodal_nodal_pos[n3_coeff] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff]) / 2;
                if (even_n3) {
                    coeff_nodal_nodal_pos[n3_coeff + 1] -=
                        (nodal_nodal_nodal_pos[n3_coeff + 1] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]) / 2;
                }
                coeff_nodal_nodal_pos += dim1_stride;
                coeff_nodal_coeff_pos += dim1_stride;
                coeff_coeff_nodal_pos += dim1_stride;
                coeff_coeff_coeff_pos += dim1_stride;
                nodal_nodal_nodal_pos += dim1_stride;
            }
            if (even_n2) {
                for (size_t k = 0; k < n3_coeff; k++) {
                    coeff_nodal_nodal_pos[k] -=
                        (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k]) / 2;
                    coeff_nodal_coeff_pos[k] -=
                        (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k] +
                         nodal_nodal_nodal_pos[k + 1] + nodal_nodal_nodal_pos[dim0_stride + k + 1]) / 4;
                }
                coeff_nodal_nodal_pos[n3_coeff] -=
                    (nodal_nodal_nodal_pos[n3_coeff] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff]) / 2;
                if (even_n3) {
                    coeff_nodal_nodal_pos[n3_coeff + 1] -=
                        (nodal_nodal_nodal_pos[n3_coeff + 1] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]) / 2;
                }
            }
            nodal_pos += dim0_stride;
            coeff_pos += dim0_stride;
        }
    }

    void decompose_level_3D(T* data_pos, size_t n1, size_t n2, size_t n3, T h, size_t dim0_stride,
                            size_t dim1_stride) {
        data_reorder_3D(data_pos, data_buffer, n1, n2, n3, dim0_stride, dim1_stride);
        compute_interpolant_difference_3D(data_pos, n1, n2, n3, dim0_stride, dim1_stride);
        size_t n1_nodal = (n1 >> 1) + 1;
        size_t n2_nodal = (n2 >> 1) + 1;
        size_t n3_nodal = (n3 >> 1) + 1;
        compute_correction_3D(data_pos, data_buffer, load_v_buffer, n1, n2, n3, n1_nodal, h, dim0_stride, dim1_stride,
                              default_batch_size);
        T* nodal_pos = data_pos;
        const T* correction_pos = data_buffer;
        for (size_t i = 0; i < n1_nodal; i++) {
            apply_correction_batched(nodal_pos, correction_pos, static_cast<int>(n2_nodal),
                                     static_cast<int>(dim1_stride), static_cast<int>(n3_nodal), true);
            nodal_pos += dim0_stride;
            correction_pos += n2_nodal * n3_nodal;
        }
    }

    void decompose_level_3D_with_hierarchical_basis(T* data_pos, size_t n1, size_t n2, size_t n3, T /*h*/,
                                                    size_t dim0_stride, size_t dim1_stride) {
        data_reorder_3D(data_pos, data_buffer, n1, n2, n3, dim0_stride, dim1_stride);
        compute_interpolant_difference_3D(data_pos, n1, n2, n3, dim0_stride, dim1_stride);
    }
};

}  // namespace MGARD
}  // namespace SZ3

#endif
