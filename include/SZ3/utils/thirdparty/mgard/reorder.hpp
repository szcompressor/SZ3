#ifndef SZ3_THIRDPARTY_MGARD_REORDER_HPP
#define SZ3_THIRDPARTY_MGARD_REORDER_HPP

#include <cstddef>
#include <cstring>

namespace SZ3 {
namespace MGARD {

template <class T>
void switch_rows_2D_by_buffer(T* data_pos, T* data_buffer, size_t n1, size_t n2, size_t stride) {
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    T* nodal_data_buffer = data_buffer + n2;
    T* coeff_data_buffer = data_buffer + n1_nodal * n2;
    T* cur_data_pos = data_pos + stride;
    for (size_t i = 0; i < n1_coeff; i++) {
        std::memcpy(coeff_data_buffer + i * n2, cur_data_pos, n2 * sizeof(T));
        cur_data_pos += stride;
        std::memcpy(nodal_data_buffer + i * n2, cur_data_pos, n2 * sizeof(T));
        cur_data_pos += stride;
    }
    if (!(n1 & 1)) {
        std::memcpy(coeff_data_buffer - n2, cur_data_pos, n2 * sizeof(T));
    }
    cur_data_pos = data_pos + stride;
    for (size_t i = 1; i < n1; i++) {
        std::memcpy(cur_data_pos, data_buffer + i * n2, n2 * sizeof(T));
        cur_data_pos += stride;
    }
}

template <class T>
void switch_rows_2D_by_buffer_reverse(T* data_pos, T* data_buffer, size_t n1, size_t n2, size_t stride) {
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    T* nodal_data_buffer = data_pos + stride;
    T* coeff_data_buffer = data_pos + n1_nodal * stride;
    T* cur_data_pos = data_buffer + n2;
    for (size_t i = 0; i < n1_coeff; i++) {
        std::memcpy(cur_data_pos, coeff_data_buffer + i * stride, n2 * sizeof(T));
        cur_data_pos += n2;
        std::memcpy(cur_data_pos, nodal_data_buffer + i * stride, n2 * sizeof(T));
        cur_data_pos += n2;
    }
    if (!(n1 & 1)) {
        std::memcpy(cur_data_pos, coeff_data_buffer - stride, n2 * sizeof(T));
    }
    cur_data_pos = data_pos + stride;
    for (size_t i = 1; i < n1; i++) {
        std::memcpy(cur_data_pos, data_buffer + i * n2, n2 * sizeof(T));
        cur_data_pos += stride;
    }
}

template <class T>
void data_reorder_1D(const T* data_pos, size_t n_nodal, size_t n_coeff, T* nodal_buffer, T* coeff_buffer) {
    T* nodal_pos = nodal_buffer;
    T* coeff_pos = coeff_buffer;
    const T* cur_data_pos = data_pos;
    for (size_t i = 0; i < n_coeff; i++) {
        *(nodal_pos++) = *(cur_data_pos++);
        *(coeff_pos++) = *(cur_data_pos++);
    }
    *(nodal_pos++) = *(cur_data_pos++);
    if (n_nodal == n_coeff + 2) {
        // even case: synthesize an extra nodal value so the last two nodals
        // interpolate to the trailing coefficient.
        *nodal_pos = 2 * cur_data_pos[0] - nodal_pos[-1];
    }
}

template <class T>
void data_reorder_2D(T* data_pos, T* data_buffer, size_t n1, size_t n2, size_t stride) {
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    T* cur_data_pos = data_pos;
    T* nodal_pos = data_buffer;
    T* coeff_pos = data_buffer + n2_nodal;
    for (size_t i = 0; i < n1; i++) {
        data_reorder_1D(cur_data_pos, n2_nodal, n2_coeff, nodal_pos, coeff_pos);
        std::memcpy(cur_data_pos, data_buffer, n2 * sizeof(T));
        cur_data_pos += stride;
    }
    if (!(n1 & 1)) {
        cur_data_pos -= stride;
        for (size_t j = 0; j < n2; j++) {
            cur_data_pos[j] = 2 * cur_data_pos[j] - cur_data_pos[-static_cast<std::ptrdiff_t>(stride) + j];
        }
    }
    switch_rows_2D_by_buffer(data_pos, data_buffer, n1, n2, stride);
}

template <class T>
void data_reorder_3D(T* data_pos, T* data_buffer, size_t n1, size_t n2, size_t n3, size_t dim0_stride,
                     size_t dim1_stride) {
    T* cur_data_pos = data_pos;
    for (size_t i = 0; i < n1; i++) {
        data_reorder_2D(cur_data_pos, data_buffer, n2, n3, dim1_stride);
        cur_data_pos += dim0_stride;
    }
    if (!(n1 & 1)) {
        cur_data_pos -= dim0_stride;
        for (size_t j = 0; j < n2; j++) {
            for (size_t k = 0; k < n3; k++) {
                cur_data_pos[k] = 2 * cur_data_pos[k] - cur_data_pos[-static_cast<std::ptrdiff_t>(dim0_stride) + k];
            }
            cur_data_pos += dim1_stride;
        }
    }
    cur_data_pos = data_pos;
    for (size_t j = 0; j < n2; j++) {
        switch_rows_2D_by_buffer(cur_data_pos, data_buffer, n1, n3, dim0_stride);
        cur_data_pos += dim1_stride;
    }
}

template <class T>
void data_reverse_reorder_1D(T* data_pos, size_t n_nodal, size_t n_coeff, const T* nodal_buffer,
                             const T* coeff_buffer) {
    const T* nodal_pos = nodal_buffer;
    const T* coeff_pos = coeff_buffer;
    T* cur_data_pos = data_pos;
    for (size_t i = 0; i < n_coeff; i++) {
        *(cur_data_pos++) = *(nodal_pos++);
        *(cur_data_pos++) = *(coeff_pos++);
    }
    *(cur_data_pos++) = *(nodal_pos++);
    if (n_nodal == n_coeff + 2) {
        *cur_data_pos = (nodal_pos[-1] + nodal_pos[0]) / 2;
    }
}

template <class T>
void data_reverse_reorder_2D(T* data_pos, T* data_buffer, size_t n1, size_t n2, size_t stride) {
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    T* cur_data_pos = data_pos;
    T* nodal_pos = data_buffer;
    T* coeff_pos = data_buffer + n2_nodal;
    switch_rows_2D_by_buffer_reverse(data_pos, data_buffer, n1, n2, stride);
    for (size_t i = 0; i < n1; i++) {
        std::memcpy(data_buffer, cur_data_pos, n2 * sizeof(T));
        data_reverse_reorder_1D(cur_data_pos, n2_nodal, n2_coeff, nodal_pos, coeff_pos);
        cur_data_pos += stride;
    }
    if (!(n1 & 1)) {
        cur_data_pos -= stride;
        for (size_t j = 0; j < n2; j++) {
            cur_data_pos[j] = (cur_data_pos[j] + cur_data_pos[-static_cast<std::ptrdiff_t>(stride) + j]) / 2;
        }
    }
}

template <class T>
void data_reverse_reorder_3D(T* data_pos, T* data_buffer, size_t n1, size_t n2, size_t n3, size_t dim0_stride,
                             size_t dim1_stride) {
    T* cur_data_pos = data_pos;
    for (size_t j = 0; j < n2; j++) {
        switch_rows_2D_by_buffer_reverse(cur_data_pos, data_buffer, n1, n3, dim0_stride);
        cur_data_pos += dim1_stride;
    }
    cur_data_pos = data_pos;
    for (size_t i = 0; i < n1; i++) {
        data_reverse_reorder_2D(cur_data_pos, data_buffer, n2, n3, dim1_stride);
        cur_data_pos += dim0_stride;
    }
    if (!(n1 & 1)) {
        cur_data_pos -= dim0_stride;
        for (size_t j = 0; j < n2; j++) {
            for (size_t k = 0; k < n3; k++) {
                cur_data_pos[k] =
                    (cur_data_pos[k] + cur_data_pos[-static_cast<std::ptrdiff_t>(dim0_stride) + k]) / 2;
            }
            cur_data_pos += dim1_stride;
        }
    }
}

}  // namespace MGARD
}  // namespace SZ3

#endif
