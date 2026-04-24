#ifndef SZ3_THIRDPARTY_MGARD_CORRECTION_HPP
#define SZ3_THIRDPARTY_MGARD_CORRECTION_HPP

#include <cstddef>
#include <vector>

namespace SZ3 {
namespace MGARD {

template <class T>
void precompute_w_and_b(T* w, T* b, size_t n_nodal) {
    b[0] = T(2) / 3;
    b[n_nodal - 1] = T(2) / 3;
    for (size_t i = 1; i < n_nodal - 1; i++) {
        b[i] = T(4) / 3;
    }
    T c = T(1) / 3;
    w[0] = 0;
    for (size_t i = 1; i < n_nodal; i++) {
        w[i] = c / b[i - 1];
        b[i] = b[i] - w[i] * c;
    }
}

constexpr double mgard_alpha = 1.0 / 12;
constexpr double mgard_beta = 0.5;
constexpr double mgard_gamma = 5.0 / 6;

template <class T>
void compute_load_vector_nodal_row(T* load_v_buffer, size_t n_nodal, size_t n_coeff, T /*h*/, const T* coeff_buffer) {
    const T* coeff = coeff_buffer;
    T ah = T(mgard_beta);
    load_v_buffer[0] = coeff[0] * ah;
    for (size_t i = 1; i < n_coeff; i++) {
        load_v_buffer[i] = (coeff[i - 1] + coeff[i]) * ah;
    }
    load_v_buffer[n_coeff] = coeff[n_coeff - 1] * ah;
    if (n_nodal == n_coeff + 2) load_v_buffer[n_coeff + 1] = 0;
}

template <class T>
void compute_load_vector_coeff_row(T* load_v_buffer, size_t n_nodal, size_t n_coeff, T /*h*/, const T* nodal_buffer,
                                   const T* coeff_buffer) {
    const T* coeff = coeff_buffer;
    const T* nodal = nodal_buffer;
    T ah = T(mgard_alpha);
    T bh = T(mgard_beta);
    T ch = T(mgard_gamma);
    load_v_buffer[0] = nodal[0] * ch / 2 + coeff[0] * bh + nodal[1] * ah;
    for (size_t i = 1; i < n_coeff; i++) {
        load_v_buffer[i] = (nodal[i - 1] + nodal[i + 1]) * ah + (coeff[i - 1] + coeff[i]) * bh + nodal[i] * ch;
    }
    load_v_buffer[n_coeff] = nodal[n_coeff - 1] * ah + coeff[n_coeff - 1] * bh + nodal[n_coeff] * ch / 2;
    if (n_nodal == n_coeff + 2) load_v_buffer[n_coeff + 1] = 0;
}

template <class T>
void compute_correction(T* correction_buffer, size_t n_nodal, T /*h*/, T* load_v_buffer) {
    size_t n = n_nodal;
    T* d = load_v_buffer;
    std::vector<T> b(n, T(4) / 3);
    b[0] = T(2) / 3;
    b[n - 1] = T(2) / 3;
    T c = T(1) / 3;
    for (size_t i = 1; i < n; i++) {
        T w = c / b[i - 1];
        b[i] = b[i] - w * c;
        d[i] = d[i] - w * d[i - 1];
    }
    correction_buffer[n - 1] = d[n - 1] / b[n - 1];
    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n) - 2; i >= 0; i--) {
        correction_buffer[i] = (d[i] - c * correction_buffer[i + 1]) / b[i];
    }
}

template <class T>
void compute_correction_precomputed(T* correction_buffer, size_t n_nodal, const T* w, const T* b, T /*h*/,
                                    T* load_v_buffer) {
    size_t n = n_nodal;
    T* d = load_v_buffer;
    T c = T(1) / 3;
    for (size_t i = 1; i < n; i++) {
        d[i] = d[i] - w[i] * d[i - 1];
    }
    correction_buffer[n - 1] = d[n - 1] / b[n - 1];
    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n) - 2; i >= 0; i--) {
        correction_buffer[i] = (d[i] - c * correction_buffer[i + 1]) / b[i];
    }
}

template <class T>
void compute_load_vector_vertical(T* load_v_buffer, const T* nodal_buffer, const T* coeff_buffer, size_t n1_nodal,
                                  size_t n1_coeff, size_t stride, T /*h*/, int batchsize) {
    T ah = T(mgard_alpha);
    T bh = T(mgard_beta);
    T ch = T(mgard_gamma);
    const T* nodal_pos = nodal_buffer;
    const T* coeff_pos = coeff_buffer;
    T* load_v_pos = load_v_buffer;
    for (int j = 0; j < batchsize; j++) {
        load_v_pos[j] = nodal_pos[j] * ch / 2 + coeff_pos[j] * bh + nodal_pos[stride + j] * ah;
    }
    load_v_pos += batchsize;
    nodal_pos += stride;
    coeff_pos += stride;
    for (size_t i = 1; i < n1_coeff; i++) {
        for (int j = 0; j < batchsize; j++) {
            load_v_pos[j] = (nodal_pos[j - stride] + nodal_pos[j + stride]) * ah +
                            (coeff_pos[j - stride] + coeff_pos[j]) * bh + nodal_pos[j] * ch;
        }
        load_v_pos += batchsize;
        nodal_pos += stride;
        coeff_pos += stride;
    }
    for (int j = 0; j < batchsize; j++) {
        load_v_pos[j] = nodal_pos[j - stride] * ah + coeff_pos[j - stride] * bh + nodal_pos[j] * ch / 2;
    }
    if (n1_nodal == n1_coeff + 2) {
        load_v_pos += batchsize;
        for (int j = 0; j < batchsize; j++) {
            load_v_pos[j] = 0;
        }
    }
}

template <class T>
void compute_correction_batched(T* correction_buffer, T /*h*/, const T* w, const T* b, size_t n_nodal, int batchsize,
                                size_t correction_stride, T* load_v_buffer) {
    size_t n = n_nodal;
    T c = T(1) / 3;
    T* load_v_pos = load_v_buffer + batchsize;
    for (size_t i = 1; i < n; i++) {
        for (int j = 0; j < batchsize; j++) {
            load_v_pos[j] -= w[i] * load_v_pos[-batchsize + j];
        }
        load_v_pos += batchsize;
    }
    T* correction_pos = correction_buffer + (n - 1) * correction_stride;
    load_v_pos -= batchsize;
    for (int j = 0; j < batchsize; j++) {
        correction_pos[j] = load_v_pos[j] / b[n - 1];
    }
    correction_pos -= correction_stride;
    load_v_pos -= batchsize;
    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n) - 2; i >= 0; i--) {
        for (int j = 0; j < batchsize; j++) {
            correction_pos[j] = (load_v_pos[j] - c * correction_pos[correction_stride + j]) / b[i];
        }
        correction_pos -= correction_stride;
        load_v_pos -= batchsize;
    }
}

template <class T>
void apply_correction_batched(T* nodal_pos, const T* correction_buffer, int n_nodal, int stride, int batchsize,
                              bool decompose) {
    const T* correction_pos = correction_buffer;
    if (decompose) {
        for (int i = 0; i < n_nodal; i++) {
            for (int j = 0; j < batchsize; j++) {
                nodal_pos[j] += correction_pos[j];
            }
            nodal_pos += stride;
            correction_pos += batchsize;
        }
    } else {
        for (int i = 0; i < n_nodal; i++) {
            for (int j = 0; j < batchsize; j++) {
                nodal_pos[j] -= correction_pos[j];
            }
            nodal_pos += stride;
            correction_pos += batchsize;
        }
    }
}

template <class T>
void compute_correction_vertical(T* /*data_pos*/, size_t n1, size_t n2, T h, T* horizontal_correction, size_t stride,
                                 T* load_v_buffer, const T* w, const T* b, int default_batch_size = 1) {
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    int batchsize = default_batch_size;
    int num_batches = static_cast<int>((n2_nodal - 1) / batchsize);
    T* nodal_pos = horizontal_correction;
    T* coeff_pos = horizontal_correction + n1_nodal * stride;
    for (int i = 0; i < num_batches; i++) {
        compute_load_vector_vertical(load_v_buffer, nodal_pos, coeff_pos, n1_nodal, n1_coeff, stride, h, batchsize);
        compute_correction_batched(nodal_pos, h, w, b, n1_nodal, batchsize, stride, load_v_buffer);
        nodal_pos += batchsize;
        coeff_pos += batchsize;
    }
    if (static_cast<int>(n2_nodal) - batchsize * num_batches > 0) {
        batchsize = static_cast<int>(n2_nodal) - batchsize * num_batches;
        compute_load_vector_vertical(load_v_buffer, nodal_pos, coeff_pos, n1_nodal, n1_coeff, stride, h, batchsize);
        compute_correction_batched(nodal_pos, h, w, b, n1_nodal, batchsize, stride, load_v_buffer);
    }
}

template <class T>
void compute_correction_2D(T* data_pos, T* correction_buffer, T* load_v_buffer, size_t n1, size_t n2,
                           size_t nodal_rows, T h, size_t stride, const T* w1, const T* b1, const T* w2, const T* b2,
                           int default_batch_size = 1) {
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    T* nodal_pos = data_pos;
    const T* coeff_pos = data_pos + n2_nodal;
    T* correction_pos = correction_buffer;
    for (size_t i = 0; i < n1; i++) {
        if (i < nodal_rows)
            compute_load_vector_nodal_row(load_v_buffer, n2_nodal, n2_coeff, h, coeff_pos);
        else
            compute_load_vector_coeff_row(load_v_buffer, n2_nodal, n2_coeff, h, nodal_pos, coeff_pos);
        compute_correction_precomputed(correction_pos, n2_nodal, w2, b2, h, load_v_buffer);
        nodal_pos += stride;
        coeff_pos += stride;
        correction_pos += n2_nodal;
    }
    compute_correction_vertical(data_pos, n1, n2, h, correction_buffer, n2_nodal, load_v_buffer, w1, b1,
                                default_batch_size);
}

template <class T>
void compute_correction_3D(T* data_pos, T* correction_buffer, T* load_v_buffer, size_t n1, size_t n2, size_t n3,
                           size_t /*nodal_rows_arg*/, T h, size_t dim0_stride, size_t dim1_stride,
                           int default_batch_size = 1) {
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n3_nodal = (n3 >> 1) + 1;
    std::vector<T> w1(n1_nodal);
    std::vector<T> b1(n1_nodal);
    std::vector<T> w2(n2_nodal);
    std::vector<T> b2(n2_nodal);
    std::vector<T> w3(n3_nodal);
    std::vector<T> b3(n3_nodal);
    precompute_w_and_b(w1.data(), b1.data(), n1_nodal);
    precompute_w_and_b(w2.data(), b2.data(), n2_nodal);
    precompute_w_and_b(w3.data(), b3.data(), n3_nodal);
    T* nodal_pos = data_pos;
    T* correction_pos = correction_buffer;
    for (size_t i = 0; i < n1; i++) {
        size_t local_nodal_rows = (i < n1_nodal) ? n2_nodal : 0;
        compute_correction_2D(nodal_pos, correction_pos, load_v_buffer, n2, n3, local_nodal_rows, h, dim1_stride,
                              w2.data(), b2.data(), w3.data(), b3.data(), default_batch_size);
        nodal_pos += dim0_stride;
        correction_pos += n2_nodal * n3_nodal;
    }
    correction_pos = correction_buffer;
    nodal_pos = data_pos;
    for (size_t i = 0; i < n2_nodal; i++) {
        compute_correction_vertical(data_pos, n1, n3, h, correction_pos, n2_nodal * n3_nodal, load_v_buffer,
                                    w1.data(), b1.data(), default_batch_size);
        nodal_pos += dim1_stride;
        correction_pos += n3_nodal;
    }
}

}  // namespace MGARD
}  // namespace SZ3

#endif
