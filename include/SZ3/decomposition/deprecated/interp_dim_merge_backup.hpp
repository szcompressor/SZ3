    double block_interpolation_1d_fastest_dim_first_3d(
        T *data,
        const std::array<size_t, N> &begin_idx,
        const std::array<size_t, N> &end_idx,
        const size_t &direction,
        std::array<size_t, N> &steps,
        const size_t &math_stride,
        const std::string &interp_func,
        const PredictorBehavior pb) { 
        for (size_t i = 0; i < N; i++) {
            if (end_idx[i] < begin_idx[i])
                return 0;
        }
        size_t math_begin_idx = begin_idx[direction], math_end_idx = end_idx[direction];
        size_t n = (math_end_idx - math_begin_idx) / math_stride + 1;
        if (n <= 1) {
            return 0;
        }
        double predict_error = 0.0;
        size_t begin = 0, global_end_idx = global_dimensions[direction];
        for (size_t i = 0; i < N; i++)
            begin += dimension_offsets[i] * begin_idx[i];
        size_t stride = math_stride * dimension_offsets[direction];
        std::array<size_t, N> begins, ends, strides;
        for (size_t i = 0; i < N; i++) {
            begins[i] = 0;
            ends[i] = end_idx[i] - begin_idx[i] + 1;
            strides[i] = dimension_offsets[i];
        }
        strides[direction] = stride;
        size_t stride2x = 2 * stride;
        if (pb == PB_predict_overwrite) {
            if (interp_func == "linear") {
                begins[direction] = 1;
                ends[direction] = n - 1;
                steps[direction] = 2;
                for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                    for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                        for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                            T *d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                            quantize(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                        }
                    }
                }
                if (n % 2 == 0) {
                    begins[direction] = n - 1;
                    ends[direction] = n;
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                T *d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                if (math_end_idx + math_stride < global_end_idx)
                                    quantize(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                                else if (n < 3)
                                    quantize(d - data, *d, *(d - stride));
                                else
                                    quantize(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)));
                            }
                        }
                    }
                }
            } else if (interp_func == "cubic") {
                size_t stride3x = 3 * stride;
                T *d;
                size_t i_start = 3;
                begins[direction] = i_start;
                ends[direction] = (n >= 3) ? (n - 3) : 0;
                steps[direction] = 2;
                for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                    for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                        for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                            d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                            quantize(d - data, *d,
                                     interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                        }
                    }
                }
                std::vector<size_t> boundary;
                boundary.push_back(1);
                if (n % 2 == 1) {
                    if (n > 3)
                        boundary.push_back(n - 2);
                } else {
                    if (n > 4)
                        boundary.push_back(n - 3);
                    if (n > 2)
                        boundary.push_back(n - 1);
                }
                for (auto ii : boundary) {
                    begins[direction] = ii;
                    ends[direction] = ii + 1;
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                if (ii >= 3) {
                                    if (ii + 3 < n)
                                        quantize(d - data, *d,
                                                 interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                    else if (ii + 1 < n)
                                        quantize(d - data, *d,
                                                 interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                    else
                                        quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                                } else {
                                    if (ii + 3 < n)
                                        quantize(d - data, *d,
                                                 interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                    else if (ii + 1 < n)
                                        quantize(d - data, *d,
                                                 interp_linear(*(d - stride), *(d + stride)));
                                    else
                                        quantize(d - data, *d, *(d - stride));
                                }
                            }
                        }
                    }
                }
            } else {
                size_t stride3x = 3 * stride;
                T *d;
                size_t i_start = 3;
                begins[direction] = i_start;
                ends[direction] = (n >= 3) ? (n - 3) : 0;
                steps[direction] = 2;
                for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                    for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                        for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                            d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                            quantize(d - data, *d,
                                     interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                        }
                    }
                }
                std::vector<size_t> boundary;
                boundary.push_back(1);
                if (n % 2 == 1) {
                    if (n > 3)
                        boundary.push_back(n - 2);
                } else {
                    if (n > 4)
                        boundary.push_back(n - 3);
                    if (n > 2)
                        boundary.push_back(n - 1);
                }
                for (auto ii : boundary) {
                    begins[direction] = ii;
                    ends[direction] = ii + 1;
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                if (ii >= 3) {
                                    if (ii + 3 < n)
                                        quantize(d - data, *d,
                                                 interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                    else if (ii + 1 < n)
                                        quantize(d - data, *d,
                                                 interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                    else
                                        quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                                } else {
                                    if (ii + 3 < n)
                                        quantize(d - data, *d,
                                                 interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                    else if (ii + 1 < n)
                                        quantize(d - data, *d,
                                                 interp_linear(*(d - stride), *(d + stride)));
                                    else
                                        quantize(d - data, *d, *(d - stride));
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (interp_func == "linear") {
                begins[direction] = 1;
                ends[direction] = n - 1;
                steps[direction] = 2;
                for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                    for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                        for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                            T *d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                            recover(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                        }
                    }
                }
                if (n % 2 == 0) {
                    begins[direction] = n - 1;
                    ends[direction] = n;
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                T *d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                if (math_end_idx + math_stride < global_end_idx)
                                    recover(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                                else if (n < 3)
                                    recover(d - data, *d, *(d - stride));
                                else
                                    recover(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)));
                            }
                        }
                    }
                }
            } else if (interp_func == "cubic") {
                size_t stride3x = 3 * stride;
                T *d;
                size_t i_start = 3;
                begins[direction] = i_start;
                ends[direction] = (n >= 3) ? (n - 3) : 0;
                steps[direction] = 2;
                for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                    for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                        for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                            d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                            recover(d - data, *d,
                                    interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                        }
                    }
                }
                std::vector<size_t> boundary;
                boundary.push_back(1);
                if (n % 2 == 1) {
                    if (n > 3)
                        boundary.push_back(n - 2);
                } else {
                    if (n > 4)
                        boundary.push_back(n - 3);
                    if (n > 2)
                        boundary.push_back(n - 1);
                }
                for (auto ii : boundary) {
                    begins[direction] = ii;
                    ends[direction] = ii + 1;
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                if (ii >= 3) {
                                    if (ii + 3 < n)
                                        recover(d - data, *d,
                                                interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                    else if (ii + 1 < n)
                                        recover(d - data, *d,
                                                interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                    else
                                        recover(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                                } else {
                                    if (ii + 3 < n)
                                        recover(d - data, *d,
                                                interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                    else if (ii + 1 < n)
                                        recover(d - data, *d,
                                                interp_linear(*(d - stride), *(d + stride)));
                                    else
                                        recover(d - data, *d, *(d - stride));
                                }
                            }
                        }
                    }
                }
            } else {
                size_t stride3x = 3 * stride;
                T *d;
                size_t i_start = 3;
                begins[direction] = i_start;
                ends[direction] = (n >= 3) ? (n - 3) : 0;
                steps[direction] = 2;
                for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                    for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                        for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                            d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                            recover(d - data, *d,
                                    interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                        }
                    }
                }
                std::vector<size_t> boundary;
                boundary.push_back(1);
                if (n % 2 == 1) {
                    if (n > 3)
                        boundary.push_back(n - 2);
                } else {
                    if (n > 4)
                        boundary.push_back(n - 3);
                    if (n > 2)
                        boundary.push_back(n - 1);
                }
                for (auto ii : boundary) {
                    begins[direction] = ii;
                    ends[direction] = ii + 1;
                    for (size_t i = begins[0]; i < ends[0]; i += steps[0]) {
                        for (size_t j = begins[1]; j < ends[1]; j += steps[1]) {
                            for (size_t k = begins[2]; k < ends[2]; k += steps[2]) {
                                d = data + begin + i * strides[0] + j * strides[1] + k * strides[2];
                                if (ii >= 3) {
                                    if (ii + 3 < n)
                                        recover(d - data, *d,
                                                interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                                    else if (ii + 1 < n)
                                        recover(d - data, *d,
                                                interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                                    else
                                        recover(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                                } else {
                                    if (ii + 3 < n)
                                        recover(d - data, *d,
                                                interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                                    else if (ii + 1 < n)
                                        recover(d - data, *d,
                                                interp_linear(*(d - stride), *(d + stride)));
                                    else
                                        recover(d - data, *d, *(d - stride));
                                }
                            }
                        }
                    }
                }
            }
        }
        return predict_error;
    }
    




    double block_interpolation_1d_fastest_dim_first(
        T *data,
        const std::array<size_t, N> &begin_idx,
        const std::array<size_t, N> &end_idx,
        const size_t &direction,
        std::array<size_t, N> &steps,
        const size_t &math_stride,
        const std::string &interp_func,
        const PredictorBehavior pb) {
        for (size_t i = 0; i < N; i++) {
            if (end_idx[i] < begin_idx[i])
                return 0;
        }
        size_t math_begin_idx = begin_idx[direction], math_end_idx = end_idx[direction];
        size_t n = (math_end_idx - math_begin_idx) / math_stride + 1;
        if (n <= 1)
            return 0;
        double predict_error = 0.0;
        size_t begin = 0;
        for (size_t i = 0; i < N; i++)
            begin += dimension_offsets[i] * begin_idx[i];

        size_t stride = math_stride * dimension_offsets[direction];
        std::array<size_t, N> begins, ends, strides;
        for (size_t i = 0; i < N; i++) {
            begins[i] = 0;
            ends[i] = end_idx[i] - begin_idx[i] + 1;
            strides[i] = dimension_offsets[i];
        }
        strides[direction] = stride;
        size_t stride2x = 2 * stride;

        auto for_each_index = [&](auto &self, size_t dim, size_t offset,
                                    const std::array<size_t, N> &begins,
                                    const std::array<size_t, N> &ends,
                                    const std::array<size_t, N> &strides,
                                    const std::array<size_t, N> &steps,
                                    auto f) -> void {
            if (dim == N) {
                f(offset);
                return;
            }
            for (size_t idx = begins[dim]; idx < ends[dim]; idx += steps[dim])
                self(self, dim + 1, offset + idx * strides[dim], begins, ends, strides, steps, f);
        };

        // 根据预测行为（量化或恢复）以及插值类型分别处理
        if (pb == PB_predict_overwrite) {
            if (interp_func == "linear") {
                begins[direction] = 1;
                ends[direction] = n - 1;
                steps[direction] = 2;
                // 主循环：调用 quantize 和 interp_linear
                for_each_index(for_each_index, 0, begin, begins, ends, strides, steps,
                               [&](size_t offset) {
                                   T *d = data + offset;
                                   quantize(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                               });
                if (n % 2 == 0) {
                    begins[direction] = n - 1;
                    ends[direction] = n;
                    for_each_index(for_each_index, 0, begin, begins, ends, strides, steps,
                                   [&](size_t offset) {
                                       T *d = data + offset;
                                       if (n < 3)
                                           quantize(d - data, *d, *(d - stride));
                                       else
                                           quantize(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)));
                                   });
                }
            } else if (interp_func == "cubic") {
                size_t stride3x = 3 * stride;
                size_t i_start = 3;
                begins[direction] = i_start;
                ends[direction] = (n >= 3) ? (n - 3) : 0;
                steps[direction] = 2;
                for_each_index(for_each_index, 0, begin, begins, ends, strides, steps,
                               [&](size_t offset) {
                                   T *d = data + offset;
                                   quantize(d - data, *d,
                                            interp_cubic(*(d - stride3x), *(d - stride), *(d + stride),
                                                         *(d + stride3x)));
                               });
                // 边界处理
                std::vector<size_t> boundary;
                boundary.push_back(1);
                if (n % 2 == 1) {
                    if (n > 3)
                        boundary.push_back(n - 2);
                } else {
                    if (n > 4)
                        boundary.push_back(n - 3);
                    if (n > 2)
                        boundary.push_back(n - 1);
                }
                for (auto ii : boundary) {
                    begins[direction] = ii;
                    ends[direction] = ii + 1;
                    for_each_index(for_each_index, 0, begin, begins, ends, strides, steps,
                                   [&](size_t offset) {
                                       T *d = data + offset;
                                       if (ii >= 3) {
                                           if (ii + 3 < n)
                                               quantize(d - data, *d,
                                                        interp_cubic(*(d - stride3x), *(d - stride),
                                                                     *(d + stride), *(d + stride3x)));
                                           else if (ii + 1 < n)
                                               quantize(d - data, *d,
                                                        interp_quad_2(*(d - stride3x), *(d - stride),
                                                                      *(d + stride)));
                                           else
                                               quantize(d - data, *d,
                                                        interp_linear1(*(d - stride3x), *(d - stride)));
                                       } else {
                                           if (ii + 3 < n)
                                               quantize(d - data, *d,
                                                        interp_quad_1(*(d - stride), *(d + stride),
                                                                      *(d + stride3x)));
                                           else if (ii + 1 < n)
                                               quantize(d - data, *d,
                                                        interp_linear(*(d - stride), *(d + stride)));
                                           else
                                               quantize(d - data, *d, *(d - stride));
                                       }
                                   });
                }
            } else { 
                size_t stride3x = 3 * stride;
                size_t i_start = 3;
                begins[direction] = i_start;
                ends[direction] = (n >= 3) ? (n - 3) : 0;
                steps[direction] = 2;
                for_each_index(for_each_index, 0, begin, begins, ends, strides, steps,
                               [&](size_t offset) {
                                   T *d = data + offset;
                                   quantize(d - data, *d,
                                            interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride),
                                                         *(d + stride3x)));
                               });
                // 边界处理
                std::vector<size_t> boundary;
                boundary.push_back(1);
                if (n % 2 == 1) {
                    if (n > 3)
                        boundary.push_back(n - 2);
                } else {
                    if (n > 4)
                        boundary.push_back(n - 3);
                    if (n > 2)
                        boundary.push_back(n - 1);
                }
                for (auto ii : boundary) {
                    begins[direction] = ii;
                    ends[direction] = ii + 1;
                    for_each_index(for_each_index, 0, begin, begins, ends, strides, steps,
                                   [&](size_t offset) {
                                       T *d = data + offset;
                                       if (ii >= 3) {
                                           if (ii + 3 < n)
                                               quantize(d - data, *d,
                                                        interp_cubic_natural(*(d - stride3x), *(d - stride),
                                                                     *(d + stride), *(d + stride3x)));
                                           else if (ii + 1 < n)
                                               quantize(d - data, *d,
                                                        interp_quad_2(*(d - stride3x), *(d - stride),
                                                                      *(d + stride)));
                                           else
                                               quantize(d - data, *d,
                                                        interp_linear1(*(d - stride3x), *(d - stride)));
                                       } else {
                                           if (ii + 3 < n)
                                               quantize(d - data, *d,
                                                        interp_quad_1(*(d - stride), *(d + stride),
                                                                      *(d + stride3x)));
                                           else if (ii + 1 < n)
                                               quantize(d - data, *d,
                                                        interp_linear(*(d - stride), *(d + stride)));
                                           else
                                               quantize(d - data, *d, *(d - stride));
                                       }
                                   });
                }
            }
        } else {
            // pb 非 PB_predict_overwrite，调用 recover 函数
            if (interp_func == "linear") {
                begins[direction] = 1;
                ends[direction] = n - 1;
                steps[direction] = 2;
                for_each_index(for_each_index, 0, begin, begins, ends, strides, steps,
                               [&](size_t offset) {
                                   T *d = data + offset;
                                   recover(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                               });
                if (n % 2 == 0) {
                    begins[direction] = n - 1;
                    ends[direction] = n;
                    for_each_index(for_each_index, 0, begin, begins, ends, strides, steps,
                                   [&](size_t offset) {
                                       T *d = data + offset;
                                       if (n < 3)
                                           recover(d - data, *d, *(d - stride));
                                       else
                                           recover(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)));
                                   });
                }
            } else if (interp_func == "cubic") {
                size_t stride3x = 3 * stride;
                size_t i_start = 3;
                begins[direction] = i_start;
                ends[direction] = (n >= 3) ? (n - 3) : 0;
                steps[direction] = 2;
                for_each_index(for_each_index, 0, begin, begins, ends, strides, steps,
                               [&](size_t offset) {
                                   T *d = data + offset;
                                   recover(d - data, *d,
                                           interp_cubic(*(d - stride3x), *(d - stride), *(d + stride),
                                                        *(d + stride3x)));
                               });
                std::vector<size_t> boundary;
                boundary.push_back(1);
                if (n % 2 == 1) {
                    if (n > 3)
                        boundary.push_back(n - 2);
                } else {
                    if (n > 4)
                        boundary.push_back(n - 3);
                    if (n > 2)
                        boundary.push_back(n - 1);
                }
                for (auto ii : boundary) {
                    begins[direction] = ii;
                    ends[direction] = ii + 1;
                    for_each_index(for_each_index, 0, begin, begins, ends, strides, steps,
                                   [&](size_t offset) {
                                       T *d = data + offset;
                                       if (ii >= 3) {
                                           if (ii + 3 < n)
                                               recover(d - data, *d,
                                                       interp_cubic(*(d - stride3x), *(d - stride),
                                                                    *(d + stride), *(d + stride3x)));
                                           else if (ii + 1 < n)
                                               recover(d - data, *d,
                                                       interp_quad_2(*(d - stride3x), *(d - stride),
                                                                     *(d + stride)));
                                           else
                                               recover(d - data, *d,
                                                       interp_linear1(*(d - stride3x), *(d - stride)));
                                       } else {
                                           if (ii + 3 < n)
                                               recover(d - data, *d,
                                                       interp_quad_1(*(d - stride), *(d + stride),
                                                                     *(d + stride3x)));
                                           else if (ii + 1 < n)
                                               recover(d - data, *d,
                                                       interp_linear(*(d - stride), *(d + stride)));
                                           else
                                               recover(d - data, *d, *(d - stride));
                                       }
                                   });
                }
            } else {
                size_t stride3x = 3 * stride;
                size_t i_start = 3;
                begins[direction] = i_start;
                ends[direction] = (n >= 3) ? (n - 3) : 0;
                steps[direction] = 2;
                for_each_index(for_each_index, 0, begin, begins, ends, strides, steps,
                               [&](size_t offset) {
                                   T *d = data + offset;
                                   recover(d - data, *d,
                                           interp_cubic_natural(*(d - stride3x), *(d - stride), *(d + stride),
                                                        *(d + stride3x)));
                               });
                std::vector<size_t> boundary;
                boundary.push_back(1);
                if (n % 2 == 1) {
                    if (n > 3)
                        boundary.push_back(n - 2);
                } else {
                    if (n > 4)
                        boundary.push_back(n - 3);
                    if (n > 2)
                        boundary.push_back(n - 1);
                }
                for (auto ii : boundary) {
                    begins[direction] = ii;
                    ends[direction] = ii + 1;
                    for_each_index(for_each_index, 0, begin, begins, ends, strides, steps,
                                   [&](size_t offset) {
                                       T *d = data + offset;
                                       if (ii >= 3) {
                                           if (ii + 3 < n)
                                               recover(d - data, *d,
                                                       interp_cubic_natural(*(d - stride3x), *(d - stride),
                                                                    *(d + stride), *(d + stride3x)));
                                           else if (ii + 1 < n)
                                               recover(d - data, *d,
                                                       interp_quad_2(*(d - stride3x), *(d - stride),
                                                                     *(d + stride)));
                                           else
                                               recover(d - data, *d,
                                                       interp_linear1(*(d - stride3x), *(d - stride)));
                                       } else {
                                           if (ii + 3 < n)
                                               recover(d - data, *d,
                                                       interp_quad_1(*(d - stride), *(d + stride),
                                                                     *(d + stride3x)));
                                           else if (ii + 1 < n)
                                               recover(d - data, *d,
                                                       interp_linear(*(d - stride), *(d + stride)));
                                           else
                                               recover(d - data, *d, *(d - stride));
                                       }
                                   });
                }
            }

        }
        return predict_error;
    }