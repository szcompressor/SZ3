#ifndef _SZ_MULTI_DIMENSIONAL_DATA_HPP
#define _SZ_MULTI_DIMENSIONAL_DATA_HPP

#include <cassert>
#include <cstddef>
#include <memory>
#include <iterator>
#include <type_traits>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <vector>
#include <iostream>
#include <array>
#include <typeinfo>

namespace SZ {
    template<class T, uint N>
    class multi_dimensional_data : public std::enable_shared_from_this<multi_dimensional_data<T, N>> {
    public:

        class block_iterator {
        public:
            block_iterator(std::shared_ptr<multi_dimensional_data> &&mddata_, size_t block_size_) noexcept:
                    mddata(mddata_), block_size(block_size_) {
                std::fill(l.begin(), l.end(), 0);
            }

            inline bool next() {
                size_t i = N - 1;
                l[i] += block_size;
                while (i && (l[i] >= mddata->dims[i])) {
                    l[i] = 0;
                    l[--i] += block_size;
                }
                return (l[0] < mddata->dims[0]);
            }

            inline std::array<std::pair<size_t, size_t>, N> get_block_range() const {
                std::array<std::pair<size_t, size_t>, N> range;
                for (int i = 0; i < N; i++) {
                    range[i].first = l[i];
                    range[i].second = std::min(l[i] + block_size, mddata->dims[i]);
                }
                return range;
            }

            template<class ... Idx>
            ALWAYS_INLINE T *get_data(Idx ... args) const {
                auto ds = get_dim_strides();
                auto idx = std::vector<size_t>{static_cast<size_t>(std::forward<Idx>(args))...};
                size_t offset = 0;
                for (int i = 0; i < N; i++) {
                    offset += (idx[i] + l[i]) * ds[i];
                }
                return mddata->dataptr() + offset;
            }

            inline std::array<size_t, N> get_dim_strides() const {
                return mddata->get_dim_strides();
            }

            template <typename Func>
            void iterate_block(Func &&func) const {
                auto range = this->get_block_range();
                auto d = this->mddata;

                if constexpr (N == 1) {
                    for (size_t i = range[0].first; i < range[0].second; i++) {
                        T *c = d->get_data(i);
                        func(c);
                    }
                } else if constexpr (N == 2) {
                    for (size_t i = range[0].first; i < range[0].second; i++) {
                        for (size_t j = range[1].first; j < range[1].second; j++) {
                            T *c = d->get_data(i, j);
                            func(c);
                        }
                    }
                } else if constexpr (N == 3) {
                    for (size_t i = range[0].first; i < range[0].second; i++) {
                        for (size_t j = range[1].first; j < range[1].second; j++) {
                            for (size_t k = range[2].first; k < range[2].second; k++) {
                                T *c = d->get_data(i, j, k);
                                func(c);
                            }
                        }
                    }
                } else if constexpr (N == 4) {
                    for (size_t i = range[0].first; i < range[0].second; i++) {
                        for (size_t j = range[1].first; j < range[1].second; j++) {
                            for (size_t k = range[2].first; k < range[2].second; k++) {
                                for (size_t t = range[3].first; t < range[3].second; t++) {
                                    T *c = d->get_data(i, j, k, t);
                                    func(c);
                                    // func(block, i, j, k, t);
                                }
                            }
                        }
                    }
                }
            }

            friend multi_dimensional_data;
            std::shared_ptr<multi_dimensional_data> mddata;
        private:
            std::array<size_t, N> l;
            size_t block_size;
        };


        /**Begin**/
        multi_dimensional_data(const T *data_,
                               const std::vector<size_t> &dims_,
                               size_t padding_ = 0
        ) : padding(padding_) {
            if (dims_.size() != N) {
                throw std::invalid_argument("#dims does not match!");
            }
            std::copy_n(dims_.begin(), N, dims.begin());
            cal_dim_strides();

            if (padding > 0) {
                internal_buffer.resize(num_padding);
                size_t offset = std::accumulate(ds_padding.begin(), ds_padding.end(), (size_t)(0));
                data = &internal_buffer[padding * offset];
                copy_data_in(data_);
            } else {
                internal_buffer.resize(num);
                data = internal_buffer.data();
                copy_data_in(data_);
            }
        }

        multi_dimensional_data(T *data_,
                               const std::vector<size_t> &dims_,
                               size_t padding_ = 0
        ) : padding(padding_) {
            if (dims_.size() != N) {
                throw std::invalid_argument("#dims does not match!");
            }
            std::copy_n(dims_.begin(), N, dims.begin());
            cal_dim_strides();

            if (padding > 0) {
                internal_buffer.resize(num_padding);
                size_t offset = std::accumulate(ds_padding.begin(), ds_padding.end(), (size_t)(0));
                data = &internal_buffer[padding * offset];
            } else {
                data = data_;
            }
        }

        block_iterator block_iter(size_t block_size) {
            return block_iterator(this->shared_from_this(), block_size);
        }

        inline void cal_dim_strides() {
            size_t cur_stride = 1, cur_stride_pading = 1;
            for (int i = N - 1; i >= 0; i--) {
                ds[i] = cur_stride;
                ds_padding[i] = cur_stride_pading;
                cur_stride *= dims[i];
                cur_stride_pading *= (dims[i] + padding);
            }
            num = cur_stride;
            num_padding = cur_stride_pading;
        }

        inline std::array<size_t, N> get_dims() const {
            return dims;
        }

        inline std::array<size_t, N> get_dim_strides() const {
            return ds_padding;
        }

        template<class ... Idx>
        inline T *get_data(Idx ... args) {
            auto idx = std::vector<size_t>{static_cast<size_t>(std::forward<Idx>(args))...};
            size_t offset = 0;
            for (int i = 0; i < N; i++) {
                offset += idx[i] * ds_padding[i];
            }
            return dataptr() + offset;
        }


        void copy_data_in(const T *input) {
            if (padding == 0) {
                std::copy(input, input + num, internal_buffer.begin());
                return;
            }
            if (N == 1) {
                memcpy(&internal_buffer[padding], &input[0], dims[0] * sizeof(T));
            } else if (N == 2) {
                for (size_t i = 0; i < dims[0]; i++) {
                    memcpy(&internal_buffer[(i + padding) * ds_padding[0] + padding], &input[i * ds[0]], dims[1] * sizeof(T));
                }
            } else if (N == 3) {
                for (size_t i = 0; i < dims[0]; i++) {
                    for (size_t j = 0; j < dims[1]; j++) {
                        memcpy(&internal_buffer[(i + padding) * ds_padding[0] + (j + padding) * ds_padding[1] + padding],
                               &input[i * ds[0] + j * ds[1]], dims[2] * sizeof(T));
                    }
                }
            } else if (N == 4) {
                for (size_t i = 0; i < dims[0]; i++) {
                    for (size_t j = 0; j < dims[1]; j++) {
                        for (size_t k = 0; k < dims[2]; k++) {
                            memcpy(&internal_buffer[(i + padding) * ds_padding[0] + (j + padding) * ds_padding[1] + (k + padding) * ds_padding[2] +
                                                    padding],
                                   &input[i * ds[0] + j * ds[1] + k * ds[2]], dims[3] * sizeof(T));
                        }
                    }
                }
            } else {
                throw std::invalid_argument("N (dimension) should be less than 5");
            }
        }

        void copy_data_out(T *output) {
            if (N == 1) {
                memcpy(&output[0], &data[0], dims[0] * sizeof(T));
            } else if (N == 2) {
                for (size_t i = 0; i < dims[0]; i++) {
                    memcpy(&output[i * ds[0]], &data[i * ds_padding[0]], dims[1] * sizeof(T));
                }
            } else if (N == 3) {
                for (size_t i = 0; i < dims[0]; i++) {
                    for (size_t j = 0; j < dims[1]; j++) {
                        memcpy(&output[i * ds[0] + j * ds[1]], &data[i * ds_padding[0] + j * ds_padding[1]], dims[2] * sizeof(T));
                    }
                }
            } else if (N == 4) {
                for (size_t i = 0; i < dims[0]; i++) {
                    for (size_t j = 0; j < dims[1]; j++) {
                        for (size_t k = 0; k < dims[2]; k++) {
                            memcpy(&output[i * ds[0] + j * ds[1] + k * ds[2]], &data[i * ds_padding[0] + j * ds_padding[1] + k * ds_padding[2]],
                                   dims[3] * sizeof(T));
                        }
                    }
                }
            } else {
                throw std::invalid_argument("N (dimension) should be less than 5");
            }
        }

    protected:
        inline T *dataptr() {
            return data;
        }

    private:
        std::array<size_t, N> dims;            // dimension
        std::array<size_t, N> ds, ds_padding; // stride
        std::vector<T> internal_buffer;
        T *data;                                  // data pointer
        size_t padding;
        size_t num, num_padding;
    };

}
#endif
