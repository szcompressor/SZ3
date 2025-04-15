#ifndef _SZ_block_data_HPP
#define _SZ_block_data_HPP

#include <cstddef>
#include <memory>
#include <type_traits>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <vector>
#include <iostream>
#include <array>

namespace SZ {
    template<class T, uint N>
    class block_data : public std::enable_shared_from_this<block_data<T, N>> {
    public:

        class block_iterator {
        public:
            block_iterator(std::shared_ptr<block_data> &&mddata_, size_t block_size_) noexcept:
                    mddata(mddata_), block_size(block_size_) {
                std::fill(offset.begin(), offset.end(), 0);
            }

            ALWAYS_INLINE bool next() {
                size_t i = N - 1;
                offset[i] += block_size;
                while (i && (offset[i] >= mddata->dims[i])) {
                    offset[i] = 0;
                    offset[--i] += block_size;
                }
                return (offset[0] < mddata->dims[0]);
            }

            ALWAYS_INLINE std::array<std::pair<size_t, size_t>, N> get_block_range() const {
                std::array<std::pair<size_t, size_t>, N> range;
                for (int i = 0; i < N; i++) {
                    range[i].first = offset[i];
                    range[i].second = std::min(offset[i] + block_size, mddata->dims[i]);
                }
                return range;
            }

            /**
             *
             * @tparam args relative index in the block
             * @return
             */
            template <class ... Idx>
            ALWAYS_INLINE T *get_block_data(Idx ... args) const {
                auto ds = get_dim_strides();
                auto idx = std::vector<size_t>{static_cast<size_t>(std::forward<Idx>(args))...};
                size_t off = 0;
                for (int i = 0; i < N; i++) {
                    off += (idx[i] + offset[i]) * ds[i];
                }
                return mddata->dataptr() + off;
            }

            ALWAYS_INLINE std::array<size_t, N> get_dim_strides() const {
                return mddata->ds_padding;
            }

            template<typename Func>
            static ALWAYS_INLINE void foreach(const block_iterator &block, Func &&func) {
                auto range = block.get_block_range();
                if constexpr (N == 1) {
                    T *d = block.get_block_data(0);
                    for (size_t i = 0; i < range[0].second - range[0].first; i++) {
                        func(d++, {i});
                    }
                } else if constexpr (N == 2) {
                    for (size_t i = 0; i < range[0].second - range[0].first; i++) {
                        T *d = block.get_block_data(i, 0);
                        for (size_t j = 0; j < range[1].second - range[1].first; j++) {
                            func(d++, {i, j});
                        }
                    }
                } else if constexpr (N == 3) {
                    for (size_t i = 0; i < range[0].second - range[0].first; i++) {
                        for (size_t j = 0; j < range[1].second - range[1].first; j++) {
                            T *d = block.get_block_data(i, j, 0);
                            for (size_t k = 0; k < range[2].second - range[2].first; k++) {
                                func(d++, {i, j, k});
                            }
                        }
                    }
                } else if constexpr (N == 4) {
                    for (size_t i = 0; i < range[0].second - range[0].first; i++) {
                        for (size_t j = 0; j < range[1].second - range[1].first; j++) {
                            for (size_t k = 0; k < range[2].second - range[2].first; k++) {
                                T *d = block.get_block_data(i, j, k, 0);
                                for (size_t t = 0; t < range[3].second - range[3].first; t++) {
                                    func(d++, {i, j, k, t});
                                }
                            }
                        }
                    }
                }
            }

            friend block_data;
            std::shared_ptr<block_data> mddata;
        private:
            std::array<size_t, N> offset;
            size_t block_size;
        };


        /**Begin**/
        block_data(const T *data_,
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

        block_data(T *data_,
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
        ALWAYS_INLINE T *dataptr() {
            return data;
        }

    private:

       ALWAYS_INLINE void cal_dim_strides() {
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

        std::array<size_t, N> dims;            // dimension
        std::array<size_t, N> ds, ds_padding; // stride
        std::vector<T> internal_buffer;
        T *data;                                  // data pointer
        size_t padding;
        size_t num, num_padding;
    };

}
#endif
