#ifndef _SZ_FILE_ITERATOR_HPP
#define _SZ_FILE_ITERATOR_HPP

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
#include <fstream>

namespace SZ {
// N-dimensional file_multi_dimensional_range
    template<class T, uint N>
    class file_multi_dimensional_range : public std::enable_shared_from_this<file_multi_dimensional_range<T, N>> {
    public:

        class file_multi_dimensional_iterator {
        public:
            using value_type = T;
            using difference_type = std::ptrdiff_t;
            using reference = T &;
            using const_reference = T const &;
            using pointer = T *;
            using iterator_category = std::bidirectional_iterator_tag;

            ~file_multi_dimensional_iterator() = default;

            file_multi_dimensional_iterator() = default;

            file_multi_dimensional_iterator(file_multi_dimensional_iterator const &) = default;

            file_multi_dimensional_iterator &operator=(file_multi_dimensional_iterator const &) = default;

            file_multi_dimensional_iterator(file_multi_dimensional_iterator &&) noexcept = default;

            file_multi_dimensional_iterator &operator=(file_multi_dimensional_iterator &&) noexcept = default;

            file_multi_dimensional_iterator(std::shared_ptr<file_multi_dimensional_range> &&range_, std::size_t current_offset_) noexcept:
                    range(range_), global_offset(current_offset_), local_index{} {

            }

            inline file_multi_dimensional_iterator &operator++() {
                size_t i = N - 1;
                local_index[i]++;
                ptrdiff_t offset = range->global_dim_strides[i];
                while (i && (local_index[i] == range->dimensions[i])) {
                    offset -= range->dimensions[i] * range->global_dim_strides[i];
                    local_index[i--] = 0;
                    offset += range->global_dim_strides[i];
                    local_index[i]++;
                }
                global_offset += offset;

                if (global_offset >= range->buffer_offset + range->buffer_size || global_offset < range->buffer_offset) {
                    range->read_buffer(global_offset);
                }

                // std::cout << "offset=" << offset << ", current_offset=" << current_offset << std::endl;
                return *this;
            }

            file_multi_dimensional_iterator operator++(int) {
                auto cpy = *this;
                ++(*this);
                return cpy;
            }

            pointer operator->() {
                return range->buffer[global_offset - range->buffer_offset];
            }

            pointer operator->() const {
                return range->buffer[global_offset - range->buffer_offset];
            }

            reference operator*() {
                return range->buffer[global_offset - range->buffer_offset];
            }

            const_reference operator*() const {
                return range->buffer[global_offset - range->buffer_offset];
            }

            bool operator==(file_multi_dimensional_iterator const &rhs) const {
                return global_offset == rhs.global_offset;
            }

            bool operator!=(file_multi_dimensional_iterator const &rhs) const {
                return global_offset != rhs.global_offset;
            }


            std::array<size_t, N> get_global_index() const {
                auto offset = global_offset;
                std::array<size_t, N> global_idx{0};
                for (int i = N - 1; i >= 0; i--) {
                    global_idx[i] = offset % range->global_dimensions[i];
                    offset /= range->global_dimensions[i];
                }
                return global_idx;
            }

            std::array<size_t, N> get_local_index() const {
                return local_index;
            }

            size_t get_local_index(size_t i) const {
                return local_index[i];
            }

            ptrdiff_t get_offset() const {
                return global_offset;
            }

            std::array<size_t, N> get_dimensions() const {
                return range->get_dimensions();
            }


            void print() {
                std::cout << "global_offset=" << global_offset << ", local_index=(";
                for (auto const &i:local_index) {
                    std::cout << i << ",";
                }
                std::cout << "),local_dim=[";
                for (auto const &i:range->dimensions) {
                    std::cout << i << ",";
                }
                std::cout << "]" << " ";
                std::cout << "global_index=(";
                for (auto const &i:get_global_index()) {
                    std::cout << i << ",";
                }
                std::cout << "),global_dim=[";
                for (auto const &i:range->global_dimensions) {
                    std::cout << i << ",";
                }
                std::cout << "]" << std::endl;
            }


        private:
            friend file_multi_dimensional_range;
            std::shared_ptr<file_multi_dimensional_range> range;
            std::array<size_t, N> local_index;        // index of current_offset position
            ptrdiff_t global_offset;
        };

        using iterator = file_multi_dimensional_iterator;
        using const_iterator = file_multi_dimensional_iterator;
        using value_type = T;
        using reference = T &;
        using pointer = T *;

        file_multi_dimensional_iterator begin() {
            return file_multi_dimensional_iterator(this->shared_from_this(), start_offset);
        }

        file_multi_dimensional_iterator end() {
            return file_multi_dimensional_iterator(this->shared_from_this(), end_offset);
        }

        template<class ForwardIt1>
        void set_dimensions(ForwardIt1 begin, ForwardIt1 end) {
            int i = 0;
            for (auto iter = begin; iter != end; ++iter) {
                dimensions[i++] = *iter;
                // std::cout << dimensions[i-1] << " ";
            }
        }

        void set_dimensions_auto() {
            // std::cout << "dimensions: ";
            for (int i = 0; i < dimensions.size(); i++) {
                // std::cout << "g[i]=" << global_dimensions[i] << ",str=" << access_stride << " ";
                dimensions[i] = (global_dimensions[i] - 1) / access_stride + 1;
                // std::cout << dimensions[i] << " ";
            }
            // std::cout << std::endl;
        }

        void set_global_dim_strides() {
            // std::cout << "strides: ";
            size_t cur_stride = 1;
            for (int i = N - 1; i >= 0; i--) {
                global_dim_strides[i] = cur_stride * access_stride;
                cur_stride *= global_dimensions[i];
                // std::cout << dim_strides[i] << " ";
            }
            // std::cout << std::endl;
        }

        void set_offsets(ptrdiff_t offset_) {
            start_offset = offset_;
            end_offset = start_offset + dimensions[0] * global_dim_strides[0];
            if (start_offset >= buffer_offset + buffer_size || start_offset < buffer_offset) {
                read_buffer(start_offset);
            }
        }

        void set_access_stride(size_t stride_) {
            access_stride = stride_;
        }

        // NOTE: did not consider the real offset for simplicity
        void set_starting_position(const std::array<size_t, N> &dims) {
            for (int i = 0; i < N; i++) {
                start_position[i] = (dims[i] == 0);
            }
        }

        ~file_multi_dimensional_range() {
            fin.close();
        }

        template<class ForwardIt1>
        file_multi_dimensional_range(
                const char *file,
                ForwardIt1 global_dims_begin,
                ForwardIt1 global_dims_end,
                size_t stride_,
                ptrdiff_t offset_
        ): start_position{false} {
            static_assert(
                    std::is_convertible<
                            typename std::iterator_traits<ForwardIt1>::value_type,
                            std::size_t>::value,
                    "ForwardIt1 must be convertible to std::size_t"
            );
            if (global_dims_end - global_dims_begin != N) {
                std::cout << global_dims_end - global_dims_begin << " " << N << std::endl;
                std::cerr << "#dimensions does not match!\n";
                exit(0);
            }


            set_access_stride(stride_);
            // set global dimensions
            int i = 0;
            for (auto iter = global_dims_begin; iter != global_dims_end; ++iter) {
                global_dimensions[i++] = *iter;
            }
//            size_t cur_stride = stride_;
//            for (int i = N - 1; i >= 0; i--) {
//                global_dim_strides[i] = cur_stride;
//                cur_stride *= global_dimensions[i];
//            }
            // set_dimensions(dims_begin, dims_end);
            set_dimensions_auto();
            set_global_dim_strides();
            set_offsets(offset_);

            data_size = std::accumulate(global_dimensions.begin(), global_dimensions.end(), 1, std::multiplies<size_t>());
            buffer.resize(buffer_size);
            fin = std::ifstream(file, std::ios::binary);
            if (!fin) {
                std::cout << " Error, Couldn't find the file" << "\n";
                exit(0);
            }
            read_buffer(0);
        }

        size_t get_dimensions(size_t i) const {
            return dimensions[i];
        }

        std::array<size_t, N> get_dimensions() const {
            return dimensions;
        }

        std::array<size_t, N> get_global_dimensions() const {
            return global_dimensions;
        }

        bool whether_global_start_position(size_t i) const {
            return start_position[i];
        }

        void read_buffer(size_t offset) {
            if (offset >= data_size) {
                return;
            }
            fin.seekg(offset * sizeof(T), std::ios::beg);
            size_t size = data_size - offset > buffer_size ? buffer_size : data_size - offset;
            fin.read(reinterpret_cast<char *>(buffer.data()), size * sizeof(T));
            buffer_offset = offset;
//            printf("buffer %lu %lu\n", buffer_offset, size);
        }

    private:
        std::array<size_t, N> global_dimensions;
        std::array<size_t, N> global_dim_strides;
        std::array<size_t, N> dimensions;              // the dimensions
//        std::array<size_t, N> dim_strides;              // strides for dimensions
        std::array<bool, N> start_position;       // indicator for starting position, used for block-wise lorenzo predictor
        size_t access_stride;                                // stride for access pattern
        ptrdiff_t start_offset;                              // offset for start point
        ptrdiff_t end_offset;                                  // offset for end point
        std::ifstream fin;                                     // file pointer
        size_t buffer_offset = 0, buffer_size = 1000 * 1000, data_size;
        std::vector<T> buffer;

    };

}
#endif
