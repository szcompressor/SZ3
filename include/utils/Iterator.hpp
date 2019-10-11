#ifndef _SZ_ITERATOR_HPP
#define _SZ_ITERATOR_HPP

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

namespace SZ{
// N-dimensional multi_dimensional_range
template <class T, size_t N>
class multi_dimensional_range: public std::enable_shared_from_this<multi_dimensional_range<T, N>> {
  public:

  class multi_dimensional_iterator {
    public:
      using value_type = T;
      using difference_type = std::ptrdiff_t;
      using reference = T&;
      using const_reference = T const&;
      using pointer = T*;
      using iterator_category = std::bidirectional_iterator_tag;

      ~multi_dimensional_iterator()=default;
      multi_dimensional_iterator()=default;
      multi_dimensional_iterator(multi_dimensional_iterator const&)=default;
      multi_dimensional_iterator& operator=(multi_dimensional_iterator const&)=default;
      multi_dimensional_iterator(multi_dimensional_iterator &&) noexcept=default;
      multi_dimensional_iterator& operator=(multi_dimensional_iterator &&) noexcept=default;
      multi_dimensional_iterator(std::shared_ptr<multi_dimensional_range>&& range_, std::size_t current_offset_) noexcept:
        range(range_), current_offset(current_offset_), current_index{} {
        }

      multi_dimensional_iterator& operator--() { 
      	size_t i = N-1;
      	current_index[i] --;
      	ptrdiff_t offset = range->dim_strides[i];
      	while(i && (current_index[i] < 0)){
      		offset += range->dimensions[i]*range->dim_strides[i];
      		current_index[i--] = range->dimensions[i];
      		offset -= range->dim_strides[i];  
      		current_index[i] --;
      	}
      	current_offset += offset;
      	return *this;
      }
      multi_dimensional_iterator operator--(int) { auto cpy = *this; --(*this); return cpy; }
      multi_dimensional_iterator& operator++() { 
      	size_t i = N-1;
      	current_index[i] ++;
      	ptrdiff_t offset = range->dim_strides[i];
      	while(i && (current_index[i] == range->dimensions[i])){
      		offset -= range->dimensions[i]*range->dim_strides[i];
      		current_index[i--] = 0;
      		offset += range->dim_strides[i]; 
      		current_index[i] ++;
      	}
      	current_offset += offset;
      	// std::cout << "offset=" << offset << ", current_offset=" << current_offset << std::endl;
      	return *this;
      }
      multi_dimensional_iterator operator++(int) { auto cpy = *this; ++(*this); return cpy; }

      pointer operator->() {
        return range->data[current_offset];
      }
      pointer operator->() const {
        return range->data[current_offset];
      }
      reference operator*() {
        return range->data[current_offset];
      }
      const_reference operator*() const {
        return range->data[current_offset];
      }

      bool operator==(multi_dimensional_iterator const& rhs) const {
        return current_offset == rhs.current_offset;
      }
      bool operator!=(multi_dimensional_iterator const& rhs) const {
        return current_offset != rhs.current_offset;
      }
      size_t get_current_index(size_t i) const{
      	return current_index[i];
      }
      size_t get_dimensions(size_t i) const{
      	return range->dimensions[i];
      }
      ptrdiff_t get_offset() const{
      	return current_offset;
      }
      template <class... Args>
      T prev(Args&&... pos) const{
        // TODO: check int type
      	static_assert(sizeof...(Args) == N, "Must have the same number of arguments");
      	auto offset = current_offset;
      	std::array<int, N> args{std::forward<Args>(pos)...};
      	for(int i=0; i<N; i++){
      		offset -= args[i] ? range->global_dim_strides[i] : 0;
      	}
      	return (offset >= 0) ? range->data[offset] : 0;
      }

    private:
      friend multi_dimensional_range;
      std::shared_ptr<multi_dimensional_range> range;
	  std::array<size_t, N> current_index; 		// index of current_offset position
      ptrdiff_t current_offset;
  };

  using iterator = multi_dimensional_iterator;
  using const_iterator = multi_dimensional_iterator;
  using value_type = T;
  using reference = T&;
  using pointer = T*;

  multi_dimensional_iterator begin() {
    return multi_dimensional_iterator(this->shared_from_this(), start_offset);
  }
  multi_dimensional_iterator end() {
    return multi_dimensional_iterator(this->shared_from_this(), end_offset);
  }
  /**
   * constructs a multi_dimensional_range
   *
   * \param[in] global_dims the dimensions of the overall array
   * \param[in] stride how many items to skip in the global array in each direction
   * \param[in] count how many items to include from the global array in each direction
   */
  template <class ForwardIt1>
  void set_global_dimensions(ForwardIt1 begin, ForwardIt1 end){
    int i = 0;
    for(auto iter = begin; iter != end; iter ++){
    	global_dimensions[i ++] = *iter;
    }
    size_t cur_stride = 1;
    for(int i=N-1; i>=0; i--){
    	global_dim_strides[i] = cur_stride;
    	cur_stride *= global_dimensions[i];
    }     
  }
  template <class ForwardIt1>
  void set_dimensions(ForwardIt1 begin, ForwardIt1 end){
    int i = 0;
    for(auto iter = begin; iter != end; iter ++){
    	dimensions[i ++] = *iter;
	    // std::cout << dimensions[i-1] << " ";
    }
  }
  void set_dimensions_auto(){
  	// std::cout << "dimensions: ";
    for(int i=0; i<dimensions.size(); i++){
    	// std::cout << "g[i]=" << global_dimensions[i] << ",str=" << access_stride << " ";
		dimensions[i] = (global_dimensions[i] - 1) / access_stride + 1;
	    // std::cout << dimensions[i] << " ";
    }
    // std::cout << std::endl;
  }
  void set_dim_strides(){
  	// std::cout << "strides: ";
    size_t cur_stride = 1;
    for(int i=N-1; i>=0; i--){
    	dim_strides[i] = cur_stride * access_stride;
    	cur_stride *= global_dimensions[i];
	    // std::cout << dim_strides[i] << " ";
    }  	
    // std::cout << std::endl;
  }
  void set_offsets(ptrdiff_t offset_){
  	start_offset = offset_;
  	end_offset = start_offset + dimensions[0] * dim_strides[0];
  }
  void set_access_stride(size_t stride_){
  	access_stride = stride_;
  }
  template <class ForwardIt1>
  multi_dimensional_range(
      T* data_,
      ForwardIt1 global_dims_begin,
      ForwardIt1 global_dims_end,
      size_t stride_,
      ptrdiff_t offset_
      ): data(data_)
  {
    static_assert(
        std::is_convertible<
          typename std::iterator_traits<ForwardIt1>::value_type,
          std::size_t>::value,
          "ForwardIt1 must be convertible to std::size_t"
    );
    if(global_dims_end - global_dims_begin != N){
    	std::cout << global_dims_end - global_dims_begin << " " << N << std::endl;
    	std::cerr << "#dimensions does not match!\n";
    	exit(0);
    }
    set_access_stride(stride_);
    set_global_dimensions(global_dims_begin, global_dims_end);
    // set_dimensions(dims_begin, dims_end);
    set_dimensions_auto();
    set_dim_strides();
    set_offsets(offset_);
    // std::cout << start_offset << " " << end_offset << std::endl;
  }

  size_t num_dims() const { return dimensions.size(); };

  private:
	std::array<size_t, N> global_dimensions;
	std::array<size_t, N> global_dim_strides;
  std::array<size_t, N> dimensions; 			  // the dimensions 
  std::array<size_t, N> dim_strides; 			  // strides for dimensions
  size_t access_stride;					          	// stride for access pattern
  ptrdiff_t start_offset;					          // offset for start point
  ptrdiff_t end_offset;						          // offset for end point
  T * data;									                // data pointer
};

}
#endif