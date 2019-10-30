#ifndef _SZ_LORENZO_PREDICTOR_HPP
#define _SZ_LORENZO_PREDICTOR_HPP

#include "utils/Iterator.hpp"
#include "quantizer/Quantizer.hpp"

namespace SZ{

// N-d lorenzo predictor
template <class T, uint N>
class RegressionPredictor{
public:
  using Range = typename multi_dimensional_range<T, N>;
  using iterator = typename multi_dimensional_range<T, N>::iterator;
  void precompress_block(const Range& range) const noexcept {
  	std::array<size_t, N> dims;
  	size_t num_elements = 1;
  	for(int i=0; i<N; i++){
  		dims[i] = range->get_dimensions(i);
  		num_elements *= dims[i];
  	}
  	std::array<T, N + 1> coeffs = compute_regression_coefficients(range, dims);
    pred_and_quantize_coefficients(coeffs);
    std::copy(coeffs, coeffs + N + 1, current_coeffs);
  }
  void preprocess(const iterator iter) const noexcept {
  }
  void postprocess(const iterator iter) const noexcept {}
  inline T predict(iterator iter) const noexcept {
  	T pred = 0;
  	for(int i=0; i<N; i++){
  		pred += iter->get_current_index(i) * current_coeffs[i];
  	}
  	pred += current_coeffs[N];
  	return pred;
  }
private:
  LinearQuantizer<T> quantizer;
  std::vector<int> regression_coeff_quant_inds;
  size_t regression_coeff_index;
  std::array<T, N + 1> current_coeffs{0};
  std::array<T, N + 1> compute_regression_coefficients(const Range& range, const std::array<size_t, N>& dims){
  	std::array<double, N + 1> sum{0};
  	for(const auto& iter : range){
  		T data = *iter;
  		for(int i=0; i<N; i++){
  			sum[i] = iter->get_current_index(i) * data;
  		}
  		sum[N] += data;
  	}
  	T num_elements_recip = 1.0 / num_elements;
  	std::array<T, N> coeffs;
  	coeffs[N] = sum[N] * num_elements_recip;
  	for(int i=0; i<N; i++){
  		coeffs[i] = (2 * sum[i] / (dims[i] - 1) - sum[N]) * 6 * num_elements_recip / (dims[i] + 1);
  		coeffs[N] -= (dims[i] - 1) * coeffs[i] / 2;
  	}
  	return coeffs;
  }
  void pred_and_quantize_coefficients(std::array<T, N + 1>& coeffs){
  	for(int i=0; i<N; i++){
  		regression_coeff_quant_inds[regression_coeff_index] = quantizer.quantize_and_overwrite(coeffs[i], current_coeffs[i]);
  	}
  }
};

}

#endif