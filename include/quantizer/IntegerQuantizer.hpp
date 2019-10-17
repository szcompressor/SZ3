#ifndef _SZ_INTEGER_QUANTIZER_HPP
#define _SZ_INTEGER_QUANTIZER_HPP
#include <string>
#include <cassert>

namespace SZ{

// data with T type
// return int
template <class T>
class PredictionBasedQuantizer{
protected:
	T error_bound;
	T error_bound_reciprocal;
	int radius; // quantization interval radius
public:
	~PredictionBasedQuantizer()=default;
	PredictionBasedQuantizer()=default;
	PredictionBasedQuantizer(PredictionBasedQuantizer const&)=default;
	PredictionBasedQuantizer(PredictionBasedQuantizer &&)=default;
	PredictionBasedQuantizer& operator=(PredictionBasedQuantizer const&)=default;
	PredictionBasedQuantizer& operator=(PredictionBasedQuantizer &&)=default;

  void precompress_data() {}
  void postcompress_data() {}
  void precompress_block() {}
  void predecompress_data() {}
  void postdecompress_data() {}
  void predecompress_block() {}


	PredictionBasedQuantizer(T eb, int r) : error_bound(eb), error_bound_reciprocal(1.0 / eb), radius(r){}
};

template <class T>
class LinearQuantizer : public PredictionBasedQuantizer<T>{
public:
	LinearQuantizer(T eb, int r = 32768) : PredictionBasedQuantizer<T>(eb, r){}

  using value_type = T;
  using reference = T&;

	// quantize the data with a prediction value, and returns the quantization index
	int quantize(T data, T pred);
	// quantize the data with a prediction value, and returns the quantization index and the decompressed data
	int quantize(T data, T pred, T& dec_data);
	// recover the data using the quantization index
	T recover(T pred, int quant_index);

  std::string save() const {
    std::string serialized(1 + sizeof(T) + sizeof(int),0);
    serialized.front() = 0b00000010;
    *reinterpret_cast<T*>(serialized[1]) = this->error_bound;
    *reinterpret_cast<T*>(serialized[1 + sizeof(T)]) = this->radius;
    return serialized;
  };

  static LinearQuantizer<T> load(const unsigned char*& c, size_t& remaining_length) {
    assert(remaining_length > (sizeof(uint8_t) + sizeof(int) + sizeof(T)));
    c += 1;
    remaining_length -= sizeof(uint8_t);
    T error_bound;
    int radius;
    return LinearQuantizer<T>{error_bound, radius};
  }

};

template <class T>
int LinearQuantizer<T>::quantize(T data, T pred){
	int radius = this->radius;
	int quant_index = (int)((data - pred) * this->error_bound_reciprocal);
	quant_index = (quant_index > 0) ? (quant_index + 1)/2 : (quant_index - 1)/2;
	return (quant_index > 0) ? (quant_index < radius ? quant_index + radius : 0) : (quant_index > -radius ? quant_index + radius : 0);
}
template <class T>
int LinearQuantizer<T>::quantize(T data, T pred, T& dec_data){
	int radius = this->radius;
	int quant_index = (int)((data - pred) * this->error_bound_reciprocal);
	quant_index = (quant_index > 0) ? (quant_index + 1)/2 : (quant_index - 1)/2;
	dec_data = pred + 2 * quant_index * this->error_bound;
	return (quant_index > 0) ? (quant_index < radius ? quant_index + radius : 0) : (quant_index > -radius ? quant_index + radius : 0);
}
template <class T>
T LinearQuantizer<T>::recover(T pred, int quant_index){
	return pred + 2 * (quant_index - this->radius) * this->error_bound;
}

}
#endif
