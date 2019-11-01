#ifndef _SZ_INTEGER_QUANTIZER_HPP
#define _SZ_INTEGER_QUANTIZER_HPP
#include <string>
#include <cassert>
#include <iostream>
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
	// int quantize(T data, T pred, T& dec_data);
	int quantize_and_overwrite(T& data, T pred);
	// recover the data using the quantization index
	T recover(T pred, int quant_index);

  void save(unsigned char*& c) const {
    // std::string serialized(sizeof(uint8_t) + sizeof(T) + sizeof(int),0);
    c[0] = 0b00000010;
    c += 1;
    std::cout << "saving eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
    *reinterpret_cast<T*>(c) = this->error_bound;
    c += sizeof(T);
    *reinterpret_cast<int*>(c) = this->radius;
    c += sizeof(int);
    *reinterpret_cast<size_t*>(c) = unpred.size();
    c += sizeof(size_t);
    memcpy(c, unpred.data(), unpred.size()*sizeof(T));
    c += unpred.size()*sizeof(T);
  };

  void load(const unsigned char*& c, size_t& remaining_length) {
    assert(remaining_length > (sizeof(uint8_t) + sizeof(T) + sizeof(int)));
    c += sizeof(uint8_t);
    remaining_length -= sizeof(uint8_t);
    this->error_bound = *reinterpret_cast<const T*>(c);
    c += sizeof(T);
    this->radius = *reinterpret_cast<const int*>(c);
    c += sizeof(int);
    size_t unpred_size = *reinterpret_cast<const size_t*>(c);
    c += sizeof(size_t);
    this->unpred = std::vector<T>(reinterpret_cast<const T*>(c), reinterpret_cast<const T*>(c) + unpred_size);
    c += unpred_size * sizeof(T);
    std::cout << "loading: eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
    // reset index
    index = 0;
    // std::cout << "unpred data size " << unpred.size() << std::endl;
    // return LinearQuantizer<T>(c);
  }

private:
	std::vector<T> unpred;
	size_t index = 0; // used in decompression only
};

template <class T>
int LinearQuantizer<T>::quantize(T data, T pred){
	int radius = this->radius;
	// compute quantization index
	int quant_index = (int)((data - pred) * this->error_bound_reciprocal);
	quant_index = (quant_index > 0) ? (quant_index + 1)/2 : (quant_index - 1)/2;
	// shift quantization index, set overbound to 0
	return (quant_index > 0) ? (quant_index < radius ? quant_index + radius : 0) : (quant_index > -radius ? quant_index + radius : 0);
}
template <class T>
int LinearQuantizer<T>::quantize_and_overwrite(T& data, T pred){
	int radius = this->radius;
	// compute quantization index
	int quant_index = (int)((data - pred) * this->error_bound_reciprocal);
	quant_index = (quant_index > 0) ? (quant_index + 1)/2 : (quant_index - 1)/2;
	// shift quantization index, set overbound to 0
	int quant_index_shifted = (quant_index > 0) ? (quant_index < radius ? quant_index + radius : 0) : (quant_index > -radius ? quant_index + radius : 0);
	if(quant_index_shifted){
		data = pred + 2 * quant_index * this->error_bound;		
	}
	else{
		unpred.push_back(data);
	}
	return quant_index_shifted;
}
template <class T>
T LinearQuantizer<T>::recover(T pred, int quant_index){
	if(quant_index){
		return pred + 2 * (quant_index - this->radius) * this->error_bound;
	}
	else{
		// if(index >= unpred.size()){
		// 	std::cout << "index = " << index << ", unpred_size = " << unpred.size() << std::endl;
		// 	exit(0);
		// }
		return unpred[index ++];
	}
}

}
#endif
