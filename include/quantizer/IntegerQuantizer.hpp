#ifndef _SZ_INTEGER_QUANTIZER_HPP
#define _SZ_INTEGER_QUANTIZER_HPP

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
	PredictionBasedQuantizer()=default;
	PredictionBasedQuantizer(T eb, int r) : error_bound(eb), error_bound_reciprocal(((T)1) / eb), radius(r){}
	~PredictionBasedQuantizer()=default;
};

template <class T>
class LinearQuantizer : public PredictionBasedQuantizer<T>{
public:
	LinearQuantizer()=default;
	LinearQuantizer(T eb, int r) : PredictionBasedQuantizer<T>(eb, r){}
	~LinearQuantizer()=default;
	int quantize(T data, T predicted);
	int quantize(T data, T predicted, T& dec_data);
};

template <class T>
int LinearQuantizer<T>::quantize(T data, T pred){
	int radius = this->radius;
	int quant_index = ((int)((pred - data) * this->error_bound_reciprocal) + 1);
	quant_index /= 2;
	return (quant_index > 0) ? (quant_index < radius ? quant_index + radius : 0) : (quant_index > -radius ? quant_index + radius : 0);
}
template <class T>
int LinearQuantizer<T>::quantize(T data, T pred, T& dec_data){
	int radius = this->radius;
	int quant_index = ((int)((pred - data) * this->error_bound_reciprocal) + 1);
	dec_data = pred - quant_index * this->error_bound;
	quant_index /= 2;
	return (quant_index > 0) ? (quant_index < radius ? quant_index + radius : 0) : (quant_index > -radius ? quant_index + radius : 0);
}


}
#endif
