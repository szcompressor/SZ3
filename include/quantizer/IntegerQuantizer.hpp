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
	int quantize(T data, T pred);
	int quantize(T data, T pred, T& dec_data);
	T recover(T pred, int quant_index);
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
