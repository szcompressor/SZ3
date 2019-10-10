#ifndef _SZ_INTEGER_QUANTIZER_HPP
#define _SZ_INTEGER_QUANTIZER_HPP

namespace SZ{

// data with T1 type
// return int
template <class T1>
class PredictionBasedQuantizer{
protected:
	T1 error_bound;
	T1 error_bound_reciprocal;
	int radius; // quantization interval radius
public:
	PredictionBasedQuantizer()=default;
	PredictionBasedQuantizer(T1 eb, int r) : error_bound(eb), error_bound_reciprocal(((T1)1) / eb), radius(r){}
	~PredictionBasedQuantizer()=default;
};

template <class T1>
class LinearQuantizer : public PredictionBasedQuantizer<T1>{
public:
	LinearQuantizer()=default;
	LinearQuantizer(T1 eb, int r) : PredictionBasedQuantizer<T1>(eb, r){}
	~LinearQuantizer()=default;
	int quantize(T1 data, T1 predicted);
	int quantize(T1 data, T1 predicted, T1& dec_data);
};

template <class T1>
int LinearQuantizer<T1>::quantize(T1 data, T1 pred){
	int radius = this->radius;
	int quant_index = ((int)((pred - data) * this->error_bound_reciprocal) + 1);
	quant_index /= 2;
	return (quant_index > 0) ? (quant_index < radius ? quant_index + radius : 0) : (quant_index > -radius ? quant_index + radius : 0);
}
template <class T1>
int LinearQuantizer<T1>::quantize(T1 data, T1 pred, T1& dec_data){
	int radius = this->radius;
	int quant_index = ((int)((pred - data) * this->error_bound_reciprocal) + 1);
	dec_data = pred - quant_index * this->error_bound;
	quant_index /= 2;
	return (quant_index > 0) ? (quant_index < radius ? quant_index + radius : 0) : (quant_index > -radius ? quant_index + radius : 0);
}
template class LinearQuantizer<float>;


}
#endif
