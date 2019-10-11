#ifndef _SZ_PREDICTOR_HPP
#define _SZ_PREDICTOR_HPP

#include "utils/Iterator.hpp"

namespace SZ{

// N-d lorenzo predictor
template <class T, uint N>
class LorenzoPredictor{
public:
	using Iterator = typename multi_dimensional_range<T, N>::multi_dimensional_iterator;
	LorenzoPredictor()=default;
	~LorenzoPredictor()=default;
	T predict(Iterator iter){return 0;};
};

template <class T>
class LorenzoPredictor<T, 1>{
public:
	using Iterator_1 = typename multi_dimensional_range<T, 1>::multi_dimensional_iterator;
	LorenzoPredictor()=default;
	~LorenzoPredictor()=default;
	T predict(Iterator_1 iter){
		return iter.prev(0);
	};
};

template <class T>
class LorenzoPredictor<T, 2>{
public:
	using Iterator_2 = typename multi_dimensional_range<T, 2>::multi_dimensional_iterator;
	LorenzoPredictor()=default;
	~LorenzoPredictor()=default;
	T predict(Iterator_2 iter){
		return iter.prev(0, 1) + iter.prev(1, 0) - iter.prev(1, 1);
	};
};

template <class T>
class LorenzoPredictor<T, 3>{
public:
	using Iterator_3 = typename multi_dimensional_range<T, 3>::multi_dimensional_iterator;
	LorenzoPredictor()=default;
	~LorenzoPredictor()=default;
	T predict(Iterator_3 iter){
		// std::cout << "3D predict\n";
		return iter.prev(0, 0, 1) + iter.prev(0, 1, 0) + iter.prev(1, 0, 0) 
				- iter.prev(0, 1, 1) - iter.prev(1, 0, 1) - iter.prev(1, 1, 0)
				+ iter.prev(1, 1, 1);
	};
};

}
#endif
