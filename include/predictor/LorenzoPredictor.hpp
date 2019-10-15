#ifndef _SZ_LORENZO_PREDICTOR_HPP
#define _SZ_LORENZO_PREDICTOR_HPP

#include "utils/Iterator.hpp"

namespace SZ{

// N-d lorenzo predictor
template <class T, uint N>
class LorenzoPredictor;

template <class T>
class LorenzoPredictor<T, 1>{
public:
	using Iterator = typename multi_dimensional_range<T, 1>::multi_dimensional_iterator;
	void preprocess(const Iterator iter) const noexcept{};
	void postprocess(const Iterator iter) const noexcept{};
	inline T predict(Iterator iter) const noexcept{
		return iter.prev(0);
	};
};

template <class T>
class LorenzoPredictor<T, 2>{
public:
	using Iterator = typename multi_dimensional_range<T, 2>::multi_dimensional_iterator;
	void preprocess(const Iterator iter) const noexcept{};
	void postprocess(const Iterator iter) const noexcept{};
	inline T predict(Iterator iter) const noexcept{
		return iter.prev(0, 1) + iter.prev(1, 0) - iter.prev(1, 1);
	};
};

template <class T>
class LorenzoPredictor<T, 3>{
public:
	using Iterator = typename multi_dimensional_range<T, 3>::multi_dimensional_iterator;
	void preprocess(const Iterator iter) const noexcept{};
	void postprocess(const Iterator iter) const noexcept{};
	inline T predict(Iterator iter) const noexcept{
		return iter.prev(0, 0, 1) + iter.prev(0, 1, 0) + iter.prev(1, 0, 0) 
				- iter.prev(0, 1, 1) - iter.prev(1, 0, 1) - iter.prev(1, 1, 0)
				+ iter.prev(1, 1, 1);
	};
};

}
#endif
