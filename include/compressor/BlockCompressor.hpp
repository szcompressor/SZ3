#ifndef _SZ_BLOCK_COMPRESSOR_HPP
#define _SZ_BLOCK_COMPRESSOR_HPP

#include "predictor/Predictor.hpp"
#include "quantizer/Quantizer.hpp"
#include "encoder/Encoder.hpp"
#include "def.hpp"

namespace SZ{

template <typename T>
class blockCompressor{
protected:
	const T * block_start_pos;
	Predictor<T> predictor;
	Quantizer<T> quantizer;
	Encoder<T, uint32> encoder;
public:
	blockCompressor();
	~blockCompressor();
	virtual void compress();
};

}
#endif