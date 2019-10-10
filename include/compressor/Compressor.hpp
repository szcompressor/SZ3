#ifndef _SZ_COMPRESSOR_HPP
#define _SZ_COMPRESSOR_HPP

#include "compressor/BlockCompressor.hpp"
#include "def.hpp"

namespace SZ{

template <typename T>
class Compressor{
protected:
	const T * data;
	uint32 * compressed_data;
	uint32 * compressed_data_pos;
	std::vector<int> dimensions;
	size_t num_elements;
	// BlockIterator<T> inter_block_iter;
	// BlockCompressor<T> block_compressor;
public:
	Compressor();
	Compressor(int r1);
	Compressor(int r1, int r2);
	Compressor(int r1, int r2, int r3);
	~Compressor();
	virtual void compress(const T * data) = 0;
	virtual void decompress(const uint32 * ) = 0;
	void setDimensions(std::vector<int> dim);
	std::vector<int> getDimensions() const;
	template <typename T1>
	void write(T1 variable);
	template <typename T1>
	void write(T1 * array, size_t size);
};

}
#endif
