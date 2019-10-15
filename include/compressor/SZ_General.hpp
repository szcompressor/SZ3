#ifndef _SZ_SZ_GENERAL_HPP
#define _SZ_SZ_GENERAL_HPP

#include "predictor/Predictor.hpp"
#include "quantizer/Quantizer.hpp"
#include "utils/Iterator.hpp"
#include "def.hpp"

namespace SZ{
// type T
template <class T, size_t N, uint block_size_ = 0, template<class, uint> class Predictor = LorenzoPredictor, 
	template <class> class Quantizer = LinearQuantizer >
class SZ_General_Compressor{
public:
	SZ_General_Compressor()=default;
	~SZ_General_Compressor()=default;
	template <class ... Args>
	SZ_General_Compressor(Args&&... dims) : block_size(block_size_), global_dimensions{static_cast<size_t>(dims)...}{
      	static_assert(sizeof...(Args) == N, "Number of arguments must be the same as N");
		num_elements = 1;
		for(const auto& d:global_dimensions){
			num_elements *= d;
		}
		if(block_size == 0){
			switch(N){
				case 1:
					block_size = 128;
					break;
				case 2:
					block_size = 16;
					break;
				default:
					// >= 3D
					block_size = 6;
					break;

			}
		}
	}
	// compress given the error bound
	uchar * compress (const T * data_, double eb, size_t& compressed_size){
		// make a copy of the data		
		std::vector<T> data = std::vector<T>(num_elements);
		memcpy(data.data(), data_, num_elements * sizeof(T));

		auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data.data(), 
		  std::begin(global_dimensions), std::end(global_dimensions), block_size, 0);
		auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data.data(),
		  std::begin(global_dimensions), std::end(global_dimensions), 1, 0);
		Predictor<T, N> predictor;
		int quant_radius = 32768;
		Quantizer<T> quantizer(eb, quant_radius);
		std::array<size_t, N> intra_block_dims;
		std::vector<int> quant_inds(num_elements);
		std::vector<T> unpred_data;
		int count = 0;
		for(auto block=inter_block_range->begin(); block!=inter_block_range->end(); block++){
		  // std::cout << *block << " " << lp.predict(block) << std::endl;
		  for(int i=0; i<intra_block_dims.size(); i++){
		  	size_t cur_index = block.get_current_index(i);
		  	size_t dims = block.get_dimensions(i);
		  	intra_block_dims[i] = (cur_index == dims - 1) ? global_dimensions[i] - cur_index * block_size : block_size;
		  }
		  intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
		  intra_block_range->set_offsets(block.get_offset());
	  	  T dec_data = 0;
		  for(auto element=intra_block_range->begin(); element!=intra_block_range->end(); element++){
		  	quant_inds[count] = quantizer.quantize(*element, predictor.predict(element), dec_data);
		  	if(quant_inds[count]){
		  		*element = dec_data;
		  	}
		  	else{
		  		unpred_data.push_back(*element);
		  	}
		  	count ++;
		  }
		}
		std::cout << "unpred_num = " << unpred_data.size() << std::endl;
		uchar * compressed_data = (uchar *) malloc(2 * num_elements * sizeof(T));
		uchar * compressed_data_pos = compressed_data;
		// TODO: serialize and record predictor, quantizer, and encoder
		// Or do these in a outer loop wrapper?
		write(global_dimensions.data(), N, compressed_data_pos);
		write(block_size, compressed_data_pos);
		write(eb, compressed_data_pos);
		write(quant_radius, compressed_data_pos);
		write(unpred_data.size(), compressed_data_pos);
		write(unpred_data.data(), unpred_data.size(), compressed_data_pos);
		write(quant_inds.data(), quant_inds.size(), compressed_data_pos);
		compressed_size = compressed_data_pos - compressed_data;
		return compressed_data;
	}
	// write array
	template <class T1>
	void write(T1 const * array, size_t num_elements, uchar *& compressed_data_pos){
		memcpy(compressed_data_pos, array, num_elements*sizeof(T1));
		compressed_data_pos += num_elements*sizeof(T1);
	}
	// write variable
	template <class T1>
	void write(T1 const var, uchar *& compressed_data_pos){
		memcpy(compressed_data_pos, &var, sizeof(T1));
		compressed_data_pos += sizeof(T1);
	}

private:
	uint block_size;
	size_t num_elements;
	std::array<size_t, N> global_dimensions;
};

template <class T, size_t N, template<class, uint> class Predictor = LorenzoPredictor, 
	template <class> class Quantizer = LinearQuantizer >
class SZ_General_Decompressor{
public:
	SZ_General_Decompressor()=default;
	~SZ_General_Decompressor()=default;
	T * decompress(uchar const * compressed_data){
		uchar const * compressed_data_pos = compressed_data;
		read(global_dimensions.data(), N, compressed_data_pos);
		num_elements = 1;
		for(const auto& d : global_dimensions){
			num_elements *= d;
		}
		uint block_size = 0;
		read(block_size, compressed_data_pos);
		double eb = 0;
		read(eb, compressed_data_pos);
		int quant_radius = 0;
		read(quant_radius, compressed_data_pos);
		size_t unpred_data_size = 0;
		read(unpred_data_size, compressed_data_pos);
		T const * unpred_data_pos = (T const *) compressed_data_pos;
		compressed_data_pos += unpred_data_size * sizeof(T);
		int const * quant_inds_pos = (int const *) compressed_data_pos;
		std::array<size_t, N> intra_block_dims;
		std::vector<T> dec_data = std::vector<T>(num_elements);
		auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(dec_data.data(), 
		  std::begin(global_dimensions), std::end(global_dimensions), block_size, 0);
		auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(dec_data.data(),
		  std::begin(global_dimensions), std::end(global_dimensions), 1, 0);
		Predictor<T, N> predictor;
		Quantizer<T> quantizer(eb, quant_radius);
		for(auto block=inter_block_range->begin(); block!=inter_block_range->end(); block++){
		  // std::cout << *block << " " << lp.predict(block) << std::endl;
		  for(int i=0; i<intra_block_dims.size(); i++){
		  	size_t cur_index = block.get_current_index(i);
		  	size_t dims = block.get_dimensions(i);
		  	intra_block_dims[i] = (cur_index == dims - 1) ? global_dimensions[i] - cur_index * block_size : block_size;
		  }
		  intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
		  intra_block_range->set_offsets(block.get_offset());
		  // preprocess
		  for(auto element=intra_block_range->begin(); element!=intra_block_range->end(); element++){
		  	int quant_index = *(quant_inds_pos ++);
		  	if(quant_index){
		  		*element = quantizer.recover(predictor.predict(element), quant_index);
		  	}
		  	else{
		  		*element = *(unpred_data_pos ++);
		  	}
		  }
		}
		return dec_data.data();
	}
	// read array
	template <class T1>
	void read(T1 * array, size_t num_elements, uchar const *& compressed_data_pos){
		memcpy(array, compressed_data_pos, num_elements*sizeof(T1));
		compressed_data_pos += num_elements*sizeof(T1);
	}
	// read variable
	template <class T1>
	void read(T1& var, uchar const *& compressed_data_pos ){
		memcpy(&var, compressed_data_pos, sizeof(T1));
		compressed_data_pos += sizeof(T1);
	}

private:
	uint block_size;
	size_t num_elements;
	std::array<size_t, N> global_dimensions;
};


}
#endif
