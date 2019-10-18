#ifndef _SZ_SZ_GENERAL_HPP
#define _SZ_SZ_GENERAL_HPP

#include "predictor/Predictor.hpp"
#include "quantizer/Quantizer.hpp"
#include "utils/Iterator.hpp"
#include "def.hpp"
#include <cstring>

namespace SZ{
// type T
template <class T, size_t N, uint block_size_ = 0, class Predictor = LorenzoPredictor<T,N>, 
	class Quantizer = LinearQuantizer<T> >
class SZ_General_Compressor{
public:
  static_assert(concepts::is_predictor<Predictor>::value, "must implement the predictor interface");
  static_assert(concepts::is_quantizer<Quantizer>::value, "must implement the quatizer interface");

	template <class ... Args>
	SZ_General_Compressor(Predictor predictor, Quantizer quantizer, Args&&... dims) : predictor(predictor), quantizer(quantizer), block_size(block_size_), global_dimensions{static_cast<size_t>(dims)...}{
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
		std::vector<T> data = std::vector<T>(data_, data_+num_elements);

		auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data.data(), 
		  std::begin(global_dimensions), std::end(global_dimensions), block_size, 0);
		auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data.data(),
		  std::begin(global_dimensions), std::end(global_dimensions), 1, 0);
		int quant_radius = 32768;
		std::array<size_t, N> intra_block_dims;
		std::vector<int> quant_inds(num_elements);
		std::vector<T> unpred_data;
		int count = 0;
    	predictor.precompress_data(inter_block_range->begin());
    	quantizer.precompress_data();
		for(auto block=inter_block_range->begin(); block!=inter_block_range->end(); block++){
		  // std::cout << *block << " " << lp.predict(block) << std::endl;
		  for(int i=0; i<intra_block_dims.size(); i++){
		  	size_t cur_index = block.get_current_index(i);
		  	size_t dims = inter_block_range->get_dimensions(i);
		  	intra_block_dims[i] = (cur_index == dims - 1) ? global_dimensions[i] - cur_index * block_size : block_size;
		  }
		  intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
		  intra_block_range->set_offsets(block.get_offset());
	  	  T dec_data = 0;
      	predictor.precompress_block(intra_block_range->begin());
      	quantizer.precompress_block();
		  for(auto element=intra_block_range->begin(); element!=intra_block_range->end(); element++){
		  	quant_inds[count ++] = quantizer.quantize_and_overwrite(*element, predictor.predict(element));
		  }
		}
    	predictor.postcompress_data(inter_block_range->begin());
    	quantizer.postcompress_data();

		std::cout << "unpred_num = " << unpred_data.size() << std::endl;
		uchar* compressed_data = new uchar[2 * num_elements * sizeof(T)];
		uchar* compressed_data_pos = compressed_data;
		// TODO: serialize and record predictor, quantizer, and encoder
		// Or do these in a outer loop wrapper?
		write(global_dimensions.data(), N, compressed_data_pos);
		write(block_size, compressed_data_pos);

    	auto serialized_predictor = predictor.save();
    	write(serialized_predictor.data(), serialized_predictor.size(), compressed_data_pos);
	    quantizer.save(compressed_data_pos);
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

	T * decompress(uchar const * compressed_data, const size_t length){
		uchar const * compressed_data_pos = compressed_data;
    size_t remaining_length = length;
		read(global_dimensions.data(), N, compressed_data_pos, remaining_length);
		num_elements = 1;
		for(const auto& d : global_dimensions){
			num_elements *= d;
		}
		uint block_size = 0;
		read(block_size, compressed_data_pos, remaining_length);

    	predictor.load(compressed_data_pos, remaining_length);
    	quantizer.load(compressed_data_pos, remaining_length);
		size_t unpred_data_size = 0;
		read(unpred_data_size, compressed_data_pos, remaining_length);
		T const * unpred_data_pos = (T const *) compressed_data_pos;
		compressed_data_pos += unpred_data_size * sizeof(T);
		int const * quant_inds_pos = (int const *) compressed_data_pos;
		std::array<size_t, N> intra_block_dims;
		std::vector<T> dec_data = std::vector<T>(num_elements);
		auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(dec_data.data(), 
		  std::begin(global_dimensions), std::end(global_dimensions), block_size, 0);
    	predictor.predecompress_data(inter_block_range->begin());
    	quantizer.predecompress_data();

		auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(dec_data.data(),
		  std::begin(global_dimensions), std::end(global_dimensions), 1, 0);
		for(auto block=inter_block_range->begin(); block!=inter_block_range->end(); block++){
		  // std::cout << *block << " " << lp.predict(block) << std::endl;
		  for(int i=0; i<intra_block_dims.size(); i++){
		  	size_t cur_index = block.get_current_index(i);
		  	size_t dims = inter_block_range->get_dimensions(i);
		  	intra_block_dims[i] = (cur_index == dims - 1) ? global_dimensions[i] - cur_index * block_size : block_size;
		  }
		  intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
		  intra_block_range->set_offsets(block.get_offset());

      	  predictor.predecompress_block(intra_block_range->begin());
      	  quantizer.predecompress_block();
		  for(auto element=intra_block_range->begin(); element!=intra_block_range->end(); element++){
	  		*element = quantizer.recover(predictor.predict(element), *(quant_inds_pos ++));
		  }
		}
    	predictor.postdecompress_data(inter_block_range->begin());
    	quantizer.postdecompress_data();
		return dec_data.data();
	}
	// read array
	template <class T1>
	void read(T1 * array, size_t num_elements, uchar const *& compressed_data_pos, size_t& remaining_length){
    assert(num_elements*sizeof(T1) < remaining_length);
		memcpy(array, compressed_data_pos, num_elements*sizeof(T1));
    	remaining_length -= num_elements*sizeof(T1);
		compressed_data_pos += num_elements*sizeof(T1);
	}
	// read variable
	template <class T1>
	void read(T1& var, uchar const *& compressed_data_pos, size_t& remaining_length ){
    assert(sizeof(T1) < remaining_length);
		memcpy(&var, compressed_data_pos, sizeof(T1));
    remaining_length -= sizeof(T1);
		compressed_data_pos += sizeof(T1);
	}

private:
  Predictor predictor;
  Quantizer quantizer;
	uint block_size;
	size_t num_elements;
	std::array<size_t, N> global_dimensions;
};


template <class T, uint block_size_ = 0, class Predictor, class Quantizer, class... Args>
SZ_General_Compressor<T,sizeof...(Args),block_size_,Predictor, Quantizer>
make_sz_general(
    Predictor predictor,
    Quantizer quantizer,
    Args... args
    ) {
  return SZ_General_Compressor<T, sizeof...(Args), block_size_, Predictor, Quantizer>(predictor, quantizer, args...);
}


}
#endif
