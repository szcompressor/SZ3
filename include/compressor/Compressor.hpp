#ifndef _SZ_COMPRESSOR_HPP
#define _SZ_COMPRESSOR_HPP

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
	~SZ_General_Compressor(){
		if(data) free(data);
		if(compressed_data) free(compressed_data);
	}
	template <class ... Args>
	SZ_General_Compressor(T const * data_, Args&&... dims) : block_size(block_size_), global_dimensions{static_cast<size_t>(dims)...}{
		// const uint num_args = sizeof...(Args);
		// std::cout << num_args << std::endl;
		// if(num_args > 3){
		// 	// reduce to 3D
		// 	std::array<size_t, num_args> args{static_cast<size_t>(dims)...};
		// 	global_dimensions = std::vector<size_t>(3, 1);
		// 	int i = num_args - 1;
		// 	global_dimensions[2] = args[i --];
		// 	global_dimensions[1] = args[i --];
		// 	while(i >= 0){
		// 		global_dimensions[0] *= args[i --];
		// 	}
		// }
		// else{
		// 	global_dimensions = std::vector<size_t>{static_cast<size_t>(dims)...};
		// }
		// N = num_args;
      	static_assert(sizeof...(Args) == N, "Number of arguments must be the same as N");
		num_elements = 1;
		for(const auto& d:global_dimensions){
			num_elements *= d;
		}
		data = (T *) malloc(num_elements * sizeof(T));
		memcpy(data, data_, num_elements * sizeof(T));
		compressed_data = nullptr;
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
	uint get_num_of_global_dimensions() const{ return global_dimensions.size();}
	size_t get_global_dimensions(uint i) const{ 
		return global_dimensions[i];
	}
	// compress given the error bound
	void compress(double eb){
		auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data, 
		  std::begin(global_dimensions), std::end(global_dimensions), block_size, 0);
		auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data,
		  std::begin(global_dimensions), std::end(global_dimensions), 1, 0);
		Predictor<T, N> predictor;
		Quantizer<T> quantizer(eb, 65536);
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
		// compressed_data = (uchar *) malloc(num_elements * sizeof(T));
	}

private:
	T * data;
	uchar * compressed_data;
	uint block_size;
	size_t num_elements;
	std::array<size_t, N> global_dimensions;
};

}
#endif
