// #include <compressor/Compressor.hpp>
#include <quantizer/IntegerQuantizer.hpp>
#include <utils/Iterator.hpp>
#include <cstdio>
void test_quantizer(){
	SZ::LinearQuantizer<float> q(1, 32);
	printf("%d\n", q.quantize(5, 0));
	float dec = 0;
	int i = q.quantize(5, 0, dec);
	printf("%d, %.4f\n", i, dec);	
}

void test_iterator(){
	std::array<int, 81> values;
	for(int i=0; i<81; i++){
		values[i] = i+1;
	}
	std::array<size_t, 4> global_dims = {3, 3, 3, 3};
	// std::array<size_t, 2> global_dims = {3, 5};
	size_t stride = 2;
	auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<int, 4>>(
	  values.data(),
	  std::begin(global_dims),
	  std::end(global_dims),
	  stride,
	  0
	  );
	auto r = inter_block_range->begin();
	std::array<size_t, 4> intra_block_dims;
	auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<int, 4>>(
	  values.data(),
	  std::begin(global_dims),
	  std::end(global_dims),
	  1,
	  0
	  );
	for(r=inter_block_range->begin(); r!=inter_block_range->end(); r++){
		  std::cout << *r << std::endl;
	  for(int i=0; i<intra_block_dims.size(); i++){
	  	size_t cur_index = r.get_current_index(i);
	  	size_t dims = r.get_dimensions(i);
	  	intra_block_dims[i] = (cur_index == dims - 1) ? global_dims[i] - cur_index * stride : stride;
	  }
	  intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
	  intra_block_range->set_offsets(r.get_offset());
	  auto r2 = intra_block_range->begin();
	  for(r2=intra_block_range->begin(); r2!=intra_block_range->end(); r2++){
	  	std::cout << *r2 << ": " << r2.prev(0) << std::endl;
	  }
	}
}
int main(){
	test_quantizer();
	test_iterator();
	return 0;
}