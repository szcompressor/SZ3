#include <compressor/Compressor.hpp>
#include <quantizer/Quantizer.hpp>
#include <utils/Iterator.hpp>
#include <predictor/Predictor.hpp>
#include <cstdio>
#include <fstream>
#include <iostream>

void test_quantizer(){
	SZ::LinearQuantizer<float> q(1, 32);
	printf("%d\n", q.quantize(5, 0));
	float dec = 0;
	int i = q.quantize(5, 0, dec);
	printf("%d, %.4f\n", i, dec);	
}

void test_iterator(){
	std::array<int, 45> values;
	for(int i=0; i<45; i++){
		values[i] = i+1;
	}
	std::array<size_t, 3> global_dims = {3, 3, 5};
	// std::array<size_t, 2> global_dims = {3, 5};
	size_t stride = 2;
	auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<int, 3>>(
	  values.data(),
	  std::begin(global_dims),
	  std::end(global_dims),
	  stride,
	  0
	  );
	std::array<size_t, 3> intra_block_dims;
	auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<int, 3>>(
	  values.data(),
	  std::begin(global_dims),
	  std::end(global_dims),
	  1,
	  0
	  );
	SZ::LorenzoPredictor<int, 3> lp;
	for(auto r=inter_block_range->begin(); r!=inter_block_range->end(); r++){
	  std::cout << *r << " " << lp.predict(r) << std::endl;
	  for(int i=0; i<intra_block_dims.size(); i++){
	  	size_t cur_index = r.get_current_index(i);
	  	size_t dims = r.get_dimensions(i);
	  	intra_block_dims[i] = (cur_index == dims - 1) ? global_dims[i] - cur_index * stride : stride;
	  }
	  intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
	  intra_block_range->set_offsets(r.get_offset());
	  auto r2 = intra_block_range->begin();
	}
}

template<typename Type>
Type * readfile(const char * file, size_t& num){
	std::ifstream fin(file, std::ios::binary);
	if(!fin){
        std::cout << " Error, Couldn't find the file" << "\n";
        return 0;
    }
    fin.seekg(0, std::ios::end);
    const size_t num_elements = fin.tellg() / sizeof(Type);
    fin.seekg(0, std::ios::beg);
    Type * data = (Type *) malloc(num_elements*sizeof(Type));
	fin.read(reinterpret_cast<char*>(&data[0]), num_elements*sizeof(Type));
	fin.close();
	num = num_elements;
	return data;
}

template<typename Type>
void writefile(const char * file, Type * data, size_t num_elements){
	std::ofstream fout(file, std::ios::binary);
	fout.write(reinterpret_cast<const char*>(&data[0]), num_elements*sizeof(Type));
	fout.close();
}

int main(int argc, char ** argv){
	// test_quantizer();
	// test_iterator();
	size_t num = 0;
	auto data = readfile<float>(argv[1], num);
	std::cout << "Read " << num << " elements\n";
	SZ::SZ_General_Compressor<float, 3> sz(data, 100, 500, 500);
	sz.compress(0.1);
	return 0;
}