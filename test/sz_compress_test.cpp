#include "quantizer/IntegerQuantizer.hpp"
#include <compressor/Compressor.hpp>
#include <quantizer/Quantizer.hpp>
#include <utils/Iterator.hpp>
#include <predictor/Predictor.hpp>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <fstream>
#include <cmath>

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
	size_t num = 0;
	// use Hurricane for testing
	auto data = readfile<float>(argv[1], num);
	std::cout << "Read " << num << " elements\n";
  auto sz = SZ::make_sz_general<float>(
    SZ::LorenzoPredictor<float, 3>(),
    SZ::LinearQuantizer<float>(0.0001f),
    100,
    500,
    500
  );
	//SZ::SZ_General_Compressor<float, 3> sz(lorenzo, linear_quantizer, 100, 500, 500);

	size_t compressed_size = 0;
	auto compressed = sz.compress(data, 0.0001, compressed_size);

	std::cout << compressed_size << std::endl;
	auto dec_data = sz.decompress(compressed, compressed_size);
	float max_err = 0;
	for(int i=0; i<num; i++){
		max_err = std::max(max_err, std::abs(data[i] - dec_data[i]));
	}
	std::cout << max_err << std::endl;
	return 0;
}
