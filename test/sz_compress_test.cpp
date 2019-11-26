#include <quantizer/IntegerQuantizer.hpp>
#include <compressor/Compressor.hpp>
#include <quantizer/Quantizer.hpp>
#include <utils/Iterator.hpp>
#include <predictor/Predictor.hpp>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>

template<typename Type>
std::unique_ptr<Type[]> readfile(const char * file, size_t& num){
	std::ifstream fin(file, std::ios::binary);
	if(!fin){
        std::cout << " Error, Couldn't find the file" << "\n";
        return 0;
    }
    fin.seekg(0, std::ios::end);
    const size_t num_elements = fin.tellg() / sizeof(Type);
    fin.seekg(0, std::ios::beg);
    auto data = SZ::compat::make_unique<Type[]>(num_elements);
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
	float eb = 10;

	auto P_l = SZ::LorenzoPredictor<float, 3>(eb);
	auto P_reg = SZ::RegressionPredictor<float, 3>(0.1*eb);
	SZ::ComposedPredictor<float, 3> cp(P_l, P_reg);
	auto sz = SZ::make_sz_general<float>(
		// SZ::LorenzoPredictor<float, 3>(),
		// SZ::RegressionPredictor<float, 3>(0.1*eb),
		&cp,
		SZ::LinearQuantizer<float>(eb),
		SZ::HuffmanEncoder<int>(),
		100,
		500,
		500
	);
	//SZ::SZ_General_Compressor<float, 3> sz(lorenzo, linear_quantizer, 100, 500, 500);

	size_t compressed_size = 0;
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    std::unique_ptr<unsigned char[]> compressed;
    compressed.reset(sz.compress(data.get(), eb, compressed_size));
    err = clock_gettime(CLOCK_REALTIME, &end);
    std::cout << "Compression time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << std::endl;
	std::cout << "Compressed size = " << compressed_size << std::endl;
	writefile("test.dat", compressed.get(), compressed_size);
    err = clock_gettime(CLOCK_REALTIME, &start);
	auto dec_data = sz.decompress(compressed.get(), compressed_size);
    err = clock_gettime(CLOCK_REALTIME, &end);
    std::cout << "Decompression time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << std::endl;
	float max_err = 0;
	for(int i=0; i<num; i++){
		max_err = std::max(max_err, std::abs(data[i] - dec_data[i]));
	}
	std::cout << "Max error = " << max_err << std::endl;
	return 0;
}
