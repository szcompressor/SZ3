#include <iostream>
#include <filesystem>
#include "utils/Verification.hpp"
#include "utils/Timer.hpp"
#include "ExaaltDef.hpp"
#include "InstantStore.hpp"
#include "PersistentStore.hpp"
#include "CompressedStore.hpp"

typedef float datatype;


int main(int argc, char **argv) {

    int block_size = BLOCK_SIZE;
    int buffersize = LEVEL2_BATCH_SIZE;
    size_t level1_capacity = LEVEL1_CAPACITY;
    size_t level2_capacity = LEVEL2_CAPACITY;
    double eb = 0.0001;

    if (argc <= 1) {
        std::cout << "Usage : " << argv[0] << " [datafile] [buffersize] [error_bound]"
                  << std::endl;
        std::cout << "Example: " << argv[0] << " ../../exaalt_data_5423x3137/xx.dat 1000 1E-4" << std::endl;
        std::exit(0);
    }

    std::string datafile(argv[1]);
    if (argc > 2) {
        buffersize = atoi(argv[2]);
        eb = atof(argv[3]);
    }

    size_t input_num;
    auto input = SZ::readfile<datatype>(datafile.c_str(), input_num);

    int testcase = input_num / block_size; //compute the number of blocks
    size_t num_element = testcase * block_size;
    std::vector<RawDataVector> rawData(testcase); //initialize the raw data test case vector
    std::vector<int64> keys(testcase); //initialize the keys vector

    //initialize the test cases
    srand(time(0));
    for (int i = 0; i < testcase; i++) {
        rawData[i] = std::vector<char>((char *) input.get() + i * block_size * sizeof(datatype),
                                       (char *) input.get() + (i + 1) * block_size * sizeof(datatype));
        keys[i] = rand() % testcase;
    }
    input.reset();

    std::cout << "data size: " << num_element
              << ", block size: " << block_size
              << ", number of blocks: " << testcase
              << ", buffer size: " << buffersize
              << ", error bound: " << eb
              << std::endl;

//    std::cout << "raw_char vector size: " << input_data[0].size() << std::endl;
//    std::cout << input[0] << "," << input[1] << std::endl;
//    std::cout << raw_double[0][0] << "," << raw_double[0][1] << std::endl;
    SZ::Timer timer;
    bool verify = false;
    std::vector<RawDataVector> uncompressed_vector(verify ? testcase : 1);
    std::vector<RawDataVector> decompressed_vector(verify ? testcase : 1);

    {
        std::cout << "==============buffer without compression===============" << std::endl;
        InstantStore instStore;
//    instStore.initialize(level1_capacity, block_size);
        timer.start();
        for (int i = 0; i < testcase; i++) {
            int64 t = i;
            instStore.put(1, t, rawData[i]);
        }
        timer.stop("put");

        timer.start();
        for (int i = 0; i < testcase; i++) {
            instStore.get(1, keys[i], uncompressed_vector[verify ? i : 0]);
        }
        timer.stop("get");
    }

    {
        std::cout << "==============buffer with compression===============" << std::endl;

        CompressedStore<datatype, InstantStore> cInstStore;
        cInstStore.initialize("", level2_capacity, block_size, eb, buffersize);

        timer.start();
        for (int i = 0; i < testcase; i++) {
            int64 t = i;
            cInstStore.put(1, t, rawData[i]);
        }

        int finalBufferSize = cInstStore.getBufferSize();
        if (finalBufferSize > 0) {
            cInstStore.compressBuffer();
        }
        //std::cout << "sum = " << sum << ", counter = " << counter << "\n";
        timer.stop("put");
//    std::cout << "final buffer size = " << finalBufferSize << std::endl;

        timer.start();
        for (int i = 0; i < testcase; i++) {
            cInstStore.get(1, keys[i], decompressed_vector[verify ? i : 0]);
        }
        timer.stop("get");
    }

    if (verify) {
        std::cout << "==============verification===============" << std::endl;
        std::vector<datatype> uncompressed(num_element), decompressed(num_element);
        for (int i = 0; i < testcase; i++) {
            memcpy(&uncompressed[i * block_size], uncompressed_vector[i].data(), block_size * sizeof(datatype));
            memcpy(&decompressed[i * block_size], decompressed_vector[i].data(), block_size * sizeof(datatype));
        }
        SZ::verify(uncompressed.data(), decompressed.data(), uncompressed.size());
    }
}
