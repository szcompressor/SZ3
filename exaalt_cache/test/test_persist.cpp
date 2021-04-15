#include <iostream>
#include <filesystem>
#include "ExaaltDef.hpp"
#include "PersistentStore.hpp"
#include "CompressedStore.hpp"
#include "utils/Timer.hpp"
#include "utils/Verification.hpp"

namespace fs = std::filesystem;
typedef float datatype;

int main(int argc, char **argv) {
    size_t block_size = BLOCK_SIZE;
    int buffersize = LEVEL2_BATCH_SIZE;
    double eb = 0.0001;

    if (argc <= 1) {
        std::cout << "Usage : " << argv[0] << " [dbdir] [datafile] [buffersize] [error_bound]" << std::endl;
        std::cout << "Example: " << argv[0] << " tmp  ../../exaalt_data_5423x3137/xx.dat 1000 1E-4" << std::endl;
        std::exit(0);
    }

    std::string dbdir(argv[1]);
    std::string datafile(argv[2]);
    if (argc > 3) {
        buffersize = atoi(argv[3]);
        eb = atof(argv[4]);
    }

//    std::filesystem::remove_all(dbdir);
//    std::filesystem::create_directories(dbdir);
//    store.initialize(dbdir, "exaalt_sz", block_size, eb, buffersize);
    size_t input_num;
//    auto input = SZ::readfile<float>(std::string(homedir + "/data/exaalt-2869440/xx.dat2").c_str(), input_num);
    auto input = SZ::readfile<float>(datafile.c_str(), input_num);

    srand(time(0));
    int testcase = input_num / block_size;
    size_t num_element = testcase * block_size;
    std::vector<RawDataVector> rawData(testcase);
    std::vector<int64> keys(testcase);
    for (int i = 0; i < testcase; i++) {
        rawData[i] = std::vector<char>((char *) input.get() + i * block_size * sizeof(datatype),
                                       (char *) input.get() + (i + 1) * block_size * sizeof(datatype));
        keys[i] = rand() % testcase;
    }
    input.reset();

    std::cout << "data size: " << input_num << ", block size: " << block_size
              << ", number of blocks: " << testcase
              << ", buffer size: " << buffersize
              << ", error bound: " << eb
              << std::endl;

    std::cout << "raw_char vector size: " << rawData[0].size() << std::endl;
//    std::cout << input[0] << "," << input[1] << std::endl;
//    std::cout << raw_double[0][0] << "," << raw_double[0][1] << std::endl;

    SZ::Timer timer;
    bool verify = false;
    std::vector<RawDataVector> uncompressed_vector(verify ? testcase : 1);
    std::vector<RawDataVector> decompressed_vector(verify ? testcase : 1);
    {
        std::cout << "==============buffer without compression===============" << std::endl;
        PersistentStore instStore;
        instStore.initialize(dbdir + "/exaalt");
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

        CompressedStore<datatype, PersistentStore> cInstStore;
        cInstStore.initialize(dbdir, LEVEL2_CAPACITY, block_size, eb, buffersize);

        timer.start();
        for (int i = 0; i < testcase; i++) {
            int64 t = i;
            cInstStore.put(1, t, rawData[i]);
        }

        int finalBufferSize = cInstStore.getBufferSize();
        if (finalBufferSize > 0) {
            cInstStore.compressBuffer();
        }
        timer.stop("put");

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
