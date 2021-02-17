#include <iostream>
#include <filesystem>
#include "ExaaltDef.hpp"
#include "MultiLevelStoreL1.hpp"
#include "MultiLevelStoreL2.hpp"
#include "utils/Timer.hpp"
#include "utils/Verification.hpp"
#include <algorithm>
#include <random>
#include <fstream>

namespace fs = std::filesystem;
typedef float datatype;

int main(int argc, char **argv) {

    if (argc <= 1) {
        std::cout << "Usage : " << argv[0]
                  << " dbdir datafile [error_bound block_size l1_capacity l2_capacity l2_batch_size]"
                  << std::endl;
        std::cout << "Example: " << argv[0] << " tmp  ../../exaalt_data_5423x3137/xx.dat 1E-4 3137 10000 20000 400"
                  << std::endl;
        std::exit(0);
    }

    std::string dbdir(argv[1]);
    std::string datafile(argv[2]);

    int argp = 3;
    double error_bound = ERROR_BOUND;
    if (argp < argc) {
        error_bound = atof(argv[argp++]);
    }
    size_t block_size = BLOCK_SIZE;
    if (argp < argc) {
        block_size = atof(argv[argp++]);
    }
    int64 l1_capacity = LEVEL1_CAPACITY;
    if (argp < argc) {
        l1_capacity = atof(argv[argp++]);
    }
    int64 l2_capacity = LEVEL2_CAPACITY;
    if (argp < argc) {
        l2_capacity = atof(argv[argp++]);
    }
    int64 l2_batch_size = LEVEL2_BATCH_SIZE;
    if (argp < argc) {
        l2_batch_size = atof(argv[argp++]);
    }

    size_t input_num;
    auto input = SZ::readfile<float>(datafile.c_str(), input_num);

    int testcase = input_num / block_size;
    size_t num_element = testcase * block_size;
    std::cout << "data size = " << input_num
              << ", testing data size = " << num_element
              << ", error bound = " << error_bound
              << std::endl
              << "block size = " << block_size
              << ", number of blocks = " << testcase
              << ", compress batch size = " << l2_batch_size
              << std::endl
              << "instant store capacity = " << l1_capacity
              << ", compressed store capacity = " << l2_capacity
              << std::endl;

    std::vector<RawDataVector> data(testcase);
    std::vector<int64> index(testcase);
    for (int i = 0; i < testcase; i++) {
        data[i] = std::vector<char>((char *) input.get() + i * block_size * sizeof(datatype),
                                    (char *) input.get() + (i + 1) * block_size * sizeof(datatype));
        index[i] = rand() % testcase;
    }
    input.reset();
//    std::random_device rd;
//    std::default_random_engine rng(rd());
//    std::shuffle(std::begin(index), std::end(index), rng);


//    std::cout << input[0] << "," << input[1] << std::endl;
//    std::cout << raw_double[0][0] << "," << raw_double[0][1] << std::endl;

    SZ::Timer timer;
    bool verify = false;
    std::vector<RawDataVector> uncompressed_vector(verify ? testcase : 1);
    std::vector<RawDataVector> decompressed_vector(verify ? testcase : 1);
    std::string homedir = getenv("HOME");

    {
        std::cout << "==============buffer without compression===============" << std::endl;
        MultiLevelStoreL1<datatype> persistStore;
        persistStore.initialize(dbdir, l1_capacity + l2_capacity / 10, 0, block_size, error_bound,
                                l2_batch_size);

        timer.start();
        SZ::Timer timer2;
        double gettime = 0, puttime = 0;

        int64 idx;
        std::string op;
        std::ifstream infile((homedir + "/code/sz-kai/exaalt/memory_trace.txt").c_str());
        while (infile >> op >> idx) {
            if (idx >= testcase) {
                continue;
            }
            if (op == "get") {
                timer2.start();
                persistStore.get(1, idx, data[idx]);
                gettime += timer2.stop();
            } else {
                timer2.start();
                persistStore.put(1, idx, data[idx]);
                puttime = timer2.stop();
            }
        }

        timer.stop("no compress");
        persistStore.print();
        printf("get = %.2f , put = %.2f\n", gettime, puttime);
    }
    {

        std::cout << "==============buffer with compression===============" << std::endl;

        MultiLevelStoreL2<datatype> cPersistStore;
        cPersistStore.initialize(dbdir, l1_capacity, l2_capacity - l2_batch_size, block_size, error_bound,
                                 l2_batch_size);

        SZ::Timer timer2;
        double gettime = 0, puttime = 0;

        timer.start();
        int64 idx;
        std::string op;
        std::ifstream infile((homedir + "/code/sz-kai/exaalt/memory_trace.txt").c_str());
        while (infile >> op >> idx) {
            if (idx >= testcase) {
                continue;
            }
            if (op == "get") {
                timer2.start();
                cPersistStore.get(1, idx, data[idx]);
                gettime += timer2.stop();
            } else {
                timer2.start();
                cPersistStore.put(1, idx, data[idx]);
                puttime = timer2.stop();
            }
        }

        timer.stop("compress");
        cPersistStore.print();
        printf("get = %.2f , put = %.2f\n", gettime, puttime);

    }
    if (verify) {
        std::cout << "==============verification===============" << std::endl;

        std::vector<datatype> uncompressed(num_element), decompressed(num_element);
        for (int i = 0; i < testcase; i++) {
            memcpy(&uncompressed[i * block_size], uncompressed_vector[i].data(), block_size * sizeof(datatype));
            memcpy(&decompressed[i * block_size], decompressed_vector[i].data(), block_size * sizeof(datatype));
        }
        SZ::verify(uncompressed.data(), decompressed.data(), uncompressed.size());
//    for (size_t t = 0; t < testcase; t++) {
//        for (size_t i = 0; i < block_size; i++) {
//            auto idx = t * block_size + i;
//            if (fabs(uncompressed[idx] - decompressed[idx]) > error_bound) {
//                std::cout << t << " " << i << " " << uncompressed[idx] << " " << decompressed[idx] << std::endl;
//            }
//        }
//    }
    }
}
