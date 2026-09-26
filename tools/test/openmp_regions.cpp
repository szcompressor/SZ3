// conf.openmp = true inside a program that runs its own OpenMP.
//
//   openmp_regions throw    one chunk throws: the call must throw, not hang or terminate
//   openmp_regions nested   compress from each thread of the program's own region; with nested
//                           parallelism on and OMP_THREAD_LIMIT set, regions get fewer threads
//                           than SZ3 asked for, and every stream must still decode inside the bound

#include <omp.h>

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "SZ3/api/sz.hpp"

int main(int argc, char** argv) {
    const std::string mode = argc > 1 ? argv[1] : "";
    if (mode == "throw") {
        // 5 rows over 4 threads: the chunk with 2 rows stays 4D, which ALGO_BIOMD refuses.
        SZ3::Config conf(5, 8, 8, 8);
        conf.cmprAlgo = SZ3::ALGO_BIOMD;
        conf.absErrorBound = 1e-3;
        conf.openmp = true;
        std::vector<float> data(conf.num, 1.0f);
        omp_set_num_threads(4);
        try {
            size_t cmpSize = 0;
            delete[] SZ_compress(conf, data.data(), cmpSize);
        } catch (const std::invalid_argument& e) {
            printf("PASS  threw: %s\n", e.what());
            return 0;
        }
        fprintf(stderr, "FAIL  a chunk that cannot be compressed did not throw\n");
        return 1;
    }
    if (mode == "nested") {
        int failures = 0;
#pragma omp parallel num_threads(4) reduction(+ : failures)
        for (int round = 0; round < 250; round++) {
            try {
                SZ3::Config conf(32, 64, 64);
                conf.absErrorBound = 1e-3;
                conf.openmp = true;
                std::vector<float> data(conf.num);
                for (size_t i = 0; i < conf.num; i++) data[i] = static_cast<float>(10 * std::sin(0.001 * i));
                size_t cmpSize = 0;
                std::unique_ptr<char[]> cmpData(SZ_compress(conf, data.data(), cmpSize));
                SZ3::Config readConf;
                std::unique_ptr<float[]> dec(SZ_decompress<float>(readConf, cmpData.get(), cmpSize));
                for (size_t i = 0; i < conf.num; i++) {
                    if (std::fabs(static_cast<double>(dec[i]) - data[i]) > conf.absErrorBound * (1 + 1e-6)) {
                        failures++;
                        break;
                    }
                }
            } catch (const std::exception& e) {
                failures++;
            }
        }
        if (failures) {
            fprintf(stderr, "FAIL  %d round trips threw or left the bound\n", failures);
            return 1;
        }
        printf("PASS  1000 round trips inside the bound\n");
        return 0;
    }
    fprintf(stderr, "usage: %s <throw|nested>\n", argv[0]);
    return 2;
}
