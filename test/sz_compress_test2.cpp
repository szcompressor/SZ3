#include <SZ.hpp>
#include <def.hpp>
#include <cstdio>
#include <iostream>
#include <memory>

int main(int argc, char **argv) {
    size_t num = 0;
    auto data = SZ::readfile<float>(argv[1], num);
    std::cout << "Read " << num << " elements\n";

    int r1 = atoi(argv[2]);
    int r2 = atoi(argv[3]);
    int r3 = atoi(argv[4]);
    float reb = atof(argv[5]);
    int block_size = atoi(argv[6]);
    int pred_dim = atoi(argv[7]);
    int lorenzo_op = atoi(argv[8]);
    int regression_op = atoi(argv[9]);

    int stride = block_size;

    float max = data[0];
    float min = data[0];
    for (int i = 1; i < num; i++) {
        if (max < data[i]) max = data[i];
        if (min > data[i]) min = data[i];
    }
    float eb = reb * (max - min);
    bool lossless = true;
    bool lorenzo_1 = lorenzo_op == 1 || lorenzo_op == 3;
    bool lorenzo_2 = lorenzo_op == 2 || lorenzo_op == 3;
    bool regression_1 = regression_op == 1 || regression_op == 3;
    bool regression_2 = regression_op == 2 || regression_op == 3;

    std::cout << "value range = " << max - min << std::endl;
//    std::cout << "abs error bound = " << eb << std::endl;

    auto ratio = SZ_Compress(data, eb, r1, r2, r3);
//    std::cerr << ratio;

    return 0;
}
