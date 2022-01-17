#include "sz_interp.hpp"

int main(int argc, char **argv) {


    int dim = atoi(argv[2] + 1);
    assert(1 <= dim && dim <= 4);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }
    double reb = atof(argv[argp++]);
    if (argp >= argc) {
        if (dim == 1) {
            interp_lorenzo_tuning<int, 1>(argv[1], reb, false, dims[0]);
        } else if (dim == 2) {
            interp_lorenzo_tuning<int, 2>(argv[1], reb, false, dims[0], dims[1]);
        } else if (dim == 3) {
            interp_lorenzo_tuning<int, 3>(argv[1], reb, true, dims[0], dims[1], dims[2]);
        } else if (dim == 4) {
            interp_lorenzo_tuning<int, 4>(argv[1], reb, false, dims[0], dims[1], dims[2], dims[3]);
        }
        return 0;
    }

    int interp_op = 0; // linear
    int direction_op = 0; // dimension high -> low
    if (argp < argc) {
        interp_op = atoi(argv[argp++]);
    }
    if (argp < argc) {
        direction_op = atoi(argv[argp++]);
    }

    int block_size = 6;
    int interp_block_size = 32;
    if (argp < argc) {
        block_size = atoi(argv[argp++]);
    }
    if (argp < argc) {
        interp_block_size = atoi(argv[argp++]);
    }
    size_t num = 0;

    auto data = SZ::readfile<int>(argv[1], num);
    double eb = reb * SZ::data_range(data.get(), num);


    if (dim == 1) {
        interp_compress_decompress<int, 1>(argv[1], data.get(), num, eb, interp_op, direction_op, block_size,
                                             interp_block_size, dims[0]);
    } else if (dim == 2) {
        interp_compress_decompress<int, 2>(argv[1], data.get(), num, eb, interp_op, direction_op, block_size,
                                             interp_block_size, dims[0], dims[1]);
    } else if (dim == 3) {
        interp_compress_decompress<int, 3>(argv[1], data.get(), num, eb, interp_op, direction_op, block_size,
                                             interp_block_size, dims[0], dims[1], dims[2]);
    } else if (dim == 4) {
        interp_compress_decompress<int, 4>(argv[1], data.release(), num, eb, interp_op,
                                             direction_op, block_size,
                                             interp_block_size, dims[0], dims[1], dims[2], dims[3]);
    }


    return 0;
}
