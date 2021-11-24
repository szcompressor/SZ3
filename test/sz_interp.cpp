#include "sz_interp.hpp"


int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "SZ-interp v" << SZ_versionString() << std::endl;
        std::cout << "usage: " << argv[0] <<
                  " data_file -num_dim dim0 .. dimn relative_eb [interp_op interp_direction]"
                  << std::endl;
        std::cout << "example: " << argv[0] <<
                  " qmcpack.dat -3 33120 69 69 1e-3 [1 0]" << std::endl;
        return 0;
    }

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
            interp_lorenzo_tuning<float, 1>(argv[1], reb, false, dims[0]);
        } else if (dim == 2) {
            interp_lorenzo_tuning<float, 2>(argv[1], reb, false, dims[0], dims[1]);
        } else if (dim == 3) {
            interp_lorenzo_tuning<float, 3>(argv[1], reb, true, dims[0], dims[1], dims[2]);
        } else if (dim == 4) {
            interp_lorenzo_tuning<float, 4>(argv[1], reb, false, dims[0], dims[1], dims[2], dims[3]);
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

    auto data = SZ::readfile<float>(argv[1], num);
    double eb = reb * SZ::data_range(data.get(), num);


    if (dim == 1) {
        interp_compress_decompress<float, 1>(argv[1], data.get(), num, eb, interp_op, direction_op, block_size,
                                             interp_block_size, dims[0]);
    } else if (dim == 2) {
        interp_compress_decompress<float, 2>(argv[1], data.get(), num, eb, interp_op, direction_op, block_size,
                                             interp_block_size, dims[0], dims[1]);
    } else if (dim == 3) {
        interp_compress_decompress<float, 3>(argv[1], data.get(), num, eb, interp_op, direction_op, block_size,
                                             interp_block_size, dims[0], dims[1], dims[2]);
    } else if (dim == 4) {
        interp_compress_decompress<float, 4>(argv[1], data.release(), num, eb, interp_op,
                                             direction_op, block_size,
                                             interp_block_size, dims[0], dims[1], dims[2], dims[3]);
    }


    return 0;
}
