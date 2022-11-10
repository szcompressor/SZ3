#include <vector>

#include <mdz.hpp>


int main(int argc, char **argv) {


    int dim = atoi(argv[2] + 1);
    assert(1 <= dim && dim <= 2);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }

    SZ::Config conf({1, dims[0]});
    if (dim == 2) {
        conf = SZ::Config({dims[0], dims[1]});
    }
    std::string input_path = argv[1];

    char *eb_op = argv[argp++] + 1;
    if (*eb_op == 'a') {
        conf.errorBoundMode = SZ::EB_ABS;
        conf.absErrorBound = atof(argv[argp++]);
    } else {
        conf.errorBoundMode = SZ::EB_REL;
        conf.relErrorBound = atof(argv[argp++]);
    }
    size_t batch_size = 0;

    if (argp < argc) {
        batch_size = atoi(argv[argp++]);
    }
    int method = 9;
    method_batch = 50; //method_batch >0 indicates ADP
    if (argp < argc) {
        int tmp = atoi(argv[argp++]);
        if (tmp <= 0) {
            method = -tmp;
            method_batch = 0;
        } else {
            method_batch = tmp;
        }
    }

    conf.blockSize = 128;
    conf.stride = 128;
    conf.quantbinCnt = 1024;
//    conf.enable_regression = false;
//    conf.quant_state_num = 4096;
    if (argp < argc) {
        conf.quantbinCnt = atoi(argv[argp++]);
    }
    MDZ_Compress<float, 2>(conf, method, batch_size, input_path);

}