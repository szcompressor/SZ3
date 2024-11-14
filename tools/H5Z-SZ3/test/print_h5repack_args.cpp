#include <cerrno>
#include <cstdio>
#include <cstdlib>

#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"

inline void usage() {
    printf("Usage: print_h5repack_args <options>\n");
    printf("Options:\n");
    printf("	-c config file path\n");
    printf("* examples: \n");
    printf(
        "	print_h5repack_args -c sz3.conf (output: -f "
        "UD=32024,0,11,4077060608,16843009,0,33554432,4054449152,1348619730,16818239,257,2147483904,2147483648,"
        "16777216)\n");
    exit(0);
}

int main(int argc, char *argv[]) {
    if (argc == 1) {
        usage();
    }
    char *conPath = nullptr;
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] != '-' || argv[i][2]) {
            usage();
        }
        switch (argv[i][1]) {
            case 'c':
                if (++i == argc) usage();
                conPath = argv[i];
                break;
            default:
                usage();
                break;
        }
    }
    SZ3::Config conf;
    conf.loadcfg(conPath);
    std::vector<unsigned int> cd_values(std::ceil(conf.size_est() / 1.0 / sizeof(int)), 0);
    auto buffer = reinterpret_cast<unsigned char *>(cd_values.data());
    auto confSizeReal = conf.save(buffer);

    int cd_nelmts = std::ceil(confSizeReal / 1.0 / sizeof(int));
    // conf.print();
    if (cd_nelmts >= 20) {
        printf("h5repack can only take 20 cd_values, but got %d\n cd_values in SZ3", cd_nelmts);
        return 0;
    }
    printf("-f UD=32024,0,%d", cd_nelmts);
    for (int i = 0; i < cd_nelmts; i++) {
        printf(",%u", cd_values[i]);
    }
    printf("\n");
}
