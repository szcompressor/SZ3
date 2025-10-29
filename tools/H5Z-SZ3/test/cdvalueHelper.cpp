#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"

inline void usage() {
    printf("Usage: print_h5repack_args <options>\n");
    printf("Options:\n");
    printf("\t-c config_file_path\n\t\t(read config file and print h5repack arguments to stdout)\n");
    printf("\t-r 'cd_values_string' -o output_config_file_path\n\t\t(read h5repack arguments and write config file)\n");
    printf("* examples: \n");
    printf(
        "\tprint_h5repack_args -c sz3.conf\n"
        "\t (output: -f UD=32024,0,11,4077060608,16843009,0,33554432,4054449152,1348619730,16818239,257,2147483904,2147483648,16777216)\n");
    printf(
        "\tprint_h5repack_args -r '32024,0,11,4077060608,16843009,0,33554432,4054449152,1348619730,16818239,257,2147483904,2147483648,16777216' -o sz3.conf\n");
    exit(0);
}

void parse_cd_values(const std::string& s, std::vector<unsigned int>& values) {
    std::stringstream ss(s);
    std::string item;
    // The format is filter_id,flags,nelmts,val1,val2,...
    // The config data is in val1, val2, ...
    // So we skip the first 3 items.
    std::getline(ss, item, ','); // filter id
    std::getline(ss, item, ','); // flags
    std::getline(ss, item, ','); // nelmts
    while (std::getline(ss, item, ',')) {
        values.push_back(std::stoul(item));
    }
}


int main(int argc, char* argv[]) {
    if (argc == 1) {
        usage();
    }
    char* conPath = nullptr;
    char* cdValuesStr = nullptr;
    char* outPath = nullptr;
    bool reverse_mode = false;

    for (int i = 1; i < argc; i++) {
        if (argv[i][0] != '-') {
            usage();
        }
        switch (argv[i][1]) {
            case 'c':
                if (++i == argc) usage();
                conPath = argv[i];
                break;
            case 'r':
                reverse_mode = true;
                if (++i == argc) usage();
                cdValuesStr = argv[i];
                break;
            case 'o':
                if (++i == argc) usage();
                outPath = argv[i];
                break;
            default:
                usage();
                break;
        }
    }

    if (reverse_mode) {
        if (!cdValuesStr || !outPath) {
            usage();
        }
        std::vector<unsigned int> cd_values;
        std::string cd_values_str(cdValuesStr);

        size_t pos = cd_values_str.find("UD=");
        if (pos != std::string::npos) {
            cd_values_str = cd_values_str.substr(pos + 3);
        }

        parse_cd_values(cd_values_str, cd_values);

        SZ3::Config conf;
        auto buffer = reinterpret_cast<const unsigned char*>(cd_values.data());
        conf.load(buffer);

        std::ofstream file(outPath);
        if (!file.is_open()) {
            printf("Error: cannot open output file %s\n", outPath);
            return 1;
        }
        file << conf.save_ini();
        file.close();
        printf("Config file saved to %s\n", outPath);
    } else {
        if (!conPath) {
            usage();
        }
        SZ3::Config conf;
        conf.loadcfg(conPath);
        std::vector<unsigned int> cd_values(std::ceil(conf.size_est() / 1.0 / sizeof(int)), 0);
        auto buffer = reinterpret_cast<unsigned char*>(cd_values.data());
        auto confSizeReal = conf.save(buffer);

        int cd_nelmts = std::ceil(confSizeReal / 1.0 / sizeof(int));
        // conf.print();
        if (cd_nelmts >= 20) {
            printf("h5repack can only take 20 cd_values, but got %d\n cd_values in SZ3", cd_nelmts);
            return 0;
        }
        printf("-f UD=32024,0,%d", cd_nelmts);
        for (int i = 0; i < cd_nelmts; i++) {
#if SZ3_BIG_ENDIAN
            printf(",%u", SZ3::byteswap(cd_values[i]));
#else
            printf(",%u", cd_values[i]);
#endif
        }
        printf("\n");
    }

    return 0;
}