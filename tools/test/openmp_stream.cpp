// One stream compressed with conf.openmp = true, decoded by a build of the other kind.
//
//   openmp_stream <with|without> write <prefix>   compress, write <prefix>.sz3 and its decode <prefix>.dec
//   openmp_stream <with|without> read <prefix>    decode <prefix>.sz3, compare with <prefix>.dec bit for bit
//
// tools/test/CMakeLists.txt builds this twice, with OpenMP and without, and has the build without it
// read what the other wrote. The first argument says which build this is meant to be, so a build that got OpenMP
// when it should not have (or the reverse) fails instead of testing nothing.

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "SZ3/api/sz.hpp"

static std::vector<char> slurp(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    return std::vector<char>(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

static void spill(const std::string& path, const char* data, size_t size) {
    std::ofstream(path, std::ios::binary).write(data, static_cast<std::streamsize>(size));
}

int main(int argc, char** argv) {
    if (argc != 4) {
        fprintf(stderr, "usage: %s <with|without> <write|read> <prefix>\n", argv[0]);
        return 2;
    }
    const std::string kind = argv[1], mode = argv[2], prefix = argv[3];
#ifdef _OPENMP
    const std::string built = "with";
#else
    const std::string built = "without";
#endif
    if (kind != built) {
        fprintf(stderr, "FAIL  this binary was built %s OpenMP, not %s it\n", built.c_str(), kind.c_str());
        return 1;
    }

    try {
        if (mode == "write") {
            SZ3::Config conf(64, 32, 32);
            conf.cmprAlgo = SZ3::ALGO_INTERP_LORENZO;
            conf.errorBoundMode = SZ3::EB_ABS;
            conf.absErrorBound = 1e-3;
            conf.openmp = true;
            std::vector<float> data(conf.num);
            for (size_t i = 0; i < conf.num; i++) {
                const double x = static_cast<double>(i);
                data[i] = static_cast<float>(std::sin(0.01 * x) + 0.001 * static_cast<double>(i % 97));
            }
            size_t cmpSize = 0;
            char* cmpData = SZ_compress(conf, data.data(), cmpSize);
            std::vector<float> dec(conf.num);
            float* decPtr = dec.data();
            SZ3::Config readConf;
            SZ_decompress(readConf, cmpData, cmpSize, decPtr);
            // The payload of a chunked stream starts, after the 16-byte header, with its chunk count.
            int32_t chunks = 0;
            if (readConf.openmp) {
                const SZ3::uchar* pos = reinterpret_cast<const SZ3::uchar*>(cmpData) + 16;
                SZ3::read(chunks, pos);
            }
            printf("wrote %s.sz3: %zu bytes, openmp flag %d, %d chunks\n", prefix.c_str(), cmpSize,
                   static_cast<int>(readConf.openmp), static_cast<int>(chunks));
#ifdef _OPENMP
            // One chunk would decode the same through either path and prove nothing.
            if (!readConf.openmp || chunks < 2) {
                fprintf(stderr, "FAIL  an OpenMP build wrote no multi-chunk stream; set OMP_NUM_THREADS above 1\n");
                delete[] cmpData;
                return 1;
            }
#endif
            spill(prefix + ".sz3", cmpData, cmpSize);
            spill(prefix + ".dec", reinterpret_cast<const char*>(dec.data()), dec.size() * sizeof(float));
            delete[] cmpData;
        } else if (mode == "read") {
            std::vector<char> cmp = slurp(prefix + ".sz3");
            std::vector<char> want = slurp(prefix + ".dec");
            if (cmp.empty() || want.empty()) {
                fprintf(stderr, "FAIL  nothing to read at %s.sz3 and %s.dec\n", prefix.c_str(), prefix.c_str());
                return 1;
            }
            SZ3::Config conf;
            float* dec = nullptr;
            SZ_decompress(conf, cmp.data(), cmp.size(), dec);
            const size_t got = conf.num * sizeof(float);
            const bool same = got == want.size() && memcmp(dec, want.data(), got) == 0;
            delete[] dec;
            if (!same) {
                fprintf(stderr, "FAIL  %s.sz3 decodes to other bits here than where it was written\n",
                        prefix.c_str());
                return 1;
            }
            printf("read %s.sz3: openmp flag %d, %zu bytes identical\n", prefix.c_str(),
                   static_cast<int>(conf.openmp), got);
        } else {
            fprintf(stderr, "unknown mode %s\n", mode.c_str());
            return 2;
        }
    } catch (const std::exception& e) {
        fprintf(stderr, "FAIL  %s\n", e.what());
        return 1;
    }
    return 0;
}
