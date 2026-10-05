// One stream compressed with conf.openmp = true, decoded by a build of the other kind.
//
//   openmp_stream <with|without> write <prefix> [biomd]   compress, write <prefix>.sz3 and its decode <prefix>.dec
//   openmp_stream <with|without> read <prefix>            decode <prefix>.sz3, compare with <prefix>.dec bit for bit
//
// biomd writes ALGO_BIOMD water {frames, atoms, 3}, which the chunks split along the atoms, inside a molecule.
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
    if (argc != 4 && argc != 5) {
        fprintf(stderr, "usage: %s <with|without> <write|read> <prefix> [biomd]\n", argv[0]);
        return 2;
    }
    const std::string kind = argv[1], mode = argv[2], prefix = argv[3];
    const bool biomd = argc == 5 && std::string(argv[4]) == "biomd";
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
            // 1003 waters and one ion, 3010 atoms: 3 of the 4 chunks (at atoms 752, 1505, 2257) start inside a water
            SZ3::Config conf = biomd ? SZ3::Config(20, 3010, 3) : SZ3::Config(64, 32, 32);
            conf.cmprAlgo = biomd ? SZ3::ALGO_BIOMD : SZ3::ALGO_INTERP_LORENZO;
            conf.errorBoundMode = SZ3::EB_ABS;
            conf.absErrorBound = biomd ? 5e-4 : 1e-3;
            conf.openmp = true;
            std::vector<float> data(conf.num);
            for (size_t i = 0; i < conf.num; i++) {
                const double x = static_cast<double>(i);
                data[i] = static_cast<float>(std::sin(0.01 * x) + 0.001 * static_cast<double>(i % 97));
            }
            if (biomd) {
                // O, then H1 and H2 0.09572 nm from it at 104.52 degrees; each frame turns and moves every water
                const size_t atoms = conf.dims[1];
                for (size_t f = 0; f < conf.dims[0]; f++) {
                    for (size_t w = 0; w * 3 + 2 < atoms; w++) {
                        float* o = &data[(f * atoms + w * 3) * 3];
                        const double a = 0.1 * static_cast<double>(f) + 0.7 * static_cast<double>(w);
                        const double half = 0.5 * 104.52 * 3.14159265358979 / 180;
                        for (int c = 0; c < 3; c++)
                            o[c] = static_cast<float>(0.31 * static_cast<double>((w * (7 + c)) % 31) +
                                                      0.002 * static_cast<double>(f) * (c + 1));
                        for (int h = 0; h < 2; h++) {
                            const double t = a + (h ? half : -half);
                            o[3 + 3 * h] = static_cast<float>(o[0] + 0.09572 * std::cos(t));
                            o[4 + 3 * h] = static_cast<float>(o[1] + 0.09572 * std::sin(t));
                            o[5 + 3 * h] = o[2];
                        }
                    }
                }
            }
            size_t cmpSize = 0;
            char* cmpData = SZ_compress(conf, data.data(), cmpSize);
            std::vector<float> dec(conf.num);
            float* decPtr = dec.data();
            SZ3::Config readConf;
            SZ_decompress(readConf, cmpData, cmpSize, decPtr);
            for (size_t i = 0; i < conf.num; i++) {
                if (std::fabs(dec[i] - data[i]) > conf.absErrorBound) {
                    fprintf(stderr, "FAIL  value %zu comes back past the bound\n", i);
                    delete[] cmpData;
                    return 1;
                }
            }
            // The payload of a chunked stream starts, after the 16-byte header, with its chunk count.
            int32_t chunks = 0;
            if (readConf.openmp) {
                const SZ3::uchar* pos = reinterpret_cast<const SZ3::uchar*>(cmpData) + 16;
                SZ3::read(chunks, pos);
                // biomd: the first chunk keeps all frames, so it was split along the atoms
                SZ3::Config first;
                size_t remaining = cmpSize - 20;
                first.load(pos, remaining);
                if (biomd && first.dims[0] != conf.dims[0]) {
                    fprintf(stderr, "FAIL  the first chunk holds %zu of %zu frames\n", first.dims[0], conf.dims[0]);
                    delete[] cmpData;
                    return 1;
                }
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
