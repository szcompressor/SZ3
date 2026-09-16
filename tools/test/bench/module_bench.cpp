/**
 * @file module_bench.cpp
 * @brief Measure a module composition's characteristics on a real dataset.
 *
 * Module characteristics -- compression ratio, achieved error, throughput -- cannot be authored,
 * only measured, and they depend on the data. This prints one CSV row per (composition, error
 * bound) so the numbers can be collected across datasets and compared over time.
 *
 * @code
 *   module_bench --input velocityx.d64 --type f64 --dims 256 384 384 \
 *                --eb 1e-2,1e-3,1e-4 --pipelines mgard-fused,mgard,sperr-fused,sperr
 * @endcode
 *
 * `--pipelines` names compositions, not `ALGO_*` values: the point is to compare a fused
 * implementation against the composable one built from the same modules.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "SZ3/api/sz_dev.hpp"

namespace {

using clk = std::chrono::steady_clock;

struct Result {
    size_t cmp_bytes = 0;
    double max_abs_err = 0;
    double compress_s = 0;
    double decompress_s = 0;
};

template <class T>
std::vector<T> readRaw(const std::string &path, size_t num) {
    std::vector<T> v(num);
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    f.read(reinterpret_cast<char *>(v.data()), static_cast<std::streamsize>(num * sizeof(T)));
    if (!f) throw std::runtime_error("short read from " + path + " (wrong dims or type?)");
    return v;
}

template <class T>
double valueRange(const std::vector<T> &v) {
    auto mm = std::minmax_element(v.begin(), v.end());
    return static_cast<double>(*mm.second) - static_cast<double>(*mm.first);
}

template <class T>
double maxAbsErr(const std::vector<T> &a, const std::vector<T> &b) {
    double m = 0;
    for (size_t i = 0; i < a.size(); i++) {
        m = std::max(m, std::fabs(static_cast<double>(a[i]) - static_cast<double>(b[i])));
    }
    return m;
}

/// Run one composition end to end. `mk` returns a fresh compressor; two are needed because a
/// composition carries state from compression that decompression must not see.
template <class T, class Make>
Result run(const SZ3::Config &conf, const std::vector<T> &original, Make mk) {
    std::vector<T> work = original;
    std::vector<SZ3::uchar> buf(conf.num * sizeof(T) * 2 + (1u << 20));

    auto c = mk();
    auto t0 = clk::now();
    const size_t n = c->compress(conf, work.data(), buf.data(), buf.size());
    auto t1 = clk::now();

    auto d = mk();
    std::vector<T> out(conf.num, T{0});
    auto t2 = clk::now();
    d->decompress(conf, buf.data(), n, out.data());
    auto t3 = clk::now();

    Result r;
    r.cmp_bytes = n;
    r.max_abs_err = maxAbsErr(original, out);
    r.compress_s = std::chrono::duration<double>(t1 - t0).count();
    r.decompress_s = std::chrono::duration<double>(t3 - t2).count();
    return r;
}

template <class T, SZ3::uint N>
bool runPipeline(const std::string &name, const SZ3::Config &conf, const std::vector<T> &data, double eb, Result &out) {
    const int radius = conf.quantbinCnt / 2;
    using Q = SZ3::LinearQuantizer<T>;

    if (name == "interp") {
        out = run(conf, data, [&] {
            return SZ3::make_compressor_sz_generic<T, N>(SZ3::InterpolationDecomposition<T, N, Q>(conf, Q(eb, radius)),
                                                         SZ3::HuffmanEncoder<int>(), SZ3::Lossless_zstd());
        });
        return true;
    }
    if (name == "lorenzo") {
        using P = SZ3::LorenzoPredictor<T, N, 1>;
        out = run(conf, data, [&] {
            return SZ3::make_compressor_sz_generic<T, N>(
                SZ3::BlockwiseDecomposition<T, N, P, Q>(conf, P(eb), Q(eb, radius)), SZ3::HuffmanEncoder<int>(),
                SZ3::Lossless_zstd());
        });
        return true;
    }
    if (name == "nopred") {
        out = run(conf, data, [&] {
            return SZ3::make_compressor_sz_generic<T, N>(SZ3::NoPredictionDecomposition<T, N, Q>(conf, Q(eb, radius)),
                                                         SZ3::HuffmanEncoder<int>(), SZ3::Lossless_zstd());
        });
        return true;
    }
    if constexpr (N == 1 || N == 2 || N == 3) {
        if (name == "mgard-fused") {
            out = run(conf, data, [&] {
                return SZ3::make_compressor_sz_generic<T, N>(SZ3::MGARDFusedDecomposition<T, N>(eb, radius),
                                                             SZ3::HuffmanEncoder<int>(), SZ3::Lossless_zstd());
            });
            return true;
        }
        if (name == "mgard") {
            out = run(conf, data, [&] {
                return SZ3::make_compressor_sz_generic<T, N>(
                    SZ3::MGARDDecomposition<T, N>(eb, [radius](double lvl) { return Q(lvl, radius); }),
                    SZ3::HuffmanEncoder<int>(), SZ3::Lossless_zstd());
            });
            return true;
        }
    }
    if constexpr (N == 3) {
        if (std::is_floating_point<T>::value) {
            const SZ3::SPERR::dims_type sd{conf.dims[2], conf.dims[1], conf.dims[0]};
            if (name == "sperr-fused") {
                out = run(conf, data, [&] {
                    return SZ3::make_compressor_sz_generic<T, N>(
                        SZ3::SPERRFusedDecomposition<T, N>(), SZ3::SPERREncoder<int64_t, N>(sd), SZ3::Lossless_zstd());
                });
                return true;
            }
            if (name == "sperr") {
                out = run(conf, data, [&] {
                    return SZ3::make_compressor_sz_generic<T, N>(
                        SZ3::SPERRDecomposition<T, N>(), SZ3::SPERREncoder<int64_t, N>(sd), SZ3::Lossless_zstd());
                });
                return true;
            }
        }
    }
    return false;
}

template <class T, SZ3::uint N>
int benchmark(const std::string &path, const std::vector<size_t> &dims, const std::vector<double> &rel_ebs,
              const std::vector<std::string> &pipelines) {
    SZ3::Config probe;
    probe.setDims(dims.begin(), dims.end());
    const auto data = readRaw<T>(path, probe.num);
    const double vr = valueRange(data);
    const size_t raw = probe.num * sizeof(T);

    printf("dataset,pipeline,rel_eb,abs_eb,raw_bytes,cmp_bytes,ratio,max_abs_err,err_over_eb,");
    printf("compress_MBps,decompress_MBps\n");

    for (double rel : rel_ebs) {
        const double eb = rel * vr;
        SZ3::Config conf;
        conf.setDims(dims.begin(), dims.end());
        conf.errorBoundMode = SZ3::EB_ABS;
        conf.absErrorBound = eb;

        for (const auto &name : pipelines) {
            Result r;
            if (!runPipeline<T, N>(name, conf, data, eb, r)) {
                fprintf(stderr, "skipping '%s': not available for this type/dimension\n", name.c_str());
                continue;
            }
            printf("%s,%s,%g,%g,%zu,%zu,%.4f,%.6g,%.6f,%.1f,%.1f\n", path.c_str(), name.c_str(), rel, eb, raw,
                   r.cmp_bytes, static_cast<double>(raw) / static_cast<double>(r.cmp_bytes), r.max_abs_err,
                   r.max_abs_err / eb, raw / r.compress_s / 1e6, raw / r.decompress_s / 1e6);
            fflush(stdout);
        }
    }
    return 0;
}

std::vector<std::string> split(const std::string &s, char sep) {
    std::vector<std::string> out;
    size_t start = 0;
    while (start <= s.size()) {
        const size_t at = s.find(sep, start);
        out.push_back(s.substr(start, at == std::string::npos ? std::string::npos : at - start));
        if (at == std::string::npos) break;
        start = at + 1;
    }
    return out;
}

void usage(const char *argv0) {
    fprintf(stderr,
            "usage: %s --input FILE --type f32|f64 --dims D0 [D1 [D2]]\n"
            "          [--eb 1e-2,1e-3,1e-4] [--pipelines interp,lorenzo,...]\n\n"
            "  --eb         relative to the value range; default 1e-2,1e-3,1e-4\n"
            "  --pipelines  interp, lorenzo, nopred, mgard-fused, mgard (1-3D),\n"
            "               sperr-fused, sperr (3D floating-point)\n\n"
            "Prints CSV to stdout, one row per (pipeline, error bound).\n",
            argv0);
}

}  // namespace

int main(int argc, char **argv) {
    std::string path, type = "f32";
    std::vector<size_t> dims;
    std::vector<double> rel_ebs{1e-2, 1e-3, 1e-4};
    std::vector<std::string> pipelines{"interp", "lorenzo"};

    for (int i = 1; i < argc; i++) {
        const std::string a = argv[i];
        if (a == "--input" && i + 1 < argc) {
            path = argv[++i];
        } else if (a == "--type" && i + 1 < argc) {
            type = argv[++i];
        } else if (a == "--dims") {
            while (i + 1 < argc && argv[i + 1][0] != '-') dims.push_back(strtoul(argv[++i], nullptr, 10));
        } else if (a == "--eb" && i + 1 < argc) {
            rel_ebs.clear();
            for (const auto &t : split(argv[++i], ',')) rel_ebs.push_back(strtod(t.c_str(), nullptr));
        } else if (a == "--pipelines" && i + 1 < argc) {
            pipelines = split(argv[++i], ',');
        } else {
            usage(argv[0]);
            return 1;
        }
    }
    if (path.empty() || dims.empty() || dims.size() > 3) {
        usage(argv[0]);
        return 1;
    }

    try {
        if (type == "f32") {
            if (dims.size() == 1) return benchmark<float, 1>(path, dims, rel_ebs, pipelines);
            if (dims.size() == 2) return benchmark<float, 2>(path, dims, rel_ebs, pipelines);
            return benchmark<float, 3>(path, dims, rel_ebs, pipelines);
        }
        if (type == "f64") {
            if (dims.size() == 1) return benchmark<double, 1>(path, dims, rel_ebs, pipelines);
            if (dims.size() == 2) return benchmark<double, 2>(path, dims, rel_ebs, pipelines);
            return benchmark<double, 3>(path, dims, rel_ebs, pipelines);
        }
        usage(argv[0]);
        return 1;
    } catch (const std::exception &e) {
        fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
