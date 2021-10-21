#include "utils/ByteUtil.h"
#include "utils/FileUtil.h"
#include "utils/Verification.hpp"
#include "def.hpp"
#include <cstdio>
#include <iostream>
#include <memory>


template<class T>
class PMCMean {
public:
    PMCMean(double _eb) : eb(_eb), sum(0), n(0),
                          dmax(std::numeric_limits<T>::min()),
                          dmin(std::numeric_limits<T>::max()) {}

    bool append(T d) {
//        return false;
        double avg = (sum + d) / (n + 1);
        T dmax_ = std::max(dmax, d);
        T dmin_ = std::min(dmin, d);
        if (fabs(dmax_ - avg) > eb || fabs(dmin_ - avg) > eb) {
            return false;
        }
        sum += d;
        n++;
        dmax = dmax_;
        dmin = dmin_;
        return true;
    }

    size_t get_size() {
        return 8;
    }

    void get_data(T *data, int len) {
        assert(len == n);
        for (int i = 0; i < n; i++) {
            data[i] = sum / n;
        }
    }

    float sum, eb;
    T dmax, dmin;
    int n;

};


template<class T>
class Swing {
public:
    class LinearFunction {
    public:
        LinearFunction() {};

        LinearFunction(int x1, T y1, int x2, T y2) {
            a = (y2 - y1) / (x2 - x1);
            b = y1 - a * x1;
        }

        T get(int x) {
            return a * x + b;
        }

        double a, b;
    };

    Swing(double _eb) : n(0), eb(_eb) {}

    bool append(T d) {
//        return false;
        if (n == 0) {
            d0 = d;
        } else if (n == 1) {
            ubound = LinearFunction(0, d0, 1, d + eb);
            lbound = LinearFunction(0, d0, 1, d - eb);
        } else {
            double uba = ubound.get(n), lba = lbound.get(n);
            if (uba + eb < d || lba - eb > d) {
                return false;
            } else {
                if (uba - eb > d) {
                    ubound = LinearFunction(0, d0, n, d + eb);
                }
                if (lba + eb < d) {
                    lbound = LinearFunction(0, d0, n, d - eb);
                }
            }
        }

        n++;
        return true;
    }

    size_t get_size() {
        return 8 + 4;
    }

    void get_data(T *data, int len) {
        assert(len == n);
        double a = (lbound.a + ubound.a) / 2.0, b = (lbound.b + ubound.b) / 2.0;
        for (int i = 0; i < n; i++) {
            data[i] = a * i + b;
        }
    }

    LinearFunction ubound, lbound;
    T d0;
    double eb;
    int n;

};


template<class T>
class Gorilla {
public:
    Gorilla() : sum(0), n(0), leadingzeros(std::numeric_limits<int>::max()), trailingzeros(0) {}

    bool append(T d) {
//        return false;
        if (n == 0) {
            bits += 4;
        } else {
            SZ::lfloat cur;
            cur.value = d;
            uint dxor = cur.ivalue ^ last.ivalue;
            if (dxor == 0) {
                bits += 1;
            } else {
                int lead = numberOfLeadingZeros(dxor), trail = numberOfTrailingZeros(dxor);
                if (lead >= 32) {
                    lead = 31;
                }
                bits += 1;
                if (lead >= leadingzeros && trail >= trailingzeros) {
                    bits += 1 + 32 - leadingzeros - trailingzeros;
                } else {
                    int bits_add = 1 + 5 + 6 + 32 - lead - trail;
//                    if (bits_add < 32) {
                    bits += bits_add;
                    leadingzeros = lead;
                    trailingzeros = trail;
//                    } else {
//                        return false;
//                    }
                }
            }
        }
        last.value = d;
        n++;
        return true;
    }

    size_t get_size() {
        return ceil(bits / 8.0);

    }

    void get_data(T *data) {
    }

    int numberOfTrailingZeros(int i) {
        i = ~i & (i - 1);
        if (i <= 0) return i & 32;
        int n = 1;
        if (i > 1 << 16) {
            n += 16;
            i >>= 16;
        }
        if (i > 1 << 8) {
            n += 8;
            i >>= 8;
        }
        if (i > 1 << 4) {
            n += 4;
            i >>= 4;
        }
        if (i > 1 << 2) {
            n += 2;
            i >>= 2;
        }
        return n + (i >> 1);
    }

    int numberOfLeadingZeros(int i) {
        if (i <= 0)
            return (i == 0 ? 32 : 0);
        int n = 31;
        if (i >= 1 << 16) {
            n -= 16;
            i >>= 16;
        }
        if (i >= 1 << 8) {
            n -= 8;
            i >>= 8;
        }
        if (i >= 1 << 4) {
            n -= 4;
            i >>= 4;
        }
        if (i >= 1 << 2) {
            n -= 2;
            i >>= 2;
        }
        return n - (i >> 1);
    }

    int leadingzeros, trailingzeros;
    size_t bits = 0;
    SZ::lfloat last;
    float sum;
    int n;

};

template<class T>
void compress(char *data_path, size_t r0, size_t r1, size_t batch_size, double reb) {

    size_t num;
    auto data = SZ::readfile<float>(data_path, num);
    std::cout << "Read " << num << " elements\n";
    assert(num == r1 * r0);
    printf("r1 = %lu , r0 = %lu \n", r1, r0);
    std::string src_file_name(data_path);

    float max = data[0];
    float min = data[0];
    for (int i = 1; i < num; i++) {
        if (max < data[i]) max = data[i];
        if (min > data[i]) min = data[i];
    }
    double eb = (max - min) * reb;
    printf("reb= %.5G, aeb= %.5G\n", reb, eb);

    std::vector<T> dec_data(num);
    size_t total_size = 0, total_pmc_n = 0, total_swing_n = 0, total_gorilla_n = 0;
    for (size_t t = 0; t < r1; t++) {
        T *d = &data[t * r0];
        T *dec = &dec_data[t * r0];
        for (size_t b = 0; b < r0;) {
            size_t br = (b + batch_size > r0 ? r0 : b + batch_size);
            size_t block_size = br - b, block_bytes = (br - b) * sizeof(T);
            double pmc_ratio = 0, swing_ratio = 0, gorilla_ratio = 0;
            size_t pmc_n = block_size, swing_n = 0, gorilla_n = 0;

            PMCMean<T> pmcMean(eb);
            Swing<T> swing(eb);
            Gorilla<T> gorilla;
            for (size_t i = b; i < br; i++) {
                if (!pmcMean.append(d[i])) {
                    pmc_ratio = (i - b) * sizeof(T) * 1.0 / pmcMean.get_size();
                    pmc_n = i - b;
                    break;
                }
            }
            if (pmc_n < block_size) {
                swing_n = block_size;
                for (size_t i = b; i < br; i++) {
                    if (!swing.append(d[i])) {
                        swing_ratio = (i - b) * sizeof(T) * 1.0 / swing.get_size();
                        swing_n = i - b;
                        break;
                    }
                }

                if (swing_n < block_size) {
                    gorilla_n = block_size;
                    for (size_t i = b; i < br; i++) {
                        if (!gorilla.append(d[i])) {
                            gorilla_ratio = (i - b) * sizeof(T) * 1.0 / gorilla.get_size();
                            gorilla_n = i - b;
                            break;
                        }
                    }
                    if (gorilla_n == 1) {
                        gorilla_n = 0;
                    }
                    if (gorilla.get_size() >= block_bytes) {
                        gorilla_n = 0;
                    }
                }
            }

            if (pmc_n == block_size ||
                (swing_n < block_size && pmc_n > 2 && pmc_ratio > swing_ratio && pmc_ratio > gorilla_ratio)) {
                pmcMean.get_data(dec + b, pmc_n);
                total_size += pmcMean.get_size();
                total_pmc_n += pmc_n;
                b += pmc_n;
            } else if (swing_n == block_size || (swing_n > 2 && swing_ratio > gorilla_ratio)) {
                swing.get_data(dec + b, swing_n);
                total_size += swing.get_size();
                total_swing_n += swing_n;
                b += swing_n;
            } else if (gorilla_n == block_size || (gorilla_n > 1)) {
                for (size_t i = b; i < b + gorilla_n; i++) {
                    dec[i] = d[i];
                }
                total_size += gorilla.get_size();
                total_gorilla_n += gorilla_n;
                b += gorilla_n;
            } else {
                //fallback
                total_size += block_bytes;
                for (size_t i = b; i < br; i++) {
                    dec[i] = d[i];
                }
                b += block_size;
            }
        }
    }
    double ratio = num * sizeof(T) * 1.0 / total_size;
    double max_diff, psnr, nrmse;
    SZ::verify<T>(data.get(), dec_data.data(), num, max_diff, psnr, nrmse);
    printf("method=mdb, file=%s, block=%lu, compression_ratio=%.3f, reb=%.1e, eb=%.6f, psnr=%.3f, nsmse=%e, compress_time=%.3f, decompress_time=%.3f, ori_bytes=%lu, compressed_bytes=%lu, pmc=%.2f%%, swing=%.2f%%, gorilla=%.2f%%\n\n",
           src_file_name.data(), batch_size, ratio, reb, max_diff, psnr, nrmse, 0.0, 0.0, num * sizeof(T), total_size,
           total_pmc_n * 100.0 / num, total_swing_n * 100.0 / num, total_gorilla_n * 100.0 / num);
}

int main(int argc, char **argv) {


    int dim = 2;
    int argp = 2;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }
    auto reb = atof(argv[argp++]);
    int batch_size = atoi(argv[argp++]);

    compress<float>(argv[1], dims[1], dims[0], batch_size, reb);

    return 0;
}