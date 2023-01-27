#include <SZ3/predictor/Predictor.hpp>
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/api/sz.hpp"
#include <cstring>
#include <cmath>

using namespace SZ;

template<class T, uint N, class Quantizer, class Encoder, class Lossless>
class SZInterpolationEstimator {
public:
    SZInterpolationEstimator(Quantizer quantizer, Encoder encoder, Lossless lossless) :
            quantizer(quantizer), encoder(encoder), lossless(lossless) {

        static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                      "must implement the quatizer interface");
        static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                      "must implement the encoder interface");
        static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                      "must implement the lossless interface");
    }

    double estimate(const Config &conf, T *data, double &cr) {

        dimension_offsets[N - 1] = 1;
        for (int i = N - 2; i >= 0; i--) {
            dimension_offsets[i] = dimension_offsets[i + 1] * conf.dims[i + 1];
        }

        dimension_sequences = std::vector<std::array<int, N>>();
        auto sequence = std::array<int, N>();
        for (int i = 0; i < N; i++) {
            sequence[i] = i;
        }
        do {
            dimension_sequences.push_back(sequence);
        } while (std::next_permutation(sequence.begin(), sequence.end()));


        double pred_error = 0;

        std::array<size_t, N> begin_idx, end_idx;
        for (int i = 0; i < N; i++) {
            begin_idx[i] = 0;
            end_idx[i] = conf.dims[i] - 1;
        }

        quant_inds.reserve(conf.num);
        quant_inds.clear();
        pred_error += block_interpolation(data, begin_idx, end_idx, PB_predict_overwrite,
                                          interpolators[conf.interpAlgo], conf.interpDirection, 1);

        encoder.preprocess_encode(quant_inds, 0);
        size_t bufferSize = 1.2 * (quantizer.size_est() + encoder.size_est() + sizeof(T) * quant_inds.size());

        uchar *buffer = new uchar[bufferSize];
        uchar *buffer_pos = buffer;

        quantizer.save(buffer_pos);
        quantizer.postcompress_data();

        encoder.save(buffer_pos);
        encoder.encode(quant_inds, buffer_pos);
        encoder.postprocess_encode();
//            timer.stop("Coding");
        assert(buffer_pos - buffer < bufferSize);

        size_t compressed_size = 0;
        uchar *lossless_data = lossless.compress(buffer,
                                                 buffer_pos - buffer,
                                                 compressed_size);
        lossless.postcompress_data(buffer);

        printf("est per = %.5f, avg pred error = %.5G \n", quant_inds.size() * 1.0 / conf.num, pred_error / quant_inds.size());
        cr = quant_inds.size() * sizeof(T) * 1.0 / compressed_size;
        return pred_error;
    }

private:

    enum PredictorBehavior {
        PB_predict_overwrite, PB_predict, PB_recover
    };


    inline double quantize(size_t idx, T &d, T pred) {
        quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
        return fabs(d - pred);
    }


    double block_interpolation_1d(T *data, size_t begin, size_t end, size_t stride,
                                  const std::string &interp_func,
                                  const PredictorBehavior pb) {
        size_t n = (end - begin) / stride + 1;
        if (n <= 1) {
            return 0;
        }
        double predict_error = 0;

        size_t stride3x = 3 * stride;
        size_t stride5x = 5 * stride;
        if (interp_func == "linear" || n < 5) {
            for (size_t i = 1; i + 1 < n; i += 2) {
                T *d = data + begin + i * stride;
                predict_error += quantize(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
            }
            if (n % 2 == 0) {
                T *d = data + begin + (n - 1) * stride;
                if (n < 4) {
                    predict_error += quantize(d - data, *d, *(d - stride));
                } else {
                    predict_error += quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                }
            }
        } else {

            T *d;
            size_t i;
            for (i = 3; i + 3 < n; i += 2) {
                d = data + begin + i * stride;
                predict_error += quantize(d - data, *d,
                                          interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
            }
            d = data + begin + stride;
            predict_error += quantize(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

            d = data + begin + i * stride;
            predict_error += quantize(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
            if (n % 2 == 0) {
                d = data + begin + (n - 1) * stride;
                predict_error += quantize(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
            }

        }

        return predict_error;
    }

    template<uint NN = N>
    typename std::enable_if<NN == 1, double>::type
    block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                        const std::string &interp_func, const int direction, size_t stride = 1) {
        return block_interpolation_1d(data, begin[0], end[0], stride, interp_func, pb);
    }

    template<uint NN = N>
    typename std::enable_if<NN == 2, double>::type
    block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                        const std::string &interp_func, const int direction, size_t stride = 1) {
        double predict_error = 0;
        size_t stride2x = stride * 2;
        const std::array<int, N> dims = dimension_sequences[direction];
        for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
            size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]];
            predict_error += block_interpolation_1d(data, begin_offset,
                                                    begin_offset +
                                                    (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
                                                    stride * dimension_offsets[dims[0]], interp_func, pb);
        }
        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
            size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]];
            predict_error += block_interpolation_1d(data, begin_offset,
                                                    begin_offset +
                                                    (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
                                                    stride * dimension_offsets[dims[1]], interp_func, pb);
        }
        return predict_error;
    }

    template<uint NN = N>
    typename std::enable_if<NN == 3, double>::type
    block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                        const std::string &interp_func, const int direction, size_t stride = 1) {
        double predict_error = 0;
        size_t stride2x = stride * 2;
        const std::array<int, N> dims = dimension_sequences[direction];
        for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
            for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                      k * dimension_offsets[dims[2]];
                predict_error += block_interpolation_1d(data, begin_offset,
                                                        begin_offset +
                                                        (end[dims[0]] - begin[dims[0]]) *
                                                        dimension_offsets[dims[0]],
                                                        stride * dimension_offsets[dims[0]], interp_func, pb);
            }
        }
        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
            for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                      k * dimension_offsets[dims[2]];
                predict_error += block_interpolation_1d(data, begin_offset,
                                                        begin_offset +
                                                        (end[dims[1]] - begin[dims[1]]) *
                                                        dimension_offsets[dims[1]],
                                                        stride * dimension_offsets[dims[1]], interp_func, pb);
            }
        }
        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                      begin[dims[2]] * dimension_offsets[dims[2]];
                predict_error += block_interpolation_1d(data, begin_offset,
                                                        begin_offset +
                                                        (end[dims[2]] - begin[dims[2]]) *
                                                        dimension_offsets[dims[2]],
                                                        stride * dimension_offsets[dims[2]], interp_func, pb);
            }
        }
        return predict_error;
    }


    template<uint NN = N>
    typename std::enable_if<NN == 4, double>::type
    block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                        const std::string &interp_func, const int direction, size_t stride = 1) {
        double predict_error = 0;
        size_t stride2x = stride * 2;
        const std::array<int, N> dims = dimension_sequences[direction];
        for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
            for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                     t <= end[dims[3]]; t += stride2x) {
                    size_t begin_offset =
                            begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                            k * dimension_offsets[dims[2]] +
                            t * dimension_offsets[dims[3]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[0]] - begin[dims[0]]) *
                                                            dimension_offsets[dims[0]],
                                                            stride * dimension_offsets[dims[0]], interp_func, pb);
                }
            }
        }
        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
            for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                     t <= end[dims[3]]; t += stride2x) {
                    size_t begin_offset =
                            i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                            k * dimension_offsets[dims[2]] +
                            t * dimension_offsets[dims[3]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[1]] - begin[dims[1]]) *
                                                            dimension_offsets[dims[1]],
                                                            stride * dimension_offsets[dims[1]], interp_func, pb);
                }
            }
        }
        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                     t <= end[dims[3]]; t += stride2x) {
                    size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                          begin[dims[2]] * dimension_offsets[dims[2]] +
                                          t * dimension_offsets[dims[3]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[2]] - begin[dims[2]]) *
                                                            dimension_offsets[dims[2]],
                                                            stride * dimension_offsets[dims[2]], interp_func, pb);
                }
            }
        }

        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride : 0); k <= end[dims[2]]; k += stride) {
                    size_t begin_offset =
                            i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                            k * dimension_offsets[dims[2]] +
                            begin[dims[3]] * dimension_offsets[dims[3]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[3]] - begin[dims[3]]) *
                                                            dimension_offsets[dims[3]],
                                                            stride * dimension_offsets[dims[3]], interp_func, pb);
                }
            }
        }
        return predict_error;
    }

    std::vector<std::string> interpolators = {"linear", "cubic"};
    std::array<size_t, N> dimension_offsets;
    std::vector<std::array<int, N>> dimension_sequences;
    Quantizer quantizer;
    Encoder encoder;
    Lossless lossless;
    std::vector<int> quant_inds;
};

template<class T, uint N>
void estimate_compress(Config conf, T *data, double abs) {
    conf.cmprAlgo = ALGO_INTERP;
    conf.interpAlgo = INTERP_ALGO_CUBIC;
    conf.interpDirection = 0;
    conf.errorBoundMode = SZ::EB_ABS;
    conf.absErrorBound = abs;


    Timer timer(true);
    double CR = 0;
    SZInterpolationEstimator<T, N, SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd> estimator(
            SZ::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd());
    double pred_error = estimator.estimate(conf, data, CR);
    timer.stop("estimator");

    timer.start();
    size_t cmpSize = 0;
    auto cmpdata = SZ_compress(conf, data, cmpSize);
    timer.stop("compressor");
    delete[] cmpdata;

    printf("Est CR=%.5f, Real CR = %.3f\n", CR, conf.num * sizeof(T) * 1.0 / cmpSize);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cout << "usage: " << argv[0] << " data_file dim dim0 .. dimn" << std::endl;
        std::cout << "example: " << argv[0] << " qmcpack.dat 3 33120 69 69" << std::endl;
        return 0;
    }

    size_t num = 0;
    auto data = SZ::readfile<float>(argv[1], num);

    int dim = atoi(argv[2]);
    if (dim == 1) {
        Config config(atoi(argv[3]));
        estimate_compress<float, 1>(config, data.get(), atof(argv[4]));
    } else if (dim == 2) {
        Config config(atoi(argv[3]), atoi(argv[4]));
        estimate_compress<float, 2>(config, data.get(), atof(argv[5]));
    } else if (dim == 3) {
        Config config(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
        estimate_compress<float, 3>(config, data.get(), atof(argv[6]));
    } else if (dim == 4) {
        Config config(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
        estimate_compress<float, 4>(config, data.get(), atof(argv[7]));
    }
}