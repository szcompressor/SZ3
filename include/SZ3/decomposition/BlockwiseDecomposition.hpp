/**
 * @file BlockwiseDecomposition.hpp
 * @ingroup Decomposition
 */

#ifndef SZ3_BLOCKWISE_DECOMPOSITION_HPP
#define SZ3_BLOCKWISE_DECOMPOSITION_HPP

#include <cstring>

#include "Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/utils/BlockwiseIterator.hpp"
#include "SZ3/utils/Config.hpp"

namespace SZ3 {

/**
 * @brief Predictor-based block-wise decomposition (one of two primary decomposition strategies).
 * 
 * This is the **predictor-based** option for the `Decomposition` stage in `SZGenericCompressor`.
 * It splits data into fixed-size blocks and applies a `PredictorInterface` implementation
 * (e.g., Lorenzo, Regression) to each data point within the block, then quantizes the prediction error.
 * 
 * This is the appropriate choice when:
 * - You have or want to implement a scalar point predictor (`PredictorInterface`).
 * - You need block-local prediction strategies (e.g., switching between Lorenzo and Regression per block).
 * 
 * The alternative decomposition approach (e.g., `InterpolationDecomposition`) operates directly on
 * the global data array without using `PredictorInterface`.
 * 
 * @tparam T Data type
 * @tparam N Dimension
 * @tparam Predictor A class implementing `PredictorInterface<T, N>`
 * @tparam Quantizer A class implementing `QuantizerInterface<T, int>`
 */
template <class T, uint N, class Predictor, class Quantizer>
class BlockwiseDecomposition : public concepts::DecompositionInterface<T, int, N> {
   public:
    using Block_iter = typename block_data<T, N>::block_iterator;

    /**
     * @brief Construct a new Blockwise Decomposition object
     * 
     * @param conf Configuration
     * @param predictor Predictor instance
     * @param quantizer Quantizer instance
     */
    BlockwiseDecomposition(const Config &conf, Predictor predictor, Quantizer quantizer)
        : predictor(predictor), quantizer(quantizer), fallback_predictor(conf.absErrorBound) {
        static_assert(std::is_base_of<concepts::PredictorInterface<T, N>, Predictor>::value,
                      "must implement the Predictor interface");
    }

    std::vector<int> compress(const Config &conf, T *data) override {
        auto data_with_padding = std::make_shared<block_data<T, N>>(data, conf.dims, predictor.get_padding(), true);
        auto block = data_with_padding->block_iter(conf.blockSize);
        std::vector<int> quant_inds;
        quant_inds.reserve(conf.num);
        do {
            concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;
            if (!predictor.precompress(block)) {
                predictor_withfallback = &fallback_predictor;
            }
            predictor_withfallback->precompress_block_commit();
            Block_iter::foreach (block, [&](T *c, const std::array<size_t, N> &index) {
                T pred = predictor_withfallback->predict(block, c, index);
                quant_inds.push_back(quantizer.quantize_and_overwrite(*c, pred));
            });

        } while (block.next());
        return quant_inds;
    }

    T *decompress(const Config &conf, std::vector<int> &quant_inds, T *dec_data) override {
        int *quant_inds_pos = &quant_inds[0];

        auto data_with_padding =
            std::make_shared<block_data<T, N>>(dec_data, conf.dims, predictor.get_padding(), false);
        auto block = data_with_padding->block_iter(conf.blockSize);
        do {
            concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;
            if (!predictor.predecompress(block)) {
                predictor_withfallback = &fallback_predictor;
            }
            Block_iter::foreach (block, [&](T *c, const std::array<size_t, N> &index) {
                T pred = predictor_withfallback->predict(block, c, index);
                *c = quantizer.recover(pred, *(quant_inds_pos++));
            });

        } while (block.next());

        return dec_data;
    }

    void save(uchar *&c) override {
        fallback_predictor.save(c);
        predictor.save(c);
        quantizer.save(c);
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        fallback_predictor.load(c, remaining_length);
        predictor.load(c, remaining_length);
        quantizer.load(c, remaining_length);
    }

    std::pair<int, int> get_out_range() override { return quantizer.get_out_range(); }

   private:
    Predictor predictor;
    Quantizer quantizer;
    LorenzoPredictor<T, N, 1> fallback_predictor;
};

template <class T, uint N, class Predictor, class Quantizer>
BlockwiseDecomposition<T, N, Predictor, Quantizer> make_decomposition_blockwise(const Config &conf, Predictor predictor,
                                                                                Quantizer quantizer) {
    return BlockwiseDecomposition<T, N, Predictor, Quantizer>(conf, predictor, quantizer);
}

}  // namespace SZ3
#endif
