#ifndef _SZ_COMPOSED_PREDICTOR_HPP
#define _SZ_COMPOSED_PREDICTOR_HPP

#include "SZ3/def.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include <cassert>
#include <iostream>
#include <memory>

namespace SZ {

    template<class T, uint N>
    class ComposedPredictor : public concepts::PredictorInterface<T, N> {
    public:
        using block_iter = typename multi_dimensional_data<T, N>::block_iterator;

        ComposedPredictor(std::vector<std::shared_ptr<concepts::PredictorInterface<T, N>>> predictors) {
            this->predictors = predictors;
            predict_error.resize(predictors.size());
        }

        bool precompress(const block_iter &block) {
            std::vector<bool> isvalid(predictors.size());
            for (int i = 0; i < predictors.size(); i++) {
                isvalid[i] = predictors[i]->precompress(block);
                predict_error[i] = (isvalid[i] ? predictors[i]->est_error(block) : std::numeric_limits<double>::max());
            }
            sid = std::distance(predict_error.begin(), std::min_element(predict_error.begin(), predict_error.end()));
            return isvalid[sid];
        }

        inline void compress(const block_iter &iter, std::vector<int> &quant_inds) {
            selection.push_back(sid);
            predictors[sid]->compress(iter, quant_inds);
        }


        inline bool predecompress(const block_iter &block) {
            sid = selection[current_index++];
            return predictors[sid]->predecompress(block);
        }

        inline void decompress(const block_iter &block, int *&quant_inds_pos) {
            return predictors[sid]->decompress(block, quant_inds_pos);
        }


        void save(uchar *&c) const {
            auto tmp = c;
            for (const auto &p: predictors) {
                p->save(c);
            }
            // store selection

            *reinterpret_cast<size_t *>(c) = (size_t) selection.size();
            c += sizeof(size_t);
            if (selection.size()) {
                HuffmanEncoder<int> selection_encoder;
                selection_encoder.preprocess_encode(selection, 0);
                selection_encoder.save(c);
                selection_encoder.encode(selection, c);
                selection_encoder.postprocess_encode();
            }
//            *reinterpret_cast<size_t *>(c) = (size_t) selection.size();
//            c += sizeof(size_t);
//            memcpy(c, selection.data(), selection.size() * sizeof(int));
//            c += selection.size() * sizeof(int);
        }

        void load(const uchar *&c, size_t &remaining_length) {
            auto tmp = c;
            for (const auto &p: predictors) {
                p->load(c, remaining_length);
            }

            // load selection
            // TODO: check correctness
            size_t selection_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            if (selection_size > 0) {
                remaining_length -= sizeof(size_t);
                HuffmanEncoder<int> selection_encoder;
                selection_encoder.load(c, remaining_length);
                this->selection = selection_encoder.decode(c, selection_size);
                selection_encoder.postprocess_decode();
            }
//            size_t selection_size = *reinterpret_cast<const size_t *>(c);
//            c += sizeof(size_t);
//            this->selection = std::vector<int>(reinterpret_cast<const int *>(c),
//                                               reinterpret_cast<const int *>(c) + selection_size);
//            c += selection_size * sizeof(int);
        }


        T est_error(const block_iter &iter) {
            return predictors[sid]->est_error(iter);
        }

        inline size_t size_est() {
            size_t size_all = 0;
            for (const auto &p: predictors) {
                size_all += p->size_est();
            }
            return size_all;
        }

        void print() {
            std::vector<size_t> cnt(predictors.size(), 0);
            size_t cnt_total = 0;
            for (auto &sel: selection) {
                cnt[sel]++;
                cnt_total++;
            }
            for (int i = 0; i < predictors.size(); i++) {
//                predictors[i]->print();
                printf("Blocks:%ld, Percentage:%.2f\n", cnt[i], 1.0 * cnt[i] / cnt_total);
            }
        }

        inline size_t get_padding() {
            size_t max_padding = 0;
            for (const auto &p: predictors) {
                max_padding = std::max(max_padding, p->get_padding());
            }
            return max_padding;
        };

        void clear() {
            for (auto &pred: predictors) {
                pred->clear();
            }
            selection.clear();
        }

    private:
        std::vector<std::shared_ptr<concepts::PredictorInterface<T, N>>> predictors;
        std::vector<int> selection;
        int sid = 0;                            // selected index
        size_t current_index = 0;            // for decompression only
        std::vector<double> predict_error;

    };

}


#endif
