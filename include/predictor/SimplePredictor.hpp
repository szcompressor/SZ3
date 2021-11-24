#ifndef _SZ_SIMPLE_PREDICTOR_HPP
#define _SZ_SIMPLE_PREDICTOR_HPP

#include "def.hpp"
#include "predictor/Predictor.hpp"
#include "utils/Iterator.hpp"
#include <cassert>

namespace SZ {

    // 1D 1-layer lorenzo predictor, aka simple predictor
    // Use previous data as predicted value
    template<class T, uint N>
    class SimplePredictor : public concepts::PredictorInterface<T, N> {
    public:
        static const uint8_t predictor_id = 0b00000001;
        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        SimplePredictor() {
            this->noise = 0;
        }

        SimplePredictor(double eb) {
            this->noise = this->noise = 0.5 * eb;
        }

        void precompress_data(const iterator &) const {}

        void postcompress_data(const iterator &) const {}

        void predecompress_data(const iterator &) const {}

        void postdecompress_data(const iterator &) const {}

        bool precompress_block(const std::shared_ptr<Range> &) { return true; }

        void precompress_block_commit() noexcept {}

        bool predecompress_block(const std::shared_ptr<Range> &) { return true; }

        /*
         * save doesn't need to store anything except the id
         */
        // std::string save() const {
        //   return std::string(1, predictor_id);
        // }
        void save(uchar *&c) const {
//            std::cout << "save Simple predictor" << std::endl;
            c[0] = predictor_id;
            c += sizeof(uint8_t);
        }

        /*
         * just verifies the ID, increments
         */
        // static SimplePredictor<T,N> load(const unsigned char*& c, size_t& remaining_length) {
        //   assert(remaining_length > sizeof(uint8_t));
        //   c += 1;
        //   remaining_length -= sizeof(uint8_t);
        //   return SimplePredictor<T,N>{};
        // }
        void load(const uchar *&c, size_t &remaining_length) {
//            std::cout << "load Simple predictor" << std::endl;
            c += sizeof(uint8_t);
            remaining_length -= sizeof(uint8_t);
        }

        void print() const {
            std::cout << "Simple predictor, noise = " << noise << "\n";
        }

        inline T estimate_error(const iterator &iter) const noexcept {
            return fabs(*iter - predict(iter)) + this->noise;
        }

        inline T predict(const iterator &iter) const noexcept {
            return do_predict(iter);
        }

        void clear() {}

    protected:
        T noise = 0;

    private:
        template<uint NN = N>
        inline typename std::enable_if<NN == 1, T>::type do_predict(const iterator &iter) const noexcept {
            return iter.prev(1);
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 2, T>::type do_predict(const iterator &iter) const noexcept {
            return iter.prev(0, 1);
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 3, T>::type do_predict(const iterator &iter) const noexcept {
            return iter.prev(0, 0, 1);

        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 4, T>::type do_predict(const iterator &iter) const noexcept {
            return iter.prev(0, 0, 0, 1);
        }
    };
}
#endif
