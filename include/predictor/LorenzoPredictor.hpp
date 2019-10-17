#ifndef _SZ_LORENZO_PREDICTOR_HPP
#define _SZ_LORENZO_PREDICTOR_HPP

#include <cassert>
#include "def.hpp"
#include "utils/Iterator.hpp"

namespace SZ{

  // N-d lorenzo predictor
  template <class T, uint N>
  class LorenzoPredictor;

  namespace {
    template <class T, uint N>
    class LorenzoBase {
      public:
      static const uint8_t predictor_id = 0b00000001;
      using iterator = typename multi_dimensional_range<T, N>::iterator;
      void precompress_data(const iterator) {}
      void postcompress_data(const iterator) {}
      void precompress_block(const iterator) {}
      void predecompress_data(const iterator) {}
      void postdecompress_data(const iterator) {}
      void predecompress_block(const iterator) {}

      /*
       * save doesn't need to store anything except the id
       */
      std::string save() const {
        return std::string(1, predictor_id);
      }

      /*
       * just verifies the ID, increments
       */
      static LorenzoPredictor<T,N> load(const unsigned char*& c, size_t& remaining_length) {
        assert(remaining_length > sizeof(uint8_t));
        c += 1;
        remaining_length -= sizeof(uint8_t);
        return LorenzoPredictor<T,N>{};
      }
    };
  }


  template <class T>
  class LorenzoPredictor<T, 1> : public LorenzoBase<T,1> {
  public:
    using iterator = typename LorenzoBase<T,1>::iterator;
    inline T predict(iterator iter) const noexcept { return iter.prev(0); };
  };

  template <class T>
  class LorenzoPredictor<T, 2>: public LorenzoBase<T,2> {
  public:
    using iterator = typename LorenzoBase<T,2>::iterator;
    inline T predict(iterator iter) const noexcept{
      return iter.prev(0, 1) + iter.prev(1, 0) - iter.prev(1, 1);
    };
  };

  template <class T>
  class LorenzoPredictor<T, 3>: public LorenzoBase<T,3> {
  public:
    using iterator = typename LorenzoBase<T,3>::iterator;
    inline T predict(const iterator iter) const noexcept{
      return iter.prev(0, 0, 1) + iter.prev(0, 1, 0) + iter.prev(1, 0, 0) 
          - iter.prev(0, 1, 1) - iter.prev(1, 0, 1) - iter.prev(1, 1, 0)
          + iter.prev(1, 1, 1);
    };
  };

}
#endif
