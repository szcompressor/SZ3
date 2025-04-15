#ifndef _SZ_LORENZO_PREDICTOR_HPP
#define _SZ_LORENZO_PREDICTOR_HPP

#include "SZ3/def.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include <cassert>

namespace SZ {

// N-dimension L-layer lorenzo predictor
template <class T, uint N, uint L, class Quantizer>
class LorenzoPredictor : public concepts::PredictorInterface<T, N> {
public:
  using block_iter = typename multi_dimensional_data<T, N>::block_iterator;

  static const uint8_t predictor_id = 0b00000001;

  LorenzoPredictor(Quantizer quantizer) : quantizer(quantizer) {
    this->noise = 0;
  }

  LorenzoPredictor(Quantizer quantizer, double eb) : quantizer(quantizer) {
    this->noise = 0;
    if (L == 1) {
      if (N == 1) {
        this->noise = 0.5 * eb;
      } else if (N == 2) {
        this->noise = 0.81 * eb;
      } else if (N == 3) {
        this->noise = 1.22 * eb;
      } else if (N == 4) {
        this->noise = 1.79 * eb;
      }
    } else if (L == 2) {
      if (N == 1) {
        this->noise = 1.08 * eb;
      } else if (N == 2) {
        this->noise = 2.76 * eb;
      } else if (N == 3) {
        this->noise = 6.8 * eb;
      }
    }
  }

  inline void save(uchar *&c) const {
    c[0] = predictor_id;
    c += sizeof(uint8_t);
    quantizer.save(c);
  }

  inline void load(const uchar *&c, size_t &remaining_length) {
    c += sizeof(uint8_t);
    remaining_length -= sizeof(uint8_t);
    quantizer.load(c, remaining_length);
  }

  void print() {
    std::cout << L << "-Layer " << N << "D Lorenzo predictor, noise = " << noise
              << "\n";
  }

  inline T est_error(const block_iter &block) {
    auto range = block.get_block_range();
    //            auto d = block.mddata;
    auto ds = block.get_dim_strides();

    size_t min_size = std::numeric_limits<size_t>::max();
    for (const auto &r : range) {
      min_size = std::min(min_size, r.second - r.first);
    }

    T err = 0;
    for (size_t i = 2; i < min_size; i++) {
      size_t bmi = min_size - i;
      if constexpr (N == 1) {
        T *c = block.get_data(i);
        T pred = lorenzo_predict(block, c);
        err += fabs(*c - pred) + this->noise;

        c = block.get_data(bmi);
        pred = lorenzo_predict(block, c);
        err += fabs(*c - pred) + this->noise;

      } else if constexpr (N == 2) {
        T *c = block.get_data(i, i);
        T pred = lorenzo_predict(block, c);
        err += fabs(*c - pred) + this->noise;

        c = block.get_data(i, bmi);
        pred = lorenzo_predict(block, c);
        err += fabs(*c - pred) + this->noise;

        c = block.get_data(bmi, i);
        pred = lorenzo_predict(block, c);
        err += fabs(*c - pred) + this->noise;

        c = block.get_data(bmi, bmi);
        pred = lorenzo_predict(block, c);
        err += fabs(*c - pred) + this->noise;
      } else if constexpr (N == 3) {
        T *c = block.get_data(i, i, i);
        T pred = lorenzo_predict(block, c);
        err += fabs(*c - pred) + this->noise;

        c = block.get_data(i, i, bmi);
        pred = lorenzo_predict(block, c);
        err += fabs(*c - pred) + this->noise;

        c = block.get_data(i, bmi, i);
        pred = lorenzo_predict(block, c);
        err += fabs(*c - pred) + this->noise;

        c = block.get_data(i, bmi, bmi);
        pred = lorenzo_predict(block, c);
        err += fabs(*c - pred) + this->noise;
      }
    }
    return err;
  }

  size_t size_est() { return quantizer.size_est(); }

  size_t get_padding() { return 2; }

  // Helper functions for Lorenzo prediction
  ALWAYS_INLINE T prev1(T *c, size_t x) { return c[-x]; }

  ALWAYS_INLINE T prev2(T *c, const std::array<size_t, N> &ds, size_t y,
                        size_t x) {
    return c[-y * ds[0] - x];
  }

  ALWAYS_INLINE T prev3(T *c, const std::array<size_t, N> &ds, size_t z,
                        size_t y, size_t x) {
    return c[-z * ds[1] - y * ds[0] - x];
  }

  ALWAYS_INLINE T prev4(T *c, const std::array<size_t, N> &ds, size_t w,
                        size_t z, size_t y, size_t x) {
    return c[-w * ds[2] - z * ds[1] - y * ds[0] - x];
  }

  // Unified Lorenzo predictor that handles all dimension and layer combinations
  ALWAYS_INLINE T lorenzo_predict(const block_iter &block, T *c) {
    auto ds = block.mddata->get_dim_strides();
    if constexpr (N == 1 && L == 1) {
      return prev1(c, 1);
    } else if constexpr (N == 2 && L == 1) {
      return prev2(c, ds, 0, 1) + prev2(c, ds, 1, 0) - prev2(c, ds, 1, 1);
    } else if constexpr (N == 3 && L == 1) {
      return prev3(c, ds, 0, 0, 1) + prev3(c, ds, 0, 1, 0) +
             prev3(c, ds, 1, 0, 0) - prev3(c, ds, 0, 1, 1) -
             prev3(c, ds, 1, 0, 1) - prev3(c, ds, 1, 1, 0) +
             prev3(c, ds, 1, 1, 1);
    } else if constexpr (N == 4 && L == 1) {
      return prev4(c, ds, 0, 0, 0, 1) + prev4(c, ds, 0, 0, 1, 0) -
             prev4(c, ds, 0, 0, 1, 1) + prev4(c, ds, 0, 1, 0, 0) -
             prev4(c, ds, 0, 1, 0, 1) - prev4(c, ds, 0, 1, 1, 0) +
             prev4(c, ds, 0, 1, 1, 1) + prev4(c, ds, 1, 0, 0, 0) -
             prev4(c, ds, 1, 0, 0, 1) - prev4(c, ds, 1, 0, 1, 0) +
             prev4(c, ds, 1, 0, 1, 1) - prev4(c, ds, 1, 1, 0, 0) +
             prev4(c, ds, 1, 1, 0, 1) + prev4(c, ds, 1, 1, 1, 0) -
             prev4(c, ds, 1, 1, 1, 1);
    } else if constexpr (N == 1 && L == 2) {
      return 2 * prev1(c, 1) - prev1(c, 2);
    } else if constexpr (N == 2 && L == 2) {
      return 2 * prev2(c, ds, 0, 1) - prev2(c, ds, 0, 2) +
             2 * prev2(c, ds, 1, 0) - 4 * prev2(c, ds, 1, 1) +
             2 * prev2(c, ds, 1, 2) - prev2(c, ds, 2, 0) +
             2 * prev2(c, ds, 2, 1) - prev2(c, ds, 2, 2);
    } else if constexpr (N == 3 && L == 2) {
      return 2 * prev3(c, ds, 0, 0, 1) - prev3(c, ds, 0, 0, 2) +
             2 * prev3(c, ds, 0, 1, 0) - 4 * prev3(c, ds, 0, 1, 1) +
             2 * prev3(c, ds, 0, 1, 2) - prev3(c, ds, 0, 2, 0) +
             2 * prev3(c, ds, 0, 2, 1) - prev3(c, ds, 0, 2, 2) +
             2 * prev3(c, ds, 1, 0, 0) - 4 * prev3(c, ds, 1, 0, 1) +
             2 * prev3(c, ds, 1, 0, 2) - 4 * prev3(c, ds, 1, 1, 0) +
             8 * prev3(c, ds, 1, 1, 1) - 4 * prev3(c, ds, 1, 1, 2) +
             2 * prev3(c, ds, 1, 2, 0) - 4 * prev3(c, ds, 1, 2, 1) +
             2 * prev3(c, ds, 1, 2, 2) - prev3(c, ds, 2, 0, 0) +
             2 * prev3(c, ds, 2, 0, 1) - prev3(c, ds, 2, 0, 2) +
             2 * prev3(c, ds, 2, 1, 0) - 4 * prev3(c, ds, 2, 1, 1) +
             2 * prev3(c, ds, 2, 1, 2) - prev3(c, ds, 2, 2, 0) +
             2 * prev3(c, ds, 2, 2, 1) - prev3(c, ds, 2, 2, 2);
    }
    // Handle unsupported cases
    else {
      static_assert(N <= 4 && L <= 2,
                    "Unsupported dimension or layer configuration");
      return T(0);
    }
  }

  void compress(const block_iter &block, std::vector<int> &quant_inds) {
    block.iterate_block([&](T *c) {
      T pred = lorenzo_predict(block, c);
      quant_inds.push_back(quantizer.quantize_and_overwrite(*c, pred));
    });
  }

  void decompress(const block_iter &block, int *&quant_inds_pos) {
    block.iterate_block([&](T *c) {
      T pred = lorenzo_predict(block, c);
      *c = quantizer.recover(pred, *(quant_inds_pos++));
    });
  }

  void clear() {}

protected:
  T noise = 0;

private:
  Quantizer quantizer;
};
} // namespace SZ
#endif
