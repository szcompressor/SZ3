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
  using block_iter = typename block_data<T, N>::block_iterator;

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
        T *c = block.get_block_data(i);
        T pred = predict(block, c, {i});
        err += fabs(*c - pred) + this->noise;

        c = block.get_block_data(bmi);
        pred = predict(block, c, {bmi});
        err += fabs(*c - pred) + this->noise;

      } else if constexpr (N == 2) {
        T *c = block.get_block_data(i, i);
        T pred = predict(block, c, {i, i});
        err += fabs(*c - pred) + this->noise;

        c = block.get_block_data(i, bmi);
        pred = predict(block, c, {i, bmi});
        err += fabs(*c - pred) + this->noise;

        c = block.get_block_data(bmi, i);
        pred = predict(block, c, {bmi, i});
        err += fabs(*c - pred) + this->noise;

        c = block.get_block_data(bmi, bmi);
        pred = predict(block, c, {bmi, bmi});
        err += fabs(*c - pred) + this->noise;
      } else if constexpr (N == 3) {
        T *c = block.get_block_data(i, i, i);
        T pred = predict(block, c, {i, i, i});
        err += fabs(*c - pred) + this->noise;

        c = block.get_block_data(i, i, bmi);
        pred = predict(block, c, {i, i, bmi});
        err += fabs(*c - pred) + this->noise;

        c = block.get_block_data(i, bmi, i);
        pred = predict(block, c, {i, bmi, i});
        err += fabs(*c - pred) + this->noise;

        c = block.get_block_data(i, bmi, bmi);
        pred = predict(block, c, {i, bmi, bmi});
        err += fabs(*c - pred) + this->noise;
      }
    }
    return err;
  }

  size_t size_est() { return quantizer.size_est(); }

  size_t get_padding() { return 2; }

  // Helper functions for Lorenzo prediction
  ALWAYS_INLINE T prev1(T *d, size_t i) { return d[-i]; }

  ALWAYS_INLINE T prev2(T *d, const std::array<size_t, N> &ds, size_t j,
                        size_t i) {
    return d[-j * ds[0] - i];
  }

  ALWAYS_INLINE T prev3(T *d, const std::array<size_t, N> &ds, size_t k,
                        size_t j, size_t i) {
    return d[-k * ds[1] - j * ds[0] - i];
  }

  ALWAYS_INLINE T prev4(T *d, const std::array<size_t, N> &ds, size_t t,
                        size_t k, size_t j, size_t i) {
    return d[-t * ds[2] - k * ds[1] - j * ds[0] - i];
  }

  ALWAYS_INLINE T predict(const block_iter &block, T *d, const std::array<size_t, N>& index) {
    auto ds = block.get_dim_strides();
    if constexpr (N == 1 && L == 1) {
      return prev1(d, 1);
    } else if constexpr (N == 2 && L == 1) {
      return prev2(d, ds, 0, 1) + prev2(d, ds, 1, 0) - prev2(d, ds, 1, 1);
    } else if constexpr (N == 3 && L == 1) {
      return prev3(d, ds, 0, 0, 1) + prev3(d, ds, 0, 1, 0) +
             prev3(d, ds, 1, 0, 0) - prev3(d, ds, 0, 1, 1) -
             prev3(d, ds, 1, 0, 1) - prev3(d, ds, 1, 1, 0) +
             prev3(d, ds, 1, 1, 1);
    } else if constexpr (N == 4 && L == 1) {
      return prev4(d, ds, 0, 0, 0, 1) + prev4(d, ds, 0, 0, 1, 0) -
             prev4(d, ds, 0, 0, 1, 1) + prev4(d, ds, 0, 1, 0, 0) -
             prev4(d, ds, 0, 1, 0, 1) - prev4(d, ds, 0, 1, 1, 0) +
             prev4(d, ds, 0, 1, 1, 1) + prev4(d, ds, 1, 0, 0, 0) -
             prev4(d, ds, 1, 0, 0, 1) - prev4(d, ds, 1, 0, 1, 0) +
             prev4(d, ds, 1, 0, 1, 1) - prev4(d, ds, 1, 1, 0, 0) +
             prev4(d, ds, 1, 1, 0, 1) + prev4(d, ds, 1, 1, 1, 0) -
             prev4(d, ds, 1, 1, 1, 1);
    } else if constexpr (N == 1 && L == 2) {
      return 2 * prev1(d, 1) - prev1(d, 2);
    } else if constexpr (N == 2 && L == 2) {
      return 2 * prev2(d, ds, 0, 1) - prev2(d, ds, 0, 2) +
             2 * prev2(d, ds, 1, 0) - 4 * prev2(d, ds, 1, 1) +
             2 * prev2(d, ds, 1, 2) - prev2(d, ds, 2, 0) +
             2 * prev2(d, ds, 2, 1) - prev2(d, ds, 2, 2);
    } else if constexpr (N == 3 && L == 2) {
      return 2 * prev3(d, ds, 0, 0, 1) - prev3(d, ds, 0, 0, 2) +
             2 * prev3(d, ds, 0, 1, 0) - 4 * prev3(d, ds, 0, 1, 1) +
             2 * prev3(d, ds, 0, 1, 2) - prev3(d, ds, 0, 2, 0) +
             2 * prev3(d, ds, 0, 2, 1) - prev3(d, ds, 0, 2, 2) +
             2 * prev3(d, ds, 1, 0, 0) - 4 * prev3(d, ds, 1, 0, 1) +
             2 * prev3(d, ds, 1, 0, 2) - 4 * prev3(d, ds, 1, 1, 0) +
             8 * prev3(d, ds, 1, 1, 1) - 4 * prev3(d, ds, 1, 1, 2) +
             2 * prev3(d, ds, 1, 2, 0) - 4 * prev3(d, ds, 1, 2, 1) +
             2 * prev3(d, ds, 1, 2, 2) - prev3(d, ds, 2, 0, 0) +
             2 * prev3(d, ds, 2, 0, 1) - prev3(d, ds, 2, 0, 2) +
             2 * prev3(d, ds, 2, 1, 0) - 4 * prev3(d, ds, 2, 1, 1) +
             2 * prev3(d, ds, 2, 1, 2) - prev3(d, ds, 2, 2, 0) +
             2 * prev3(d, ds, 2, 2, 1) - prev3(d, ds, 2, 2, 2);
    }
    // Handle unsupported cases
    else {
      static_assert(N <= 4 && L <= 2,
                    "Unsupported dimension or layer configuration");
      return T(0);
    }
  }



  void compress(const block_iter &block, std::vector<int> &quant_inds) {
    block_iter::foreach (block, [&](T *c, const std::array<size_t, N>& index) {
      T pred = predict(block, c, index);
      quant_inds.push_back(quantizer.quantize_and_overwrite(*c, pred));
    });
  }

  void decompress(const block_iter &block, int *&quant_inds_pos) {
    block_iter::foreach (block, [&](T *c, const std::array<size_t, N>& index) {
      T pred = predict(block, c, index);
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
