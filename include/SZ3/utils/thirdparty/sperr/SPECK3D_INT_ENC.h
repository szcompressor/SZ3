#ifndef SPECK3D_INT_ENC_H
#define SPECK3D_INT_ENC_H

#include "SPECK3D_INT.h"

namespace SZ3 {
namespace SPERR {

//
// Main SPECK3D_INT_ENC class
//
template <typename T>
class SPECK3D_INT_ENC final : public SPECK3D_INT<T> {
 private:
  //
  // Consistant with the base class.
  //
  using uint_type = T;
  using vecui_type = std::vector<uint_type>;

  //
  // Bring members from parent classes to this derived class.
  //
  using SPECK_INT<T>::m_LIP_mask;
  using SPECK_INT<T>::m_dims;
  using SPECK_INT<T>::m_LSP_new;
  using SPECK_INT<T>::m_threshold;
  using SPECK_INT<T>::m_coeff_buf;
  using SPECK_INT<T>::m_bit_buffer;
  using SPECK_INT<T>::m_sign_array;
  using SPECK3D_INT<T>::m_LIS;
  using SPECK3D_INT<T>::m_partition_S_XYZ;
  using SPECK3D_INT<T>::m_code_S;

  void m_process_S(size_t idx1, size_t idx2, size_t& counter, bool output) final;
  void m_process_P(size_t idx, size_t morton, size_t& counter, bool output) final;
  void m_process_P_lite(size_t idx) final;
  void m_additional_initialization() final;

  // Data structures and functions for morton data layout.
  vecui_type m_morton_buf;
  void m_deposit_set(Set3D);
};

}  // namespace SPERR
}  // namespace SZ3
#endif
