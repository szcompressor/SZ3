#ifndef SZ3_SPERR_TYPES_HPP
#define SZ3_SPERR_TYPES_HPP

#include <cstddef>
#include <cstdint>
#include <variant>
#include <vector>

#include "SZ3/utils/thirdparty/sperr/SPERRHeaderOnly.hpp"

namespace SZ3 {

/// Quality control mode interpreted from SZ3 `errorBoundMode`.
enum class SPERRMode { PSNR_MODE, PWE_MODE };

/// 3D shape in SPERR ordering `[x, y, z]`.
struct SPERR3DLayout {
    size_t dimx = 0;
    size_t dimy = 0;
    size_t dimz = 0;
};

using SPERRUIntVec = std::variant<std::vector<uint8_t>, std::vector<uint16_t>, std::vector<uint32_t>,
                                  std::vector<uint64_t>>;

/// Intermediate state exchanged between decomposition and encoding modules.
struct SPERRFrame {
    SZ3::SPERR::dims_type dims = {0, 0, 0};
    SZ3::SPERR::condi_type conditioner_header{};
    bool constant_field = false;

    SPERRMode mode = SPERRMode::PWE_MODE;
    double quality = 0.0;
    double conditioned_range = 0.0;
    double q = 0.0;

    SZ3::SPERR::UINTType uint_flag = SZ3::SPERR::UINTType::UINT8;
    SPERRUIntVec coeffs_ui;
    SZ3::SPERR::Bitmask sign_array;
    std::vector<double> wavelet_coeffs;
    std::vector<double> conditioned_values;
    std::vector<SZ3::SPERR::Outlier> outliers;

    SPERRFrame() : coeffs_ui(std::vector<uint8_t>{}) {}
};

}  // namespace SZ3

#endif
