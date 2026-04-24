#ifndef SZ3_THIRDPARTY_MGARD_HEADER_ONLY_HPP
#define SZ3_THIRDPARTY_MGARD_HEADER_ONLY_HPP

/**
 * @file MGARDHeaderOnly.hpp
 * @brief Single include for the bundled MGARD multigrid transform.
 *
 * Bundled (and trimmed) from https://github.com/lxAltria/MGARDx — provides only
 * the forward/inverse multigrid transform; quantization, entropy coding, and
 * lossless backends are layered on top by SZ3.
 */

#include "SZ3/utils/thirdparty/mgard/decompose.hpp"
#include "SZ3/utils/thirdparty/mgard/recompose.hpp"

#endif
