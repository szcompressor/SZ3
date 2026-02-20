#ifndef SZ3_SPERR_HEADER_ONLY_HPP
#define SZ3_SPERR_HEADER_ONLY_HPP

// Avoid collision with sz3c's global macro `PSNR`.
#ifdef PSNR
#pragma push_macro("PSNR")
#undef PSNR
#define SZ3_SPERR_RESTORE_PSNR_MACRO
#endif

#include "Bitmask.h"
#include "Bitstream.h"

#include "CDF97.h"
#include "Conditioner.h"
#include "Outlier_Coder.h"
#include "SPECK1D_INT.h"
#include "SPECK1D_INT_DEC.h"
#include "SPECK1D_INT_ENC.h"
#include "SPECK3D_INT.h"
#include "SPECK3D_INT_DEC.h"
#include "SPECK3D_INT_ENC.h"
#include "SPECK_INT.h"
#include "sperr_helper.h"

#include "sperr_core_impl.hpp"

#ifdef SZ3_SPERR_RESTORE_PSNR_MACRO
#pragma pop_macro("PSNR")
#undef SZ3_SPERR_RESTORE_PSNR_MACRO
#endif

#endif
