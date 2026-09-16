/**
 * @file sz_dev.hpp
 * @ingroup API
 * @brief Developer header that contains all modules in SZ3.
 *
 * Most users should use sz.hpp, which provides the standard SZ3
 * compress/decompress API along with popular modules and utilities.
 *
 * Include this header only if you need to use modules listed below.
 *
 * @code
 * #include "SZ3/api/sz_dev.hpp"
 * @endcode
 */

#ifndef SZ3_SZ_DEV_HPP
#define SZ3_SZ_DEV_HPP

// Core API
#include "SZ3/api/sz.hpp"
#include "SZ3/def.hpp"

// --- Decompositions ---
#include "SZ3/decomposition/MultiLevelDecomposition.hpp"
#include "SZ3/decomposition/NoPredictionDecomposition.hpp"
#include "SZ3/decomposition/PaSTRIDecomposition.hpp"
#include "SZ3/decomposition/SPERRDecomposition.hpp"
#include "SZ3/decomposition/TimeSeriesDecomposition.hpp"
#include "SZ3/decomposition/ZFPDecomposition.hpp"

// --- Compressors ---
#include "SZ3/compressor/specialized/SZExaaltCompressor.hpp"
#include "SZ3/compressor/specialized/SZTruncateCompressor.hpp"

// --- Predictors ---

// --- Quantizers ---
#include "SZ3/quantizer/ClusterQuantizer.hpp"
#include "SZ3/quantizer/GranularBitRoundQuantizer.hpp"
#include "SZ3/quantizer/LevelQuantizer.hpp"
#include "SZ3/quantizer/LogDomainQuantizer.hpp"
#include "SZ3/quantizer/OutlierQuantizer.hpp"

// --- Encoders ---
#include "SZ3/encoder/ArithmeticEncoder.hpp"
#include "SZ3/encoder/BitshuffleEncoder.hpp"
#include "SZ3/encoder/BypassEncoder.hpp"
#include "SZ3/encoder/HuffmanEncoderV2.hpp"
#include "SZ3/encoder/RunlengthEncoder.hpp"
#include "SZ3/encoder/XtcBasedEncoder.hpp"
#include "SZ3/encoder/ZFPEncoder.hpp"

// --- Lossless ---
#include "SZ3/lossless/Lossless_bypass.hpp"

// --- Preprocessors ---
#include "SZ3/preprocessor/MGARDTransform.hpp"
#include "SZ3/preprocessor/PreFilter.hpp"
#include "SZ3/preprocessor/PreProcessor.hpp"
#include "SZ3/preprocessor/SPERRTransform.hpp"
#include "SZ3/preprocessor/Transpose.hpp"
#include "SZ3/preprocessor/Wavelet.hpp"

// --- Utilities ---
#include "SZ3/utils/Collections.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/KmeansUtil.hpp"
#include "SZ3/utils/MultiLevelErrorBound.hpp"
#include "SZ3/utils/MultiLevelQuantization.hpp"
#include "SZ3/utils/QuantOptimization.hpp"

#endif  // SZ3_SZ_DEV_HPP
