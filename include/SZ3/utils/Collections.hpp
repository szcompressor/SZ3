#ifndef SZ3_COLLECTIONS_HPP
#define SZ3_COLLECTIONS_HPP

// INTPTR_MAX and INT64_MAX come from here. Without it the preprocessor reads both as 0, the test
// below is 0 == 0, and a 32-bit build silently takes the 64-bit branch.
#include <cstdint>

// SZ3_USE_SKA_HASH is baked in at configure time, not passed as an INTERFACE define: it decides
// which type SZ3::unordered_map names, and two translation units that disagree give every class
// template holding one -- HuffmanEncoder among them -- two different layouts in one program.
#include "SZ3/options.hpp"

// One condition for the include and the alias together. Split, the else branch named
// std::unordered_map without ever including <unordered_map>.
#if SZ3_USE_SKA_HASH && INTPTR_MAX == INT64_MAX

// Vendored, and kept byte-identical with upstream skarupke: its parameters are named after the
// members they initialise, and it offsets a null pointer to form its end sentinel. Silence both
// here rather than patching the copy, so a re-vendor has nothing to re-apply.
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif
#if defined(__clang__)
#pragma GCC diagnostic ignored "-Wnull-pointer-arithmetic"
#endif
#include "SZ3/utils/ska_hash/unordered_map.hpp"
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace SZ3 {
using ska::unordered_map;
}

#else

#include <unordered_map>

namespace SZ3 {
using std::unordered_map;
}

#endif

#endif  // SZ3_COLLECTIONS_HPP
