// The other include order for test_zstd_decl.cpp: the real <zstd.h> first, then SZ3's own
// declarations. A translation unit can only have one order, so the pair needs two files.
//
// The orders are not interchangeable. A declaration that adds an attribute the other lacks is
// diagnosed in one direction and accepted in the other, and a consumer is as likely to include
// zstd.h before SZ3 as after. Both GCC and Clang accept the pair either way today, including
// under -fvisibility=hidden; this file is what keeps that true.
//
// It lives beside tools/test/CMakeLists.txt rather than in modules/, because the glob there makes
// one test executable per .cpp and this file is a second translation unit of an existing one, not
// a test of its own.
// clang-format off
// As in test_zstd_decl.cpp, the order here is the test; this file is the opposite order of that
// one. Do not remove these markers or reorder the two includes below.
#include <zstd.h>                          // the real header FIRST
#include "SZ3/lossless/Lossless_zstd.hpp"  // then SZ3's own declarations
// clang-format on

unsigned sz3_zstd_decl_header_first_version() { return static_cast<unsigned>(ZSTD_versionNumber()); }
