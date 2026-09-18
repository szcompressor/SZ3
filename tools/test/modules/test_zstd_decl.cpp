// Lossless_zstd.hpp declares the four Zstd Simple-API functions SZ3 calls instead of including
// <zstd.h>. Nothing else checks that those declarations still match the library SZ3 links: a
// divergence compiles cleanly on both sides and goes wrong at the call boundary. This translation
// unit includes both so a changed signature is a redeclaration conflict, and
// zstd_decl_header_first.cpp includes them the other way round.
// clang-format off
// The order is the test: clang-format would sort <zstd.h> first, which is what
// zstd_decl_header_first.cpp already covers, leaving this file testing nothing.
#include "SZ3/lossless/Lossless_zstd.hpp"  // SZ3's own declarations FIRST, no zstd.h
#include <zstd.h>                          // then the real header; the two must agree
// clang-format on

#include <cstdint>
#include <numeric>
#include <random>
#include <vector>

#include "gtest/gtest.h"

// defined in zstd_decl_header_first.cpp, which includes the same two headers in the other order
unsigned sz3_zstd_decl_header_first_version();

namespace {

// The declarations bind to the library SZ3 actually linked, not to some other Zstd.
TEST(ZstdDecl, DeclaredBoundMatchesHeaderMacro) {
    // ZSTD_COMPRESSBOUND is zstd.h's macro form of the function we declare ourselves. If our
    // declaration named a different function, or the library disagreed with the header it shipped
    // with, these would part company.
    for (size_t n : {size_t(0), size_t(1), size_t(17), size_t(128), size_t(1u << 10), size_t(1u << 16),
                     size_t(1u << 20), size_t((1u << 20) + 7)}) {
        EXPECT_EQ(ZSTD_compressBound(n), static_cast<size_t>(ZSTD_COMPRESSBOUND(n))) << "srcSize=" << n;
    }
}

TEST(ZstdDecl, BothIncludeOrdersSeeTheSameLibrary) {
    EXPECT_EQ(sz3_zstd_decl_header_first_version(), static_cast<unsigned>(ZSTD_versionNumber()));
}

// A round trip through the wrapper, sized by the wrapper's own bound. This is what the three
// api/impl headers do, so if compress_bound ever under-reported, compress() would throw here.
TEST(ZstdDecl, RoundTripThroughDeclaredFunctions) {
    std::mt19937 rng(20260918);
    std::vector<SZ3::uchar> src(1u << 18);
    // Half structured, half random: compressible enough to exercise the encoder, incompressible
    // enough that the bound is not trivially slack.
    std::iota(src.begin(), src.begin() + src.size() / 2, SZ3::uchar(0));
    for (size_t i = src.size() / 2; i < src.size(); i++) {
        src[i] = static_cast<SZ3::uchar>(rng() & 0xff);
    }

    SZ3::Lossless_zstd zstd;
    std::vector<SZ3::uchar> dst(SZ3::Lossless_zstd::compress_bound(src.size()) + sizeof(size_t));
    size_t cmp_size = zstd.compress(src.data(), src.size(), dst.data(), dst.size());
    ASSERT_GT(cmp_size, sizeof(size_t));
    ASSERT_LE(cmp_size, dst.size());

    std::vector<SZ3::uchar> out(src.size());
    SZ3::uchar *out_p = out.data();
    size_t out_size = zstd.decompress(dst.data(), cmp_size, out_p, out.size());
    ASSERT_EQ(out_size, src.size());
    EXPECT_EQ(out, src);
}

}  // namespace
