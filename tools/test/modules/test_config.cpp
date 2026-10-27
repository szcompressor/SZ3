#include <vector>

#include "SZ3/utils/Config.hpp"
#include "gtest/gtest.h"

namespace {

// A config as a later version would write it: three more bytes at the end, counted in its size byte.
std::vector<SZ3::uchar> with_appended_fields(const SZ3::Config &conf) {
    std::vector<SZ3::uchar> buf(conf.size_est() + 16);
    SZ3::uchar *c = buf.data();
    const size_t size = conf.save(c);
    buf.resize(size);
    buf.insert(buf.end(), {7, 7, 7});
    buf[0] = static_cast<SZ3::uchar>(size + 3);
    return buf;
}

}  // namespace

TEST(SZ3_Config, LoadSkipsFieldsAddedAtTheEnd) {
    SZ3::Config conf(10, 20, 30);
    conf.absErrorBound = 1e-3;
    conf.openmp = true;
    auto buf = with_appended_fields(conf);
    buf.push_back(0xAB);  // what follows the config, such as the next chunk's config in an OpenMP stream

    SZ3::Config loaded;
    const SZ3::uchar *c = buf.data();
    size_t remaining = buf.size();
    loaded.load(c, remaining);
    EXPECT_EQ(*c, 0xAB);
    EXPECT_EQ(remaining, 1u);
    EXPECT_EQ(loaded.num, conf.num);
    EXPECT_EQ(loaded.absErrorBound, conf.absErrorBound);
    EXPECT_TRUE(loaded.openmp);
}

TEST(SZ3_Config, LoadRefusesASizeShorterThanItsFields) {
    SZ3::Config conf(10, 20, 30);
    std::vector<SZ3::uchar> buf(conf.size_est() + 16);
    SZ3::uchar *c = buf.data();
    buf.resize(conf.save(c));
    buf[0] = 4;

    SZ3::Config loaded;
    const SZ3::uchar *p = buf.data();
    size_t remaining = buf.size();
    EXPECT_THROW(loaded.load(p, remaining), std::out_of_range);
}

TEST(SZ3_Config, ConstructorRejectsZeroDimensions) {
    EXPECT_THROW((SZ3::Config(0, 17, 3)), std::invalid_argument);
    EXPECT_THROW((SZ3::Config(17, 0, 3)), std::invalid_argument);
    EXPECT_THROW((SZ3::Config(17, 3, 0)), std::invalid_argument);
    EXPECT_THROW((SZ3::Config(0)), std::invalid_argument);
}

TEST(SZ3_Config, SetDimsRejectsZeroWithoutChangingTheConfig) {
    SZ3::Config conf(2, 17, 3);
    const auto original_dims = conf.dims;
    const auto original_num = conf.num;
    const auto original_n = conf.N;
    for (const std::vector<size_t> dims : {std::vector<size_t>{0, 17, 3}, {17, 0, 3}, {17, 3, 0}, {0}}) {
        EXPECT_THROW(conf.setDims(dims.begin(), dims.end()), std::invalid_argument);
        EXPECT_EQ(conf.dims, original_dims);
        EXPECT_EQ(conf.num, original_num);
        EXPECT_EQ(conf.N, original_n);
    }
}

TEST(SZ3_Config, SingletonDimensionsStillNormalize) {
    SZ3::Config atom(1, 1, 3);
    EXPECT_EQ(atom.dims, (std::vector<size_t>{3}));
    EXPECT_EQ(atom.num, 3u);
    SZ3::Config scalar(1, 1, 1);
    EXPECT_EQ(scalar.dims, (std::vector<size_t>{1}));
    EXPECT_EQ(scalar.num, 1u);
}
