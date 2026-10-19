#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <vector>

#include "SZ3/utils/FileUtil.hpp"
#include "gtest/gtest.h"

template <class T>
void textRoundTrip(const char *path) {
    std::vector<T> data = {T(10000.125), T(-10000.125), T(1.25), T(0), std::nextafter(T(1), T(2)),
                           std::numeric_limits<T>::min(), std::numeric_limits<T>::denorm_min(),
                           std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max()};
    SZ3::writeTextFile(path, data.data(), data.size());

    std::ifstream in(path);
    ASSERT_TRUE(in.is_open());
    for (T expected : data) {
        T actual = 0;
        ASSERT_TRUE(bool(in >> actual));
        EXPECT_EQ(actual, expected);
    }
    in >> std::ws;
    EXPECT_TRUE(in.eof());
}

TEST(FileUtilTest, TextRoundTrip) {
    const char *path = "file_util_round_trip.txt";
    textRoundTrip<float>(path);
    textRoundTrip<double>(path);
    std::remove(path);
}
