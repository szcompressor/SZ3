// Predictor modules under include/SZ3/predictor/ are tightly coupled to the
// BlockwiseDecomposition driver and to `block_data<T,N>::block_iterator`
// fixtures (see include/SZ3/utils/BlockwiseIterator.hpp). Building a real
// Block iterator in isolation requires reproducing Decomposition's plumbing
// (data_with_padding, block_iter ctor args, foreach), which is brittle.
//
// The behavioral roundtrip for `LorenzoPredictor` is already exercised end-to-
// end by `SZ3_DecompositionTest.BlockwiseDecomposition2D_Lorenzo` in
// test_decomposition.cpp. For `RegressionPredictor` and `ComposedPredictor`,
// we keep header-compile sanity tests here, mirroring `test_sz_dev.cpp` --
// this catches signature drift (virtual overrides, constructor changes) at CI
// time without trying to drive the predictor outside its block harness.

#include <memory>
#include <stdexcept>
#include <vector>

#include "SZ3/predictor/ComposedPredictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/predictor/RegressionPredictor.hpp"
#include "gtest/gtest.h"

TEST(SZ3_PredictorTest, LorenzoPredictorConstructs) {
    SZ3::LorenzoPredictor<float, 1, 1> p1d_l1(1e-3);
    SZ3::LorenzoPredictor<float, 2, 1> p2d_l1(1e-3);
    SZ3::LorenzoPredictor<float, 3, 1> p3d_l1(1e-3);
    SZ3::LorenzoPredictor<float, 3, 2> p3d_l2(1e-3);
    EXPECT_EQ(p1d_l1.get_padding(), 2u);
    EXPECT_EQ(p3d_l2.get_padding(), 2u);
}

TEST(SZ3_PredictorTest, RegressionPredictorConstructs) {
    SZ3::RegressionPredictor<float, 2> p2d(8, 1e-3);
    SZ3::RegressionPredictor<float, 3> p3d(8, 1e-3);
    SZ3::RegressionPredictor<float, 3> default_ctor;
    SZ3::uchar buf[256] = {0};
    SZ3::uchar* sp = buf;
    p2d.save(sp);
    EXPECT_GT(static_cast<size_t>(sp - buf), 0u);  // header (selection size) saved
}

TEST(SZ3_PredictorTest, ComposedPredictorConstructs) {
    using Pred = SZ3::concepts::PredictorInterface<float, 2>;
    std::vector<std::shared_ptr<Pred>> preds;
    preds.push_back(std::make_shared<SZ3::LorenzoPredictor<float, 2, 1>>(1e-3));
    preds.push_back(std::make_shared<SZ3::RegressionPredictor<float, 2>>(8, 1e-3));
    SZ3::ComposedPredictor<float, 2> composed(preds);
    EXPECT_EQ(composed.get_padding(), 2u);
}

TEST(SZ3_PredictorTest, ComposedPredictorRejectsEmptyList) {
    using Pred = SZ3::concepts::PredictorInterface<float, 2>;
    std::vector<std::shared_ptr<Pred>> preds;
    auto make = [&]() { return SZ3::ComposedPredictor<float, 2>(preds); };
    EXPECT_THROW(make(), std::invalid_argument);
}
