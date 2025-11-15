/**
 * @file stdafx.hpp
 * @brief Precompiled header for SZ3.
 *
 * This file includes large and frequently used headers, such as Eigen,
 * to be pre-compiled. Pre-compiling these headers significantly speeds up
 * the overall build process by avoiding repeated parsing of these large files.
 */
#ifndef SZ3_STDAFX_HPP
#define SZ3_STDAFX_HPP

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <unsupported/Eigen/CXX11/Tensor>

#endif //SZ3_STDAFX_HPP
