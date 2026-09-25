#ifndef SZ3_H5Z_SZ3_H
#define SZ3_H5Z_SZ3_H

/* The SZ3 HDF5 filter. Everything outside the __cplusplus block is plain C, so C programs can
 * include this header too; the functions that take SZ3::Config are C++ only. */

#include "hdf5.h"
#include "H5PLextern.h"

#define H5Z_FILTER_SZ3 32024

#ifdef _WIN32
#ifdef hdf5sz3_EXPORTS
#define HDF5SZ3_EXPORT __declspec(dllexport)
#else
#define HDF5SZ3_EXPORT __declspec(dllimport)
#endif
#else
#define HDF5SZ3_EXPORT
#endif

/* SZ3::EB and SZ3::ALGO, repeated because C cannot see them. */
#define H5Z_SZ3_EB_ABS 0
#define H5Z_SZ3_EB_REL 1
#define H5Z_SZ3_EB_PSNR 2
#define H5Z_SZ3_EB_L2NORM 3
#define H5Z_SZ3_EB_ABS_AND_REL 4
#define H5Z_SZ3_EB_ABS_OR_REL 5

#define H5Z_SZ3_ALGO_LORENZO_REG 0
#define H5Z_SZ3_ALGO_INTERP_LORENZO 1
#define H5Z_SZ3_ALGO_INTERP 2
#define H5Z_SZ3_ALGO_NOPRED 3
#define H5Z_SZ3_ALGO_LOSSLESS 4
#define H5Z_SZ3_ALGO_BIOMD 5
#define H5Z_SZ3_ALGO_BIOMDXTC 6

#ifdef __cplusplus
extern "C" {
#endif

/* Put the SZ3 filter on a dataset-creation property list. Only the bounds the mode uses are read;
 * pass 0 for the others. Returns 1 on success, -1 on failure with the reason on the HDF5 error
 * stack. */
HDF5SZ3_EXPORT herr_t H5Pset_sz3(hid_t plist, int algo, int eb_mode, double abs_bound, double rel_bound,
                                 double psnr_bound, double l2norm_bound);

#ifdef __cplusplus
}

#include "SZ3/api/sz.hpp"

extern "C" {

/* Put the SZ3 filter on a dataset-creation property list with a full SZ3::Config.
 * Returns 1 on success, 0 on failure. */
HDF5SZ3_EXPORT herr_t set_SZ3_conf_to_H5(const hid_t propertyList, SZ3::Config &conf);

/* Load the SZ3::Config this property list carries. Returns 1 if a Config was loaded, 0 if the list
 * carries no SZ3 filter, -1 if reading it failed. conf is left alone unless 1 is returned. */
HDF5SZ3_EXPORT herr_t get_SZ3_conf_from_H5(const hid_t propertyList, SZ3::Config &conf);
}
#endif

#endif /* SZ3_H5Z_SZ3_H */
