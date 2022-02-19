//
// Created by arham23 on 2/8/22.
//

#ifndef SZ3_H5Z_SZ3_H
#define SZ3_H5Z_SZ3_H

#define H5Z_FILTER_SZ3 32034
#include "SZ3/api/sz.hpp"
#include "hdf5.h"
#include "../test/include/rw.h"

static size_t H5Z_filter_sz3(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t* buf_size, void** buf);

void SZ_cdArrayToMetaData(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, int* cmp_algo, int* interp_algo, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1);

void SZ_cdArrayToMetaDataErr(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, int* cmp_algo, int* interp_algo, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1,
                            int* error_bound_mode, float* abs_error, float* rel_error, float* l2norm_error, float* psnr);

void SZ_metaDataToCdArray(size_t* cd_nelmts, unsigned int **cd_values, int dataType, int cmp_algo, int interp_algo, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

void SZ_metaDataErrToCdArray(size_t* cd_nelmts, unsigned int **cd_values, int dataType, int cmp_algo, int interp_algo, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1,
                             int error_bound_mode, float abs_error, float rel_error, float pw_rel_error, float psnr);

void SZ_copymetaDataToCdArray(size_t* cd_nelmts, unsigned int *cd_values, int dataType, int cmp_algo, int interp_algo, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

int checkCDValuesWithErrors(size_t cd_nelmts, const unsigned int cd_values[]);
size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
void init_dims_chunk(int dim, hsize_t dims[5], hsize_t chunk[5], size_t nbEle, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

#endif //SZ3_H5Z_SZ3_H