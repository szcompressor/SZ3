//
// Created by arham23 on 2/8/22.
//

#ifndef SZ3_H5Z_SZ3_H
#define SZ3_H5Z_SZ3_H

#define H5Z_FILTER_SZ3 32024
#include "SZ3/api/sz.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "hdf5.h"


#define ERROR(FNAME)                                              \
do {                                                              \
    int _errno = errno;                                           \
    fprintf(stderr, #FNAME " failed at line %d, errno=%d (%s)\n", \
        __LINE__, _errno, _errno?strerror(_errno):"ok");          \
    return 1;                                                     \
} while(0)

/*#define H5Z_SZ_PUSH_AND_GOTO(MAJ, MIN, RET, MSG)     \
do                                                    \
{                                                     \
//	H5Epush(H5E_DEFAULT,__FILE__,_funcname_,__LINE__, 
//		H5Z_SZ_ERRCLASS,MAJ,MIN,MSG);                
	retval = RET;                                     \
	goto done;                                        \
} while(0)*/


#define LITTLE_ENDIAN_SYSTEM 0
#define BIG_ENDIAN_SYSTEM 1

int sysEndianType = LITTLE_ENDIAN_SYSTEM;
int dataEndianType = LITTLE_ENDIAN_SYSTEM;

static herr_t H5Z_sz3_set_local(hid_t dcpl_id, hid_t type_id, hid_t chunk_space_id);
static size_t H5Z_filter_sz3(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t* buf_size, void** buf);

void SZ_cdArrayToMetaData(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1);

void SZ_cdArrayToMetaDataErr(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1,
                            int* error_bound_mode, float* abs_error, float* rel_error, float* l2norm_error, float* psnr);

void SZ_metaDataToCdArray(size_t* cd_nelmts, unsigned int **cd_values, int dataType,size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

void SZ_metaDataErrToCdArray(size_t* cd_nelmts, unsigned int **cd_values, int dataType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1,
                             int error_bound_mode, float abs_error, float rel_error, float pw_rel_error, float psnr);

void SZ_copymetaDataToCdArray(size_t* cd_nelmts, unsigned int *cd_values, int dataType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

int checkCDValuesWithErrors(size_t cd_nelmts, const unsigned int cd_values[]);
size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
void init_dims_chunk(int dim, hsize_t dims[5], hsize_t chunk[5], size_t nbEle, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

extern float bytesToFloat(unsigned char* bytes);
extern void floatToBytes(unsigned char *b, float num);

void detectSysEndianType();

#endif //SZ3_H5Z_SZ3_H
