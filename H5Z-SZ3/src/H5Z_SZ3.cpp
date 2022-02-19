//
// Created by arham23 on 2/8/22.
//

#include "H5Z_SZ3.h"
#include "H5PLextern.h"
#include <ByteToolkit.h>

//h5repack -f UD=32034,0,7,3,0,1,1,8,8,128 /home/arham23/Software/SZ3/test/testfloat_8_8_128.dat.h5 tf_8_8_128.dat.sz.h5

int MAX_CHUNK_SIZE = 2E32 - 1;

//fitler definition
const H5Z_class2_t H5Z_SZ3[1] = {{
    H5Z_CLASS_T_VERS,              /* H5Z_class_t version */
    (H5Z_filter_t)H5Z_FILTER_SZ3, /* Filter id number */
    1,              /* encoder_present flag (set to true) */
    1,              /* decoder_present flag (set to true) */
    "SZ3 compressor/decompressor for floating-point data.", /* Filter name for debugging */
    NULL,                          /* The "can apply" callback */
    NULL,                          /* The "set local" callback */
    (H5Z_func_t)H5Z_filter_sz3,   /* The actual filter function */
}};

H5PL_type_t H5PLget_plugin_type(void) {return H5PL_TYPE_FILTER;}
const void *H5PLget_plugin_info(void) {return H5Z_SZ3;}

/*FILTER FUNCTIONS*/

static size_t H5Z_filter_sz3(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t* buf_size, void** buf)
{

    printf("Entering filter FN");
    //store dimensions, num_dimensions, and data type
    size_t r1 = 0, r2 = 0, r3 = 0, r4 = 0, r5 = 0;
    int dimSize = 0, dataType = 0;

    //1 if error info included else 0
    int withErrInfo = checkCDValuesWithErrors(cd_nelmts, cd_values);
    int error_mode = 0;
    int cmp_algo = 0;
    int interp_algo = 0;
    float abs_error = 0, rel_error = 0, l2norm_error = 0, psnr = 0;
    if(withErrInfo)
        SZ_cdArrayToMetaDataErr(cd_nelmts, cd_values, &dimSize, &dataType, &cmp_algo, &interp_algo, &r5, &r4, &r3, &r2, &r1, &error_mode, &abs_error, &rel_error, &l2norm_error, &psnr);
    else
        SZ_cdArrayToMetaData(cd_nelmts, cd_values, &dimSize, &dataType, &cmp_algo, &interp_algo,&r5, &r4, &r3, &r2, &r1);


    if(flags & H5Z_FLAG_REVERSE){
        /*decompress data*/

        printf("Decompressing w/ SZ3 ");
        SZ::Config conf;
        size_t nbEle = computeDataLength(r5,r4,r3,r2,r1);

        switch(dataType){
            case 0: //FLOAT
            {
                float *f_decompressedData = new float[nbEle];
                SZ_decompress(conf, (char*) *buf, nbytes, f_decompressedData);
                free(*buf);
                *buf = f_decompressedData;
                *buf_size = nbEle * sizeof(float);
                break;
            }

            case 1: //DOUBLE
            {
                double *d_decompressedData = new double[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, d_decompressedData);
                free(*buf);
                *buf = d_decompressedData;
                *buf_size = nbEle * sizeof(double);
                break;
            }

            case 2: //INT 8
            {
                char *c_decompressedData = new char[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, c_decompressedData);
                free(*buf);
                *buf = c_decompressedData;
                *buf_size = nbEle * sizeof(char);
                break;
            }

            case 3: //UINT 8
            {
                unsigned char *uc_decompressedData = new unsigned char[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, uc_decompressedData);
                free(*buf);
                *buf = uc_decompressedData;
                *buf_size = nbEle * sizeof(unsigned char);
                break;
            }

            case 4: //INT 16
            {
                short *s_decompressedData = new short[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, s_decompressedData);
                free(*buf);
                *buf = s_decompressedData;
                *buf_size = nbEle * sizeof(short);
                break;
            }

            case 5: //UINT 16
            {
                unsigned short *us_decompressedData = new unsigned short[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, us_decompressedData);
                free(*buf);
                *buf = us_decompressedData;
                *buf_size = nbEle * sizeof(unsigned short);
                break;
            }

            case 6: //INT 32
            {
                int *i_decompressedData = new int[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, i_decompressedData);
                free(*buf);
                *buf = i_decompressedData;
                *buf_size = nbEle * sizeof(int);
                break;
            }

            case 7: //UINT 32
            {
                unsigned int *ui_decompressedData = new unsigned int[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, ui_decompressedData);
                free(*buf);
                *buf = ui_decompressedData;
                *buf_size = nbEle * sizeof(unsigned int);
                break;
            }

            case 8: //INT 64
            {
                long *l_decompressedData = new long[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, l_decompressedData);
                free(*buf);
                *buf = l_decompressedData;
                *buf_size = nbEle * sizeof(long);
                break;
            }

            case 9: //UINT 64
            {
                unsigned long *ul_decompressedData = new unsigned long[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, ul_decompressedData);
                free(*buf);
                *buf = ul_decompressedData;
                *buf_size = nbEle * sizeof(unsigned long);
                break;
            }

            default:
            {
                printf("Decompression Error: Unknown Datatype");
                exit(0);
            }
        }

    }
    else{
        /*compress data*/
        printf("Compressing w/ SZ3");
        //based on # dimensions, get relevant dimensions and load config object with them
        if(dimSize <= 0){
            printf("Error: Number of Dimensions is <= 0");
            exit(0);
        }

        SZ::Config conf;
        printf("\nDIMS_CMP:\n");
        printf("r1 %u r2 %u r3 %u r4 %u r5 %u\n", r1,r2,r3,r4,r5);
        if (r2 == 0) {
            conf = SZ::Config(r1);
        } else if (r3 == 0) {
            conf = SZ::Config(r2, r1);
        } else if (r4 == 0) {
            conf = SZ::Config(r3, r2, r1);
        } else if (r5 == 0){
            conf = SZ::Config(r4, r3, r2, r1);
        }
        else {
            conf = SZ::Config(r5, r4, r3, r2, r1);
        }

        if(error_mode < 0 || error_mode > 5){
            printf("Invalid error mode: %i, error mode should be in [0,5]", error_mode);
            exit(0);
        }

        conf.errorBoundMode = error_mode;
        conf.absErrorBound = abs_error;
        conf.relErrorBound = rel_error;
        conf.l2normErrorBound = l2norm_error;
        conf.psnrErrorBound = psnr;

        if(cmp_algo < 0 || cmp_algo > 2){
            printf("Invalid compression algo: %i, should be in [0,2]", cmp_algo);
            exit(0);
        }

        conf.cmprAlgo = cmp_algo;

        if(interp_algo < 0 || interp_algo > 1){
            printf("Invalid interpolation algo: %i, should be either 0 or 1", interp_algo);
            exit(0);
        }

        conf.interpAlgo = interp_algo;

        size_t outSize = 0;
        char* compressedData = NULL;


        switch(dataType){
            case 0: //FLOAT
            {
                printf("Compressing Float Data");
                float* f_data = (float*) *buf;

                if(f_data == NULL){
                    printf("F_DATANULL");
                }

                compressedData = SZ_compress(conf, f_data, outSize);
                printf("Leaving compress float case");
                break;
            }

            case 1: //DOUBLE
            {
                compressedData = SZ_compress(conf, (double*) *buf, outSize);
                break;
            }

            case 2: //INT 8
            {
                compressedData = SZ_compress(conf, (char*) *buf, outSize);
                break;
            }

            case 3: //UINT 8
            {
                compressedData = SZ_compress(conf, (unsigned char*) *buf, outSize);
                break;
            }

            case 4: //INT 16
            {
                compressedData = SZ_compress(conf, (short*) *buf, outSize);
                break;
            }

            case 5: //UINT 16
            {
                compressedData = SZ_compress(conf, (unsigned short*) *buf, outSize);
                break;
            }

            case 6: //INT 32
            {
                compressedData = SZ_compress(conf, (int*) *buf, outSize);
                break;
            }

            case 7: //UINT 32
            {
                compressedData = SZ_compress(conf, (unsigned int*) *buf, outSize);
                break;
            }

            case 8: //INT 64
            {
                compressedData = SZ_compress(conf, (long*) *buf, outSize);
                break;
            }

            case 9: //UINT 64
            {
                compressedData = SZ_compress(conf, (unsigned long*) *buf, outSize);
                break;
            }

            default:
            {
                printf("Compression Error: Unknown Datatype");
                exit(0);
            }
        }

        free(*buf);
        *buf = compressedData;
        *buf_size = outSize;

        printf("Ending Compression Routine");

    }

    return *buf_size;

}



/*HELPER FUNCTIONS*/
//use to convert HDF5 cd_array to SZ params inside filter
void SZ_cdArrayToMetaData(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, int* cmp_algo, int* interp_algo, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1){
    assert(cd_nelmts >= 6);
    unsigned char bytes[8];
    *dimSize = cd_values[0];
    *dataType = cd_values[1];
    *cmp_algo = cd_values[2];
    *interp_algo = cd_values[3];

    switch(*dimSize)
    {
        case 1:
            intToBytes_bigEndian(bytes, cd_values[4]);
            intToBytes_bigEndian(&bytes[4], cd_values[5]);
            if(sizeof(size_t)==4)
                *r1 = (unsigned int)bytesToLong_bigEndian(bytes);
            else
                *r1 = (unsigned long)bytesToLong_bigEndian(bytes);
            *r2 = *r3 = *r4 = *r5 = 0;
            break;
        case 2:
            *r3 = *r4 = *r5 = 0;
            *r2 = cd_values[5];
            *r1 = cd_values[4];
            break;
        case 3:
            *r4 = *r5 = 0;
            *r3 = cd_values[6];
            *r2 = cd_values[5];
            *r1 = cd_values[4];
            break;
        case 4:
            *r5 = 0;
            *r4 = cd_values[7];
            *r3 = cd_values[6];
            *r2 = cd_values[5];
            *r1 = cd_values[4];
            break;
        default:
            *r5 = cd_values[8];
            *r4 = cd_values[7];
            *r3 = cd_values[6];
            *r2 = cd_values[5];
            *r1 = cd_values[4];
    }
}

void SZ_cdArrayToMetaDataErr(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, int* cmp_algo, int* interp_algo, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1,
                             int* error_bound_mode, float* abs_error, float* rel_error, float* l2norm_error, float* psnr)
{
    //get dimension, datatype metadata from cd_values
    SZ_cdArrayToMetaData(cd_nelmts, cd_values, dimSize, dataType, cmp_algo, interp_algo, r5, r4, r3, r2, r1);
    //read in error bound value information
    int dim = *dimSize;
    int k = dim==1?4:dim+2;
    unsigned char b[4];
    int32ToBytes_bigEndian(b, cd_values[k++]);
    *error_bound_mode = bytesToInt32_bigEndian(b);
    int32ToBytes_bigEndian(b, cd_values[k++]);
    *abs_error = bytesToFloat(b);
    int32ToBytes_bigEndian(b, cd_values[k++]);
    *rel_error = bytesToFloat(b);
    int32ToBytes_bigEndian(b, cd_values[k++]);
    *l2norm_error = bytesToFloat(b);
    int32ToBytes_bigEndian(b, cd_values[k++]);
    *psnr = bytesToFloat(b);
}

//use to convert params for SZ into HDF5 cd_array
void SZ_metaDataToCdArray(size_t* cd_nelmts, unsigned int **cd_values, int dataType, int cmp_algo, int interp_algo, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
    *cd_values = (unsigned int*)malloc(sizeof(unsigned int)*7);
    SZ_copymetaDataToCdArray(cd_nelmts, *cd_values, dataType, cmp_algo, interp_algo, r5, r4, r3, r2, r1);
}

void SZ_metaDataErrToCdArray(size_t* cd_nelmts, unsigned int **cd_values, int dataType, int cmp_algo, int interp_algo, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1,
                             int error_bound_mode, float abs_error, float rel_error, float pw_rel_error, float psnr)
{
    *cd_values = (unsigned int*)malloc(sizeof(unsigned int)*12);
    SZ_copymetaDataToCdArray(cd_nelmts, *cd_values, dataType, cmp_algo, interp_algo, r5, r4, r3, r2, r1);
    int k = *cd_nelmts;
    (*cd_values)[k++] = error_bound_mode;
    unsigned char b[4];
    floatToBytes(b, abs_error);
    (*cd_values)[k++] = bytesToInt32_bigEndian(b);
    floatToBytes(b, rel_error);
    (*cd_values)[k++] = bytesToInt32_bigEndian(b);
    floatToBytes(b, pw_rel_error);
    (*cd_values)[k++] = bytesToInt32_bigEndian(b);
    floatToBytes(b, psnr);
    (*cd_values)[k++] = bytesToInt32_bigEndian(b);
    *cd_nelmts = k;
}

void SZ_copymetaDataToCdArray(size_t* cd_nelmts, unsigned int *cd_values, int dataType, int cmp_algo, int interp_algo, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
    unsigned char bytes[8] = {0};
    unsigned long size;
    int dim = computeDimension(r5, r4, r3, r2, r1);
    cd_values[0] = dim;
    cd_values[1] = dataType;	//0: FLOAT ; 1: DOUBLE ; 2,3,4,....: INTEGER....
    cd_values[2] = cmp_algo;
    cd_values[3] = interp_algo;

    switch(dim)
    {
        case 1:
            size = (unsigned long)r1;
            longToBytes_bigEndian(bytes, size);
            cd_values[4] = bytesToInt_bigEndian(bytes);
            cd_values[5] = bytesToInt_bigEndian(&bytes[4]);
            *cd_nelmts = 6;
            break;
        case 2:
            cd_values[4] = (unsigned int) r2;
            cd_values[5] = (unsigned int) r1;
            *cd_nelmts = 6;
            break;
        case 3:
            cd_values[4] = (unsigned int) r3;
            cd_values[5] = (unsigned int) r2;
            cd_values[6] = (unsigned int) r1;
            *cd_nelmts = 7;
            break;
        case 4:
            cd_values[4] = (unsigned int) r4;
            cd_values[5] = (unsigned int) r3;
            cd_values[6] = (unsigned int) r2;
            cd_values[7] = (unsigned int) r1;
            *cd_nelmts = 8;
            break;
        default:
            cd_values[4] = (unsigned int) r5;
            cd_values[5] = (unsigned int) r4;
            cd_values[6] = (unsigned int) r3;
            cd_values[7] = (unsigned int) r2;
            cd_values[8] = (unsigned int) r1;
            *cd_nelmts = 9;
    }
}

int checkCDValuesWithErrors(size_t cd_nelmts, const unsigned int cd_values[])
{
    int result = 0; //0 means no-error-information-in-cd_values; 1 means cd_values contains error information
    int dimSize = cd_values[0];
    //printf("nc_nelmts = %d\n", cd_nelmts);
    switch(dimSize)
    {
        case 1:
            if(cd_nelmts>6)
                result = 1;
            break;
        case 2:
            if(cd_nelmts>6)
                result = 1;
            break;
        case 3:
            if(cd_nelmts>7)
                result = 1;
            break;
        case 4:
            if(cd_nelmts>8)
                result = 1;
            break;
        case 5:
            if(cd_nelmts>9)
                result = 1;
            break;
    }
    return result;
}

size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
    size_t dataLength;
    if(r1==0)
    {
        dataLength = 0;
    }
    else if(r2==0)
    {
        dataLength = r1;
    }
    else if(r3==0)
    {
        dataLength = r1*r2;
    }
    else if(r4==0)
    {
        dataLength = r1*r2*r3;
    }
    else if(r5==0)
    {
        dataLength = r1*r2*r3*r4;
    }
    else
    {
        dataLength = r1*r2*r3*r4*r5;
    }
    return dataLength;
}

int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
    int dimension;
    if(r1==0)
    {
        dimension = 0;
    }
    else if(r2==0)
    {
        dimension = 1;
    }
    else if(r3==0)
    {
        dimension = 2;
    }
    else if(r4==0)
    {
        dimension = 3;
    }
    else if(r5==0)
    {
        dimension = 4;
    }
    else
    {
        dimension = 5;
    }
    return dimension;
}

void init_dims_chunk(int dim, hsize_t dims[5], hsize_t chunk[5], size_t nbEle, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
    switch(dim)
    {
        case 1:
            dims[0] = r1;
            if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
                chunk[0] = r1;
            else
                chunk[0] = 2147483648;//2^31
            break;
        case 2:
            dims[0] = r2;
            dims[1] = r1;
            if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
            {
                chunk[0] = r2;
                chunk[1] = r1;
            }
            else
            {
                printf("Error: size is too big!\n");
                exit(0);
            }
            break;
        case 3:
            dims[0] = r3;
            dims[1] = r2;
            dims[2] = r1;
            if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
            {
                chunk[0] = r3;
                chunk[1] = r2;
                chunk[2] = r1;
            }
            else
            {
                printf("Error: size is too big!\n");
                exit(0);
            }
            break;
        case 4:
            dims[0] = r4;
            dims[1] = r3;
            dims[2] = r2;
            dims[3] = r1;
            if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
            {
                chunk[0] = r4;
                chunk[1] = r3;
                chunk[2] = r2;
                chunk[3] = r1;
            }
            else
            {
                printf("Error: size is too big!\n");
                exit(0);
            }
            break;
        default:
            dims[0] = r5;
            dims[1] = r4;
            dims[2] = r3;
            dims[3] = r2;
            dims[4] = r1;
            if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
            {
                chunk[0] = r5;
                chunk[1] = r4;
                chunk[2] = r3;
                chunk[3] = r2;
                chunk[4] = r1;
            }
            else
            {
                printf("Error: size is too big!\n");
                exit(0);
            }
    }
}