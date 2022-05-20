//
// Created by arham23 on 2/8/22.
//

#include "H5Z_SZ3.hpp"
#include <fstream>
#include "H5PLextern.h"

#define CONFIG_PATH "sz3.config"

using namespace SZ;

//h5repack -f UD=32024,0 /home/arham23/Software/SZ3/test/testfloat_8_8_128.dat.h5 tf_8_8_128.dat.sz.h5

//load from "sz3.config" in local directory if 1 else use default values or cd values

int loadConfigFile = 0;
int freshCdValues = 0;
int MAX_CHUNK_SIZE = 2E32 - 1;

//filter definition
const H5Z_class2_t H5Z_SZ3[1] = {{
    H5Z_CLASS_T_VERS,              /* H5Z_class_t version */
    (H5Z_filter_t)H5Z_FILTER_SZ3, /* Filter id number */
    1,              /* encoder_present flag (set to true) */
    1,              /* decoder_present flag (set to true) */
    "SZ3 compressor/decompressor for floating-point data.", /* Filter name for debugging */
    NULL,                          /* The "can apply" callback */
    H5Z_sz3_set_local,                          /* The "set local" callback */
    (H5Z_func_t)H5Z_filter_sz3,   /* The actual filter function */
}};

H5PL_type_t H5PLget_plugin_type(void) {return H5PL_TYPE_FILTER;}
const void *H5PLget_plugin_info(void) {return H5Z_SZ3;}

/*FILTER FUNCTIONS*/

static herr_t H5Z_sz3_set_local(hid_t dcpl_id, hid_t type_id, hid_t chunk_space_id){

	detectSysEndianType();

	//printf("start in H5Z_sz3_set_local, dcpl_id = %d\n", dcpl_id);
	size_t r5=0,r4=0,r3=0,r2=0,r1=0, dsize;
	static char const *_funcname_ = "H5Z_sz3_set_local";
	int i, ndims, ndims_used = 0;	
	hsize_t dims[H5S_MAX_RANK], dims_used[5] = {0,0,0,0,0};	
	herr_t retval = 0;
	H5T_class_t dclass;
	H5T_sign_t dsign;
	unsigned int flags = 0;
	int dataType = 0; //SZ_FLOAT;

	//check if sz3.config in current directory
	std::ifstream f(CONFIG_PATH);
	if(f.good()){
		printf("sz3.config found!\n");
		loadConfigFile = 1;
	}
	else
		printf("sz3.config not found, using default parameters\n");

	f.close();

	//herr_t ret = H5Zregister(H5Z_SZ3);
	//printf("REGISTER: %i\n", ret);

	//printf("DC\n");
	if (0 > (dclass = H5Tget_class(type_id)))
		H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a datatype");

	//printf("DS\n");
	if (0 == (dsize = H5Tget_size(type_id)))
		H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "size is smaller than 0!");

	//printf("ND\n");
	if (0 > (ndims = H5Sget_simple_extent_dims(chunk_space_id, dims, 0)))
		H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a data space");
		
	for (i = 0; i < ndims; i++)
	{
		if (dims[i] <= 1) continue;
		dims_used[ndims_used] = dims[i];
		ndims_used++;
	}


	//printf("NDIM: %i\n", ndims);
	//printf("N_USE: %i\n", ndims_used);
	//printf("DCLASS: %i\n", dclass);
	//printf("DSIZE: %zu\n", dsize);

	//for(i = 0; i < ndims_used; i++){
	//	printf("DIMS[%i] : %zu\n", i, dims_used[i]);
	//}
	//printf("\nDCEQ\n");

	if(dclass == H5T_FLOAT)
	{
		dataType = dsize == 4 ? 0 : 1; //set float or double
	}
	else if(dclass == H5T_INTEGER)
	{
		if (0 > (dsign = H5Tget_sign(type_id)))
			H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "Error in calling H5Tget_sign(type_id)....");		
		if(dsign == H5T_SGN_NONE) //unsigned
		{
			switch(dsize)
			{
			case 1:
				dataType = SZ_UINT8;
				break;
			case 2:
				dataType = SZ_UINT16;
				break;
			case 4:
				dataType = SZ_UINT32;
				break;
			case 8:
				dataType = SZ_UINT64;
				break;
			}
		}
		else
		{
			switch(dsize)
			{
			case 1:
				dataType = SZ_INT8;
				break;
			case 2:
				dataType = SZ_INT16;
				break;
			case 4:
				dataType = SZ_INT32;
				break;
			case 8:
				dataType = SZ_INT64;
				break;
			}			
		}
	}
	else{
		H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_CANTGET, 0, "datatype class must be H5T_FLOAT or H5T_INTEGER");
	}

	//printf("GETFILT\n");

	size_t mem_cd_nelmts = 12*sizeof(unsigned int), cd_nelmts = 0;
	unsigned int *mem_cd_values = (unsigned int*) malloc(12*sizeof(unsigned int));

	//printf("MEM\n");

	if (0 > H5Pget_filter_by_id(dcpl_id, H5Z_FILTER_SZ3, &flags, &mem_cd_nelmts, mem_cd_values, 0, NULL, NULL)){
		H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_CANTGET, 0, "unable to get current SZ3 cd_values");
	}

	//for(int i = 0; i < mem_cd_nelmts; i++)
	//	printf("%u ", mem_cd_values[i]);
	//printf("\n");
	//printf("GETCD");

	freshCdValues = 0;
	switch(ndims_used)
	{
	case 1: 
		r1 = dims_used[0];
		if(mem_cd_nelmts<=4)
			freshCdValues = 1;
		break;
	case 2:
		r1 = dims_used[0];
		r2 = dims_used[1];
		if(mem_cd_nelmts<=4)
			freshCdValues = 1;
		break;
	case 3:
		r1 = dims_used[0];
		r2 = dims_used[1];
		r3 = dims_used[2];		
		if(mem_cd_nelmts<=5)
			freshCdValues = 1;
		break;
	case 4:
		r1 = dims_used[0];
		r2 = dims_used[1];
		r3 = dims_used[2];	
		r4 = dims_used[3];
		if(mem_cd_nelmts<=6)
			freshCdValues = 1;
		break;
	default: 
		H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, 0, "requires chunks w/1,2,3 or 4 non-unity dims"); 
	}

	//printf("set_local dims: %i, %i, %i, %i, %i\n", r1,r2,r3,r4,r5); 

	//printf("SETCD\n");

	if(freshCdValues)
	{
		//printf("SETCDVALS\n");
		

		unsigned int* cd_values = NULL;
		SZ_metaDataToCdArray(&cd_nelmts, &cd_values, dataType, r5, r4, r3, r2, r1);
		/* Now, update cd_values for the filter */
		if (0 > H5Pmodify_filter(dcpl_id, H5Z_FILTER_SZ3, flags, cd_nelmts, cd_values))
		{
			free(cd_values);
			free(mem_cd_values);
			H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, 0, "failed to modify cd_values");	
			
		}

		
		free(cd_values);
		free(mem_cd_values);
			
	}
	else
	{
	
		//printf("REPCDVALS\n");
		unsigned int* cd_values = NULL;
		unsigned int* final_cd_values = (unsigned int*) malloc(mem_cd_nelmts * sizeof(unsigned int));
		size_t tmp = 0;
		SZ_metaDataToCdArray(&tmp, &cd_values, dataType, r5, r4, r3, r2, r1);
		std::memcpy(final_cd_values, cd_values, tmp * sizeof(unsigned int));
		std::memcpy(final_cd_values+tmp, mem_cd_values+tmp, (mem_cd_nelmts-tmp) * sizeof(unsigned int)); 

		/*printf("%u %u\nCD\n", tmp, mem_cd_nelmts);

		for(int i =0; i < tmp; i++)
			printf("%u ", cd_values[i]);
		
		printf("\n");

		for(int i = 0; i < mem_cd_nelmts; i++)
			printf("%u ", mem_cd_values[i]);

		printf("\n");

		for(int i = 0; i < mem_cd_nelmts; i++)
			printf("%u ", final_cd_values[i]);

		printf("\n");*/

		if(0 > H5Pmodify_filter(dcpl_id, H5Z_FILTER_SZ3, flags, mem_cd_nelmts, final_cd_values))
		{

			free(final_cd_values);
			free(cd_values);
			H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, 0, "failed to modify cd_values");

		}
		
		free(mem_cd_values);
		free(final_cd_values);
		free(cd_values);

	}


	
	retval = 1;

	return retval;


}





static size_t H5Z_filter_sz3(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t* buf_size, void** buf)
{

    //store dimensions, num_dimensions, and data type
    size_t r1 = 0, r2 = 0, r3 = 0, r4 = 0, r5 = 0;
    int dimSize = 0, dataType = 0;

    //1 if error info included else 0
    int withErrInfo = checkCDValuesWithErrors(cd_nelmts, cd_values);
    int error_mode = 0;
    int cmp_algo = 1;
    int interp_algo = 1;
    float abs_error = 1E-4, rel_error = 1E-3, l2norm_error = 1E-4, psnr = 1E-3;
    if(withErrInfo)
        SZ_cdArrayToMetaDataErr(cd_nelmts, cd_values, &dimSize, &dataType, &r5, &r4, &r3, &r2, &r1, &error_mode, &abs_error, &rel_error, &l2norm_error, &psnr);
    else
        SZ_cdArrayToMetaData(cd_nelmts, cd_values, &dimSize, &dataType, &r5, &r4, &r3, &r2, &r1);

    //printf("\nwithErr: %i\n", withErrInfo);

    if(flags & H5Z_FLAG_REVERSE){
        /*decompress data*/

        SZ::Config conf;
        size_t nbEle = computeDataLength(r5,r4,r3,r2,r1);

        switch(dataType){
            case SZ_FLOAT: //FLOAT
            {
		float *f_decompressedData = new float[nbEle];
                SZ_decompress(conf, (char*) *buf, nbytes, f_decompressedData);
                free(*buf);
                *buf = f_decompressedData;
                *buf_size = nbEle * sizeof(float);
                break;
            }

            case SZ_DOUBLE: //DOUBLE
            {
                double *d_decompressedData = new double[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, d_decompressedData);
                free(*buf);
                *buf = d_decompressedData;
                *buf_size = nbEle * sizeof(double);
                break;
            }

            case SZ_INT8: //INT 8
            {
                int8_t *c_decompressedData = new int8_t[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, c_decompressedData);
                free(*buf);
                *buf = c_decompressedData;
                *buf_size = nbEle * sizeof(int8_t);
                break;
            }

            case SZ_UINT8: //UINT 8
            {
                uint8_t *uc_decompressedData = new uint8_t[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, uc_decompressedData);
                free(*buf);
                *buf = uc_decompressedData;
                *buf_size = nbEle * sizeof(uint8_t);
                break;
            }

            case SZ_INT16: //INT 16
            {
                int16_t *s_decompressedData = new int16_t[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, s_decompressedData);
                free(*buf);
                *buf = s_decompressedData;
                *buf_size = nbEle * sizeof(int16_t);
                break;
            }

            case SZ_UINT16: //UINT 16
            {
                uint16_t *us_decompressedData = new uint16_t[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, us_decompressedData);
                free(*buf);
                *buf = us_decompressedData;
                *buf_size = nbEle * sizeof(uint16_t);
                break;
            }

            case SZ_INT32: //INT 32
            {
                int32_t *i_decompressedData = new int32_t[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, i_decompressedData);
                free(*buf);
                *buf = i_decompressedData;
                *buf_size = nbEle * sizeof(int32_t);
                break;
            }

            case SZ_UINT32: //UINT 32
            {
                uint32_t *ui_decompressedData = new uint32_t[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, ui_decompressedData);
                free(*buf);
                *buf = ui_decompressedData;
                *buf_size = nbEle * sizeof(uint32_t);
                break;
            }

            case SZ_INT64: //INT 64
            {
                int64_t *l_decompressedData = new int64_t[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, l_decompressedData);
                free(*buf);
                *buf = l_decompressedData;
                *buf_size = nbEle * sizeof(int64_t);
                break;
            }

            case SZ_UINT64: //UINT 64
            {
                uint64_t *ul_decompressedData = new uint64_t[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, ul_decompressedData);
                free(*buf);
                *buf = ul_decompressedData;
                *buf_size = nbEle * sizeof(uint64_t);
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
        //based on # dimensions, get relevant dimensions and load config object with them
        if(dimSize <= 0){
            printf("Error: Number of Dimensions is <= 0");
            exit(0);
        }

        SZ::Config conf;
        //printf("\nDIMS_CMP:\n");
        //printf("r1 %u r2 %u r3 %u r4 %u r5 %u\n", r1,r2,r3,r4,r5);
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

	//if config file found and no user defined params, read the config file
	if(loadConfigFile && freshCdValues){
	//	printf("Loading sz3.config ...\n");
		conf.loadcfg(CONFIG_PATH);
	}
	else{

		conf.errorBoundMode = error_mode;
		conf.absErrorBound = abs_error;
		conf.relErrorBound = rel_error;
		conf.l2normErrorBound = l2norm_error;
		conf.psnrErrorBound = psnr;

		//printf("PARAMS: mode|%i, abs_eb|%f, rel_eb|%f, l2_eb|%f, psnr_eb|%f\n", error_mode, abs_error, rel_error, l2norm_error, psnr);

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

	}

        size_t outSize = 0;
        char* compressedData = NULL;


        switch(dataType){
            case SZ_FLOAT: //FLOAT
            {
                compressedData = SZ_compress(conf, (float*) *buf, outSize);
                break;
            }

            case SZ_DOUBLE: //DOUBLE
            {
                compressedData = SZ_compress(conf, (double*) *buf, outSize);
                break;
            }

            case SZ_INT8: //INT 8
            {
                compressedData = SZ_compress(conf, (int8_t*) *buf, outSize);
                break;
            }

            case SZ_UINT8: //UINT 8
            {
                compressedData = SZ_compress(conf, (uint8_t*) *buf, outSize);
                break;
            }

            case SZ_INT16: //INT 16
            {
                compressedData = SZ_compress(conf, (int16_t*) *buf, outSize);
                break;
            }

            case SZ_UINT16: //UINT 16
            {
                compressedData = SZ_compress(conf, (uint16_t*) *buf, outSize);
                break;
            }

            case SZ_INT32: //INT 32
            {
                compressedData = SZ_compress(conf, (int32_t*) *buf, outSize);
                break;
            }

            case SZ_UINT32: //UINT 32
            {
                compressedData = SZ_compress(conf, (uint32_t*) *buf, outSize);
                break;
            }

            case SZ_INT64: //INT 64
            {
                compressedData = SZ_compress(conf, (int64_t*) *buf, outSize);
                break;
            }

            case SZ_UINT64: //UINT 64
            {
                compressedData = SZ_compress(conf, (uint64_t*) *buf, outSize);
                break;
            }

            default:
            {
                printf("Compression Error: Unknown Datatype");
                exit(0);
            }
        }

	//printf("\nOS: %u \n", outSize);
        free(*buf);
        *buf = compressedData;
        *buf_size = outSize;


    }

    return *buf_size;

}



/*HELPER FUNCTIONS*/
//use to convert HDF5 cd_array to SZ params inside filter
void SZ_cdArrayToMetaData(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1){
    assert(cd_nelmts >= 4);
    unsigned char bytes[8];
    *dimSize = cd_values[0];
    *dataType = cd_values[1];

    switch(*dimSize)
    {
        case 1:
            SZ::int32ToBytes_bigEndian(bytes, cd_values[2]);
            SZ::int32ToBytes_bigEndian(&bytes[4], cd_values[3]);
            if(sizeof(size_t)==4)
                *r1 = (unsigned int) SZ::bytesToInt64_bigEndian(bytes);
            else
                *r1 = (unsigned long) SZ::bytesToInt64_bigEndian(bytes);
            *r2 = *r3 = *r4 = *r5 = 0;
            break;
        case 2:
            *r3 = *r4 = *r5 = 0;
            *r2 = cd_values[3];
            *r1 = cd_values[2];
            break;
        case 3:
            *r4 = *r5 = 0;
            *r3 = cd_values[4];
            *r2 = cd_values[3];
            *r1 = cd_values[2];
            break;
        case 4:
            *r5 = 0;
            *r4 = cd_values[5];
            *r3 = cd_values[4];
            *r2 = cd_values[3];
            *r1 = cd_values[2];
            break;
        default:
            *r5 = cd_values[6];
            *r4 = cd_values[5];
            *r3 = cd_values[4];
            *r2 = cd_values[3];
            *r1 = cd_values[2];
    }
}

void SZ_cdArrayToMetaDataErr(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1, int* error_bound_mode, float* abs_error, float* rel_error, float* l2norm_error, float* psnr)
{
    //get dimension, datatype metadata from cd_values
    SZ_cdArrayToMetaData(cd_nelmts, cd_values, dimSize, dataType, r5, r4, r3, r2, r1);
    //read in error bound value information
    int dim = *dimSize;
    int k = dim==1?4:dim+2;
    unsigned char b[4];
    SZ::int32ToBytes_bigEndian(b, cd_values[k++]);
    *error_bound_mode = SZ::bytesToInt32_bigEndian(b);
    SZ::int32ToBytes_bigEndian(b, cd_values[k++]);
    *abs_error = bytesToFloat(b);
    SZ::int32ToBytes_bigEndian(b, cd_values[k++]);
    *rel_error = bytesToFloat(b);
    SZ::int32ToBytes_bigEndian(b, cd_values[k++]);
    *l2norm_error = bytesToFloat(b);
    SZ::int32ToBytes_bigEndian(b, cd_values[k++]);
    *psnr = bytesToFloat(b);
}

//use to convert params for SZ into HDF5 cd_array
void SZ_metaDataToCdArray(size_t* cd_nelmts, unsigned int **cd_values, int dataType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
    *cd_values = (unsigned int*)malloc(sizeof(unsigned int)*7);
    SZ_copymetaDataToCdArray(cd_nelmts, *cd_values, dataType, r5, r4, r3, r2, r1);
}

void SZ_metaDataErrToCdArray(size_t* cd_nelmts, unsigned int **cd_values, int dataType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, int error_bound_mode, float abs_error, float rel_error, float l2norm_error, float psnr)
{
    *cd_values = (unsigned int*)malloc(sizeof(unsigned int)*12);
    SZ_copymetaDataToCdArray(cd_nelmts, *cd_values, dataType, r5, r4, r3, r2, r1);
    int k = *cd_nelmts;
    (*cd_values)[k++] = error_bound_mode;
    unsigned char b[4];
    floatToBytes(b, abs_error);
    (*cd_values)[k++] = SZ::bytesToInt32_bigEndian(b);
    floatToBytes(b, rel_error);
    (*cd_values)[k++] = SZ::bytesToInt32_bigEndian(b);
    floatToBytes(b, l2norm_error);
    (*cd_values)[k++] = SZ::bytesToInt32_bigEndian(b);
    floatToBytes(b, psnr);
    (*cd_values)[k++] = SZ::bytesToInt32_bigEndian(b);
    *cd_nelmts = k;
}

void SZ_copymetaDataToCdArray(size_t* cd_nelmts, unsigned int *cd_values, int dataType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
    unsigned char bytes[8] = {0};
    unsigned long size;
    int dim = computeDimension(r5, r4, r3, r2, r1);
    cd_values[0] = dim;
    cd_values[1] = dataType;	//0: FLOAT ; 1: DOUBLE ; 2,3,4,....: INTEGER....

    switch(dim)
    {
        case 1:
            size = (unsigned long)r1;
            SZ::int64ToBytes_bigEndian(bytes, size);
            cd_values[2] = SZ::bytesToInt32_bigEndian(bytes);
            cd_values[3] = SZ::bytesToInt32_bigEndian(&bytes[4]);
            *cd_nelmts = 4;
            break;
        case 2:
            cd_values[2] = (unsigned int) r2;
            cd_values[3] = (unsigned int) r1;
            *cd_nelmts = 4;
            break;
        case 3:
            cd_values[2] = (unsigned int) r3;
            cd_values[3] = (unsigned int) r2;
            cd_values[4] = (unsigned int) r1;
            *cd_nelmts = 5;
            break;
        case 4:
            cd_values[2] = (unsigned int) r4;
            cd_values[3] = (unsigned int) r3;
            cd_values[4] = (unsigned int) r2;
            cd_values[5] = (unsigned int) r1;
            *cd_nelmts = 6;
            break;
        default:
            cd_values[2] = (unsigned int) r5;
            cd_values[3] = (unsigned int) r4;
            cd_values[4] = (unsigned int) r3;
            cd_values[5] = (unsigned int) r2;
            cd_values[6] = (unsigned int) r1;
            *cd_nelmts = 7;
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
            if(cd_nelmts>4)
                result = 1;
            break;
        case 2:
            if(cd_nelmts>4)
                result = 1;
            break;
        case 3:
            if(cd_nelmts>5)
                result = 1;
            break;
        case 4:
            if(cd_nelmts>6)
                result = 1;
            break;
        case 5:
            if(cd_nelmts>7)
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

//detect sys endian type
void detectSysEndianType()
{
    //get sys endian type
    int x_temp = 1;
    char *y_temp = (char*)&x_temp;

    if(*y_temp==1)
        sysEndianType = LITTLE_ENDIAN_SYSTEM;
    else //=0
        sysEndianType = BIG_ENDIAN_SYSTEM;
}


//the byte to input is in the big-endian format
inline float bytesToFloat(unsigned char* bytes)
{
	lfloat buf;
		memcpy(buf.byte, bytes, 4);
			if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
					symTransform_4bytes(buf.byte);	
						return buf.value;
						}

						inline void floatToBytes(unsigned char *b, float num)
						{
							lfloat buf;
								buf.value = num;
									memcpy(b, buf.byte, 4);
										if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
												symTransform_4bytes(b);								}
