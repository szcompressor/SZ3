#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#define ABS 0
#define REL 1
#define NORM2 2
#define PSNR 3


#define LITTLE_ENDIAN_SYSTEM 0
#define BIG_ENDIAN_SYSTEM 1
#define LITTLE_ENDIAN_DATA 0
#define BIG_ENDIAN_DATA 1


int sysEndianType = LITTLE_ENDIAN_SYSTEM;
int dataEndianType = LITTLE_ENDIAN_DATA;

typedef union lfloat
{
    float value;
    unsigned int ivalue;
    unsigned char byte[4];
} lfloat;

lfloat buf;

void usage()
{
	printf("Usage: print_h5repack_args <options>\n");
	printf("Options:\n");
	printf("	-M <error bound mode> : 10 options as follows. \n");
	printf("		ABS (absolute error bound)\n");
	printf("		REL (value range based error bound, so a.k.a., VR_REL)\n");
	printf("		PSNR (peak signal-to-noise ratio)\n");
	printf("		NORM2 (norm2)\n");
	printf("	-A <absolute error bound>: specifying absolute error bound\n");
	printf("	-R <value_range based relative error bound>: specifying relative error bound\n");
	printf("	-N <norm2>: specifying norm2 error bound\n");
	printf("	-S <PSNR>: specifying PSNR\n");
	printf("* examples: \n");
	printf("	print_h5repack_args -M ABS -A 1E-3 (output: -f UD=32024,0,5,0,981668463,0,0,0)\n");
	printf("	print_h5repack_args -M REL -R 1E-4 (output: -f UD=32024,0,5,1,0,953267991,0,0)\n");
	exit(0);
}

void symTransform_4bytes(unsigned char data[4])
{
	unsigned char tmp = data[0];
	data[0] = data[3];
	data[3] = tmp;

	tmp = data[1];
	data[1] = data[2];
	data[2] = tmp;
}

void floatToBytes(unsigned char *b, float num)
{
	lfloat buf;
	buf.value = num;
	memcpy(b, buf.byte, 4);
	if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
		symTransform_4bytes(b);		
}

int bytesToInt32_bigEndian(unsigned char* bytes)
{
	int temp = 0;
	int res = 0;

	res <<= 8;
	temp = bytes[0] & 0xff;
	res |= temp;

	res <<= 8;
	temp = bytes[1] & 0xff;
	res |= temp;

	res <<= 8;
	temp = bytes[2] & 0xff;
	res |= temp;

	res <<= 8;
	temp = bytes[3] & 0xff;
	res |= temp;

	return res;
}

int main(int argc, char* argv[])
{
	char* errBoundMode = NULL;
	char* absErrorBound = NULL;
	char* relErrorBound = NULL;
	char* norm2ErrorBound = NULL;
	char* psnr_ = NULL;
	
	if(argc==1)
		usage();
	
	int i = 0;
	
	for(i=1;i<argc;i++)
	{
		if (argv[i][0] != '-' || argv[i][2])
			usage();
		switch (argv[i][1])
		{
		case 'h':
			usage();
			exit(0);
			break;
		case 'M':
			if (++i == argc)
				usage();
			errBoundMode = argv[i];
			break;
		case 'A':
			if (++i == argc)
				usage();
			absErrorBound = argv[i];
			break;
		case 'R':
			if (++i == argc)
				usage();
			relErrorBound = argv[i];
			break;
		case 'N':
			if (++i == argc)
				usage();
			norm2ErrorBound = argv[i];
			break;
		case 'S': 
			if (++i == argc)
				usage();
			psnr_ = argv[i];
			break;
		default: 
			usage();
			break;
		}
	}
	unsigned char b[4];
	unsigned int cd_values[5];
	if(strcmp(errBoundMode, "ABS")==0)
	{
		cd_values[0] = ABS;
		floatToBytes(b, atof(absErrorBound));
		cd_values[1]  = bytesToInt32_bigEndian(b);
		cd_values[2] = 0;
		cd_values[3] = 0;
		cd_values[4] = 0;
	}
	else if(strcmp(errBoundMode, "REL")==0)
	{
		cd_values[0] = REL;
		cd_values[1] = 0;
		floatToBytes(b, atof(relErrorBound));
		cd_values[2]  = bytesToInt32_bigEndian(b);
		cd_values[3] = 0;
		cd_values[4] = 0;		
	}
	else if(strcmp(errBoundMode, "NORM2")==0)
	{
		cd_values[0] = NORM2;
		cd_values[1] = 0;
		cd_values[2] = 0;
		floatToBytes(b, atof(norm2ErrorBound));
		cd_values[3]  = bytesToInt32_bigEndian(b);
		cd_values[4] = 0;		
	}
	else if(strcmp(errBoundMode, "PSNR")==0)
	{
		cd_values[0] = PSNR;
		cd_values[1] = 0;
		cd_values[2] = 0;
		cd_values[3] = 0;
		floatToBytes(b, atof(psnr_));
		cd_values[4]  = bytesToInt32_bigEndian(b);
	}
	else
	{
		printf("Error: wrong errBoundMode setting.\n");
		exit(0);
	}
	printf("-f UD=32024,0,5");

	for(i=0;i<5;i++)
		printf(",%d",cd_values[i]);
	printf("\n");
}
