#ifndef _DEF_HPP
#define _DEF_HPP

#include <vector>
#include <string>
#include <cmath>

/*------   Version   ------*/
#define SZ_VERSION_MAJOR    3
#define SZ_VERSION_MINOR    1
#define SZ_VERSION_RELEASE  0

#define SZ_LIB_VERSION SZ_VERSION_MAJOR.SZ_VERSION_MINOR.SZ_VERSION_RELEASE
#define SZ_QUOTE(str) #str
#define SZ_EXPAND_AND_QUOTE(str) SZ_QUOTE(str)
#define SZ_VERSION_STRING SZ_EXPAND_AND_QUOTE(SZ_LIB_VERSION)

const char *SZ_versionString() { return SZ_VERSION_STRING; };

namespace SZ {

    typedef unsigned int uint;
    typedef unsigned char uchar;
}


#define ABS 0
#define REL 1
//#define ABS_AND_REL 2
//#define ABS_OR_REL 3
//#define PSNR 4
//#define NORM 5
//#define PW_REL 10
//#define ABS_AND_PW_REL 11
//#define ABS_OR_PW_REL 12
//#define REL_AND_PW_REL 13
//#define REL_OR_PW_REL 14

#define METHOD_LORENZO_REG 0
#define METHOD_INTERP_LORENZO 1
#define METHOD_INTERP 2

#endif
