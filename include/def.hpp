#ifndef _DEF_HPP
#define _DEF_HPP

#include <vector>
#include <string>
#include <cmath>

/*------   Version   ------*/
#define SZ_VERSION_MAJOR    3
#define SZ_VERSION_MINOR    0
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
#endif
