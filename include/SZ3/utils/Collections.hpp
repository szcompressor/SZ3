/**
 * @file Collections.hpp
 * @ingroup UtilsStructure
 * @brief Wrapper for hash map implementation.
 */

#if INTPTR_MAX == INT64_MAX  // use ska for 64bit system
#include "SZ3/utils/thirdparty/ska_hash/unordered_map.hpp"
#else   // most likely 32bit system
#include <unordered_map>
#endif  // INTPTR_MAX == INT64_MAX

namespace SZ3 {

#if (SZ3_USE_SKA_HASH) && (INTPTR_MAX == INT64_MAX)  // use ska for 64bit system
    using ska::unordered_map;
#else   // most likely 32bit system
    using std::unordered_map;
#endif  // INTPTR_MAX == INT64_MAX

}
