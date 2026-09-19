#if INTPTR_MAX == INT64_MAX  // use ska for 64bit system
#include "SZ3/utils/ska_hash/unordered_map.hpp"
#else   // most likely 32bit system
#include <unordered_map>
#endif  // INTPTR_MAX == INT64_MAX

namespace SZ3 {

// SZ3_USE_SKA_HASH is an INTERFACE define from CMakeLists.txt; a consumer that includes this header
// without it gets std::unordered_map, which is what an undefined macro already selected.
#if defined(SZ3_USE_SKA_HASH) && (SZ3_USE_SKA_HASH) && (INTPTR_MAX == INT64_MAX)  // use ska for 64bit system
    using ska::unordered_map;
#else   // most likely 32bit system
    using std::unordered_map;
#endif  // INTPTR_MAX == INT64_MAX

}
