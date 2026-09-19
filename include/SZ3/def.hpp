#ifndef SZ3_DEF_HPP
#define SZ3_DEF_HPP

namespace SZ3 {

typedef unsigned int uint;
typedef unsigned char uchar;
#define SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH \
    "The buffer for compressed data is not large enough."

namespace concepts {
// Carries the virtual destructor for the interfaces below, so that none of them declares one of its
// own: a declared destructor deprecates the implicit copy, and the stages are copied by value.
class Interface {
   public:
    virtual ~Interface() = default;
    Interface() = default;
    Interface(const Interface &) = default;
    Interface &operator=(const Interface &) = default;
};
}  // namespace concepts
}  // namespace SZ3

#ifdef _MSC_VER
#define ALWAYS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define ALWAYS_INLINE inline
#endif

#endif
