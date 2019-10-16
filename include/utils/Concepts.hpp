#include <type_traits>
#include <utility>

namespace SZ {
  namespace concepts {
    template<typename... Ts> struct make_void { typedef void type;};
    template<typename... Ts> using void_t = typename make_void<Ts...>::type;
    typedef std::integral_constant<bool, true>  true_type;
    typedef std::integral_constant<bool, false> false_type;

  }
}
