#ifndef _SZ_PREDICTOR_HPP
#define _SZ_PREDICTOR_HPP

#include "utils/Concepts.hpp"
#include "LorenzoPredictor.hpp"

namespace SZ {
  namespace concepts {

    template <typename T, typename = void>
    struct is_predictor : false_type{};
    template <typename T>
    struct is_predictor<T, void_t<
      typename T::iterator,
      decltype(std::declval<T>().predict(std::declval<typename T::iterator const>())),
      decltype(std::declval<T>().preprocess(std::declval<typename T::iterator const>())),
      decltype(std::declval<T>().postprocess(std::declval<typename T::iterator const>()))
      >> : true_type{
        //we must remove_reference otherwise we get the const-ness of the reference not the underlying type
        static_assert(
          std::is_const<typename std::remove_reference<
              decltype(*std::declval<typename T::iterator const>())
            >::type
          >::value, "const iterators must not be writable");
      };

  }
}

#endif
