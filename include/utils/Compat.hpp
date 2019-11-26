#ifndef SZ_COMPAT_HEADER
#define SZ_COMPAT_HEADER
#include <memory>
#include <type_traits>

namespace SZ {
namespace compat {
#ifndef __cpp_lib_make_unique
  /**
   * early adopt std::make_unique from C++14.  This version is adapted from LLVM's libc++.
   * It was modified to not use names reserved for the standard.
   * the license header for this file is:
   *
   * // Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
   * // See https://llvm.org/LICENSE.txt for license information.
   * // SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
   */
  namespace detail {
    template<class Tp>
    struct unique_if
    {
        typedef std::unique_ptr<Tp> unique_single;
    };

    template<class Tp>
    struct unique_if<Tp[]>
    {
        typedef std::unique_ptr<Tp[]> unique_array_unknown_bound;
    };

    template<class Tp, size_t Np>
    struct unique_if<Tp[Np]>
    {
        typedef void unique_array_known_bound;
    };
  }

  template<class Tp, class... Args>
  inline 
  typename detail::unique_if<Tp>::unique_single
  make_unique(Args&&... args)
  {
      return std::unique_ptr<Tp>(new Tp(std::forward<Args>(args)...));
  }

  template<class Tp>
  inline
  typename detail::unique_if<Tp>::unique_array_unknown_bound
  make_unique(size_t n)
  {
      typedef typename std::remove_extent<Tp>::type Up;
      return std::unique_ptr<Tp>(new Up[n]());
  }

  template<class Tp, class... Args>
      typename detail::unique_if<Tp>::unique_array_known_bound
      make_unique(Args&&...) = delete;
#else
  using std::make_unique;
#endif

}
}

#endif
