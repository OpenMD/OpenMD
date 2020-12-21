/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

/**
 * @file MemoryUtils.hpp
 * @author    tlin
 * @date  10/25/2004
 * @version 1.0
 */

#ifndef UTILS_MEMORYUTILS_HPP
#define UTILS_MEMORYUTILS_HPP

#include <memory>
#include <type_traits>
#include <utility>

namespace OpenMD {

  namespace MemoryUtils {

    // Remove in favor of std::make_unique<> when we switch to C++14
    template<typename T, typename... TArgs,
             typename = typename std::enable_if<!std::is_array<T>::value>::type>
    std::unique_ptr<T> make_unique(TArgs&&... args) {
      return std::unique_ptr<T> {new T(std::forward<TArgs>(args)...)};
    }

    namespace details {

      // Remove in favor of std::void_t<> when we switch to C++17
      template<typename...>
      using void_t = void;

      template<typename, typename = void_t<>>
      struct is_container : std::integral_constant<bool, false> {};

      template<typename T>
      struct is_container<T, void_t<typename T::value_type,
                                    typename T::reference,
                                    typename T::const_reference,
                                    typename T::iterator,
                                    typename T::const_iterator,
                                    typename T::difference_type,
                                    typename T::size_type,
                                    decltype(std::declval<T>().begin()),
                                    decltype(std::declval<T>().end()),
                                    decltype(std::declval<T>().cbegin()),
                                    decltype(std::declval<T>().cend()),
                                    decltype(std::declval<T>().max_size()),
                                    decltype(std::declval<T>().empty())>> : std::integral_constant<bool, true> {};

      template<typename T>
      void lifted_deleter(T val) {
        delete val;
      }

      template<typename T1, typename T2>
      void lifted_deleter(std::pair<T1, T2>& val) {
        lifted_deleter(val.second);
      }
    }

    // Base case - Combine using constexpr-if when we switch to C++17
    template<typename Container,
             typename = typename std::enable_if<!details::is_container<typename Container::value_type>::value>::type>
    void deletePointers(Container& container) {

      for (auto& elem : container)
        details::lifted_deleter(elem);
    }

    // Recursively deallocate the previously allocated memory from each container
    template<typename Container,
             typename = typename std::enable_if<details::is_container<typename Container::value_type>::value>::type,
             typename = int /* Dummy parameter to allow multiple default template parameters in overloads */>
    void deletePointers(Container& containerOfContainers) {

      for (auto& container : containerOfContainers)
        deletePointers(container);
    }
  }
}

#endif // UTILS_MEMORYUTILS_HPP
