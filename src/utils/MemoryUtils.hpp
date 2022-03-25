/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
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
 * @file    MemoryUtils.hpp
 * @authors tlin crdrisko
 * @date    10/25/2004
 * @version 1.0
 */

#ifndef OPENMD_UTILS_MEMORYUTILS_HPP
#define OPENMD_UTILS_MEMORYUTILS_HPP

#include <memory>
#include <type_traits>
#include <utility>

namespace OpenMD::Utils {
  namespace details {

    template<typename, typename = std::void_t<>>
    struct is_container : std::false_type {};

    template<typename T>
    struct is_container<T, std::void_t<typename T::value_type,
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
                                       decltype(std::declval<T>().empty())>> : std::true_type {};

    // Convenience variable template for ease-of-use
    template<typename T>
    constexpr bool is_container_v = is_container<T>::value;

    template<typename T>
    void lifted_deleter(T val) {
      delete val;
    }

    template<typename T1, typename T2>
    void lifted_deleter(std::pair<T1, T2>& val) {
      lifted_deleter(val.second);
    }
  }  // namespace details

  // Recursively deallocate the previously allocated memory from each container
  template<typename Container>
  void deletePointers(Container& container) {
    if constexpr(details::is_container_v<typename Container::value_type>) {
      for (auto& elem : container)
        deletePointers(elem);
    } else {
      for (auto& elem : container)
        details::lifted_deleter(elem);
    }
  }
}  // namespace OpenMD::Utils

#endif  // OPENMD_UTILS_MEMORYUTILS_HPP
