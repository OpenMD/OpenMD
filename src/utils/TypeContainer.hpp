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
 * @file TypeContainer.hpp
 * @author tlin
 * @date 11/04/2004
 * @version 1.0
 */

#ifndef UTILS_TYPECONTAINER_HPP
#define UTILS_TYPECONTAINER_HPP

#include <algorithm>
#include <map>
#include <vector>

#include "utils/MemoryUtils.hpp"
#include "utils/Utility.hpp"

namespace OpenMD {

  /**
   * @class TypeContainer TypeContainer.hpp "utils/TypeContainer.hpp"
   */
  template<class ElemType, int SIZE>
  class TypeContainer {
  public:
    using ElemPtr         = ElemType*;
    using KeyType         = std::vector<std::string>;
    using KeyTypeIterator = typename KeyType::iterator;
    using ValueType       = std::pair<int, ElemPtr>;
    using MapType         = typename std::map<KeyType, ValueType>;
    using MapTypeIterator = typename std::map<KeyType, ValueType>::iterator;
    using value_type      = typename MapType::value_type;
    using MutableValues   = typename std::vector<int>;

    TypeContainer() : index_(0) {}

    ~TypeContainer() { Utils::deletePointers(data_); }

    bool add(KeyType& keys, ElemPtr elem) {
      assert(keys.size() == SIZE);
      assert(elem);
      return data_.insert(value_type(keys, std::make_pair(index_++, elem)))
          .second;
    }

    bool replace(KeyType& keys, ElemPtr elem) {
      assert(keys.size() == SIZE);
      assert(elem);

      MapTypeIterator i;
      i = data_.find(keys);
      if (i != data_.end()) {
        data_[keys] = std::make_pair((i->second).first, elem);
        return true;
      } else {
        return data_.insert(value_type(keys, std::make_pair(index_++, elem)))
            .second;
      }
    }

    /** Exact Match */
    ElemPtr find(KeyType& keys) {
      assert(keys.size() == SIZE);
      MapTypeIterator i;

      i = data_.find(keys);
      if (i != data_.end()) { return (i->second).second; }

      KeyType reversedKeys = keys;
      std::reverse(reversedKeys.begin(), reversedKeys.end());

      i = data_.find(reversedKeys);
      if (i != data_.end()) {
        return (i->second).second;
      } else {
        return NULL;
      }
    }

    ElemPtr permutedFindSkippingFirstElement(KeyType& keys) {
      assert(keys.size() == SIZE);
      MapTypeIterator i;

      KeyType permutedKeys = keys;

      // skip the first element:
      KeyTypeIterator start;
      start = permutedKeys.begin();
      ++start;

      std::sort(start, permutedKeys.end());

      do {
        i = data_.find(permutedKeys);
        if (i != data_.end()) { return (i->second).second; }
      } while (std::next_permutation(start, permutedKeys.end()));

      return NULL;
    }

    /**
     * @todo
     */
    ElemPtr find(KeyType& keys, const std::string& wildCard) {
      assert(keys.size() == SIZE);

      std::vector<KeyTypeIterator> iterCont;
      KeyType replacedKey;
      MapTypeIterator i;
      std::vector<ValueType> foundTypes;

      while (replaceWithWildCard(iterCont, keys, replacedKey, wildCard)) {
        i = data_.find(replacedKey);
        if (i != data_.end()) { foundTypes.push_back(i->second); }
      }

      // reverse the order of keys
      KeyType reversedKeys = keys;
      std::reverse(reversedKeys.begin(), reversedKeys.end());

      // if the reversedKeys is the same as keys, just skip it
      if (reversedKeys != keys) {
        // empty the iterator container
        iterCont.clear();

        while (replaceWithWildCard(iterCont, reversedKeys, replacedKey,
                                   wildCard)) {
          i = data_.find(replacedKey);
          if (i != data_.end()) { foundTypes.push_back(i->second); }
        }

        // replaceWithWildCard can not generate this particular sequence, we
        // have to do it manually
        KeyType allWildCards(SIZE, wildCard);
        i = data_.find(replacedKey);
        if (i != data_.end()) { foundTypes.push_back(i->second); }
      }

      typename std::vector<ValueType>::iterator j;
      j = std::min_element(foundTypes.begin(), foundTypes.end());

      return j == foundTypes.end() ? NULL : j->second;
    }

    size_t size() { return data_.size(); }

    ElemPtr beginType(MapTypeIterator& i) {
      i = data_.begin();
      return i != data_.end() ? (i->second).second : NULL;
    }

    ElemPtr nextType(MapTypeIterator& i) {
      ++i;
      return i != data_.end() ? (i->second).second : NULL;
    }

    KeyType getKeys(MapTypeIterator& i) { return i->first; }

  private:
    int index_;
    MapType data_;
  };
}  // namespace OpenMD

#endif  // UTILS_TYPECONTAINER_HPP
