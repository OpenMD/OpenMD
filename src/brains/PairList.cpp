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

#include "brains/PairList.hpp"

#include <functional>
#include <iterator>
#include <utility>

namespace OpenMD {

  int* PairList::getPairList() {
    if (modified_) {
      pairList_.clear();

      for (std::set<std::pair<int, int>>::iterator i = pairSet_.begin();
           i != pairSet_.end(); ++i) {
        pairList_.push_back(i->first + 1);
        pairList_.push_back(i->second + 1);
      }
      modified_ = false;
    }

    return pairList_.empty() ? NULL : &(pairList_[0]);
  }

  void PairList::addPair(int i, int j) {
    if (i == j) {
      return;
    } else if (i > j) {
      std::swap(i, j);
    }

    std::set<std::pair<int, int>>::iterator iter =
        pairSet_.find(std::make_pair(i, j));

    if (iter == pairSet_.end()) {
      pairSet_.insert(std::make_pair(i, j));
      modified_ = true;
    }
  }

  void PairList::addPairs(std::set<int>& set1, std::set<int>& set2) {
    for (std::set<int>::iterator iter1 = set1.begin(); iter1 != set1.end();
         ++iter1) {
      for (std::set<int>::iterator iter2 = set2.begin(); iter2 != set2.end();
           ++iter2) {
        this->addPair(*iter1, *iter2);
      }
    }
  }

  template<typename IterType1, typename IterType2>
  void PairList::addPairs(IterType1 iter1_first, IterType1 iter1_last,
                          IterType2 iter2_first, IterType2 iter2_last) {
    for (IterType1 iter1 = iter1_first; iter1 != iter1_last; ++iter1) {
      for (IterType2 iter2 = iter2_first; iter2 != iter2_last; ++iter2) {
        this->addPair(*iter1, *iter2);
      }
    }
  }

  void PairList::removePair(int i, int j) {
    if (i == j) {
      return;
    } else if (i > j) {
      std::swap(i, j);
    }

    std::set<std::pair<int, int>>::iterator iter =
        pairSet_.find(std::make_pair(i, j));

    if (iter != pairSet_.end()) {
      pairSet_.erase(iter);
      modified_ = true;
    }
  }

  void PairList::removePairs(std::set<int>& set1, std::set<int>& set2) {
    for (std::set<int>::iterator iter1 = set1.begin(); iter1 != set1.end();
         ++iter1) {
      for (std::set<int>::iterator iter2 = set2.begin(); iter2 != set2.end();
           ++iter2) {
        this->removePair(*iter1, *iter2);
      }
    }
  }

  template<typename IterType1, typename IterType2>
  void PairList::removePairs(IterType1 iter1_first, IterType1 iter1_last,
                             IterType2 iter2_first, IterType2 iter2_last) {
    for (IterType1 iter1 = iter1_first; iter1 != iter1_last; ++iter1) {
      for (IterType2 iter2 = iter2_first; iter2 != iter2_last; ++iter2) {
        this->removePair(*iter1, *iter2);
      }
    }
  }

  bool PairList::hasPair(int i, int j) {
    if (i == j) {
      return false;
    } else if (i > j) {
      std::swap(i, j);
    }

    std::set<std::pair<int, int>>::iterator iter =
        pairSet_.find(std::make_pair(i, j));
    return iter == pairSet_.end() ? false : true;
  }

  int PairList::getSize() { return pairSet_.size(); }

  std::ostream& operator<<(std::ostream& o, PairList& e) {
    std::set<std::pair<int, int>>::iterator i;

    int index;

    index = 0;

    for (i = e.pairSet_.begin(); i != e.pairSet_.end(); ++i) {
      o << "pairList[" << index << "] i, j: " << (*i).first << " - "
        << (*i).second << "\n";
      index++;
    }

    return o;
  }

}  // namespace OpenMD
