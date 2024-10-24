/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
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
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

/**
 * @file LocalIndexManager.hpp
 * @author    tlin
 * @date  11/02/2004
 * @version 1.0
 */

#ifndef UTILS_LOCALINDEXMANAGER_HPP
#define UTILS_LOCALINDEXMANAGER_HPP

#include <algorithm>
#include <cassert>
#include <iostream>
#include <list>
#include <utility>

namespace OpenMD {

  /**
   * @class IndexListContainer
   * @brief
   * @todo documentation
   */
  class IndexListContainer {
  public:
    static const int MAX_INTEGER = 2147483647;

    using IndexListContainerIterator = std::list<std::pair<int, int>>::iterator;

    IndexListContainer(int minIndex = 0, int maxIndex = MAX_INTEGER) :
        maxIndex_(maxIndex) {
      indexContainer_.push_back(std::make_pair(minIndex, maxIndex));
    }

    int pop() {
      if (indexContainer_.empty()) { std::cerr << "" << std::endl; }

      int result;

      result = indexContainer_.front().first;

      if (indexContainer_.front().first == indexContainer_.front().second) {
        indexContainer_.pop_front();
      } else if (indexContainer_.front().first <
                 indexContainer_.front().second) {
        ++indexContainer_.front().first;
      } else {
        std::cerr << "" << std::endl;
      }

      return result;
    }

    /**
     *
     */
    void insert(int index) {
      IndexListContainerIterator insertPos = internalInsert(index, index);
      merge(insertPos);
    }

    /**
     * Reclaims an index range
     * @param beginIndex
     * @param endIndex
     */
    void insert(int beginIndex, int endIndex) {
      IndexListContainerIterator insertPos =
          internalInsert(beginIndex, endIndex);
      merge(insertPos);
    }

    /**
     * Reclaims an index array.
     * @param indices
     */
    void insert(std::vector<int>& indices) {
      if (indices.empty()) { return; }

      std::sort(indices.begin(), indices.end());
      auto last = std::unique(indices.begin(), indices.end());
      indices.erase(last, indices.end());

      std::vector<int>::iterator i;
      int beginIndex;

      beginIndex = indices[0];

      for (i = indices.begin() + 1; i != indices.end(); ++i) {
        if (*i != *(i - 1) + 1) {
          insert(beginIndex, *(i - 1));
          beginIndex = *i;
        }
      }
    }

    std::vector<int> getIndicesBefore(int index) {
      std::vector<int> result;
      IndexListContainerIterator i;

      for (i = indexContainer_.begin(); i != indexContainer_.end(); ++i) {
        if ((*i).first > index) {
          // we locate the node whose minimum index is greater that index
          indexContainer_.erase(indexContainer_.begin(), i);
          break;
        } else if ((*i).second < index) {
          // The biggest index current node hold is less than the index we want
          for (int j = (*i).first; j <= (*i).second; ++j) {
            result.push_back(j);
          }
          continue;
        } else if ((*i).first == (*i).second) {
          // the index happen to equal to a node which only contains one index
          result.push_back((*i).first);
          indexContainer_.erase(indexContainer_.begin(), i);
          break;
        } else {
          for (int j = (*i).first; j < index; ++j) {
            result.push_back(j);
          }
          (*i).first = index;
          indexContainer_.erase(indexContainer_.begin(), i);
          break;
        }
      }
      return result;
    }

    int getMaxIndex() { return maxIndex_; }

  private:
    IndexListContainerIterator internalInsert(int beginIndex, int endIndex) {
      if (beginIndex > endIndex) {
        std::swap(beginIndex, endIndex);
        std::cerr << "" << std::endl;
      }

      if (endIndex > maxIndex_) { std::cerr << "" << std::endl; }

      IndexListContainerIterator i = indexContainer_.begin();
      for (; i != indexContainer_.end(); ++i) {
        if ((*i).first > endIndex) {
          indexContainer_.insert(i, std::make_pair(beginIndex, endIndex));
          return --i;
        } else if ((*i).second < beginIndex) {
          continue;
        } else {
          std::cerr << "" << std::endl;
        }
      }

      indexContainer_.push_back(std::make_pair(beginIndex, endIndex));
      return --indexContainer_.end();
    }

    void merge(IndexListContainerIterator i) {
      IndexListContainerIterator j;

      // check whether current node can be merged with its previous node
      if (i != indexContainer_.begin()) {
        j = i;
        --j;
        if (j != indexContainer_.begin() && (*j).second + 1 == (*i).first) {
          (*i).first = (*j).first;
          indexContainer_.erase(j);
        }
      }

      // check whether current node can be merged with its next node
      if (i != indexContainer_.end()) {
        j = i;
        ++j;

        if (j != indexContainer_.end() && (*i).second + 1 == (*j).first) {
          (*i).second = (*j).second;
          indexContainer_.erase(j);
        }
      }
    }

    int maxIndex_;
    std::list<std::pair<int, int>> indexContainer_;
  };

  /**
   * @class LocalIndexManager LocalIndexManager.hpp
   * "utils/LocalIndexManager.hpp"
   * @brief
   */
  class LocalIndexManager {
  public:
    int getNextAtomIndex() { return atomIndexContainer_.pop(); }

    std::vector<int> getAtomIndicesBefore(int index) {
      return atomIndexContainer_.getIndicesBefore(index);
    }

    void releaseAtomIndex(int index) { atomIndexContainer_.insert(index); }

    void releaseAtomIndex(int beginIndex, int endIndex) {
      atomIndexContainer_.insert(beginIndex, endIndex);
    }

    void releaseAtomIndex(std::vector<int> indices) {
      atomIndexContainer_.insert(indices);
    }

    int getNextRigidBodyIndex() { return rigidBodyIndexContainer_.pop(); }

    std::vector<int> getRigidBodyIndicesBefore(int index) {
      return rigidBodyIndexContainer_.getIndicesBefore(index);
    }

    void releaseRigidBodyIndex(int index) {
      rigidBodyIndexContainer_.insert(index);
    }

    void releaseRigidBodyIndex(int beginIndex, int endIndex) {
      rigidBodyIndexContainer_.insert(beginIndex, endIndex);
    }

    void releaseRigidBodyIndex(std::vector<int> indices) {
      rigidBodyIndexContainer_.insert(indices);
    }

    int getNextCutoffGroupIndex() { return cutoffGroupIndexContainer_.pop(); }

    std::vector<int> getCutoffGroupIndicesBefore(int index) {
      return cutoffGroupIndexContainer_.getIndicesBefore(index);
    }

    void releaseCutoffGroupIndex(int index) {
      cutoffGroupIndexContainer_.insert(index);
    }

    void releaseCutoffGroupIndex(int beginIndex, int endIndex) {
      cutoffGroupIndexContainer_.insert(beginIndex, endIndex);
    }

    void releaseCutoffGroupIndex(std::vector<int> indices) {
      cutoffGroupIndexContainer_.insert(indices);
    }

    int getNextBondIndex() { return bondIndexContainer_.pop(); }

    std::vector<int> getBondIndicesBefore(int index) {
      return bondIndexContainer_.getIndicesBefore(index);
    }

    void releaseBondIndex(int index) { bondIndexContainer_.insert(index); }

    void releaseBondIndex(int beginIndex, int endIndex) {
      bondIndexContainer_.insert(beginIndex, endIndex);
    }

    void releaseBondIndex(std::vector<int> indices) {
      bondIndexContainer_.insert(indices);
    }

    int getNextBendIndex() { return bendIndexContainer_.pop(); }

    std::vector<int> getBendIndicesBefore(int index) {
      return bendIndexContainer_.getIndicesBefore(index);
    }

    void releaseBendIndex(int index) { bendIndexContainer_.insert(index); }

    void releaseBendIndex(int beginIndex, int endIndex) {
      bendIndexContainer_.insert(beginIndex, endIndex);
    }

    void releaseBendIndex(std::vector<int> indices) {
      bendIndexContainer_.insert(indices);
    }

    int getNextTorsionIndex() { return torsionIndexContainer_.pop(); }

    std::vector<int> getTorsionIndicesBefore(int index) {
      return torsionIndexContainer_.getIndicesBefore(index);
    }

    void releaseTorsionIndex(int index) {
      torsionIndexContainer_.insert(index);
    }

    void releaseTorsionIndex(int beginIndex, int endIndex) {
      torsionIndexContainer_.insert(beginIndex, endIndex);
    }

    void releaseTorsionIndex(std::vector<int> indices) {
      torsionIndexContainer_.insert(indices);
    }

    int getNextInversionIndex() { return inversionIndexContainer_.pop(); }

    std::vector<int> getInversionIndicesBefore(int index) {
      return inversionIndexContainer_.getIndicesBefore(index);
    }

    void releaseInversionIndex(int index) {
      inversionIndexContainer_.insert(index);
    }

    void releaseInversionIndex(int beginIndex, int endIndex) {
      inversionIndexContainer_.insert(beginIndex, endIndex);
    }

    void releaseInversionIndex(std::vector<int> indices) {
      inversionIndexContainer_.insert(indices);
    }

  private:
    IndexListContainer atomIndexContainer_;
    IndexListContainer rigidBodyIndexContainer_;
    IndexListContainer cutoffGroupIndexContainer_;
    IndexListContainer bondIndexContainer_;
    IndexListContainer bendIndexContainer_;
    IndexListContainer torsionIndexContainer_;
    IndexListContainer inversionIndexContainer_;
  };
}  // namespace OpenMD

#endif  // UTILS_LOCALINDEXMANAGER_HPP
