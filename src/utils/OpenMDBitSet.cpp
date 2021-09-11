/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
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

#include "utils/OpenMDBitSet.hpp"

#include <algorithm>
#include <cassert>
#include <functional>
#include <iterator>
#include <string>

#ifdef IS_MPI
#include <mpi.h>
#endif

namespace OpenMD {

  int OpenMDBitSet::countBits() {
    return std::count(bitset_.begin(), bitset_.end(), true);
  }

  void OpenMDBitSet::flip(int fromIndex, int toIndex) {
    assert(fromIndex <= toIndex);
    assert(fromIndex >= 0);
    assert(toIndex <= size());
    std::vector<bool>::iterator first = bitset_.begin() + fromIndex;
    std::vector<bool>::iterator last  = bitset_.begin() + toIndex;

    std::transform(first, last, first, std::logical_not<bool>());
  }

  OpenMDBitSet OpenMDBitSet::get(int fromIndex, int toIndex) {
    assert(fromIndex <= toIndex);
    assert(fromIndex >= 0);
    assert(toIndex <= size());
    std::vector<bool>::iterator first = bitset_.begin() + fromIndex;
    std::vector<bool>::iterator last  = bitset_.begin() + toIndex;

    OpenMDBitSet result;
    std::copy(first, last, std::back_inserter(result.bitset_));
    return result;
  }

  bool OpenMDBitSet::none() {
    std::vector<bool>::iterator i =
        std::find(bitset_.begin(), bitset_.end(), true);
    return i == bitset_.end() ? true : false;
  }

  int OpenMDBitSet::nextOffBit(int fromIndex) const {
    if (fromIndex <= -1) {
      // in case -1 or other negative number is passed to this function
      return -1;
    }

    ++fromIndex;
    while (fromIndex < size()) {
      if (!bitset_[fromIndex]) { return fromIndex; }
      ++fromIndex;
    }

    return -1;
  }

  int OpenMDBitSet::nthOffBit(unsigned long int fromIndex,
                              unsigned long int n) const {
    std::vector<int> indices;
    for (int i = nextOffBit(fromIndex); i != -1; i = nextOffBit(i + 1)) {
      indices.push_back(i);
    }

    if (n <= indices.size()) return indices[n];
    return -1;
  }

  int OpenMDBitSet::nextOnBit(int fromIndex) const {
    if (fromIndex <= -1) {
      // in case -1 or other negative number is passed to this function
      return -1;
    }

    ++fromIndex;
    while (fromIndex < size()) {
      if (bitset_[fromIndex]) { return fromIndex; }
      ++fromIndex;
    }

    return -1;
  }

  int OpenMDBitSet::nthOnBit(unsigned long int fromIndex,
                             unsigned long int n) const {
    std::vector<int> indices;
    for (int i = nextOnBit(fromIndex); i != -1; i = nextOnBit(i + 1)) {
      indices.push_back(i);
    }

    if (n <= indices.size()) return indices[n];
    return -1;
  }

  void OpenMDBitSet::andOperator(const OpenMDBitSet& bs) {
    assert(size() == bs.size());

    std::transform(bs.bitset_.begin(), bs.bitset_.end(), bitset_.begin(),
                   bitset_.begin(), std::logical_and<bool>());
  }

  void OpenMDBitSet::orOperator(const OpenMDBitSet& bs) {
    assert(size() == bs.size());
    std::transform(bs.bitset_.begin(), bs.bitset_.end(), bitset_.begin(),
                   bitset_.begin(), std::logical_or<bool>());
  }

  void OpenMDBitSet::xorOperator(const OpenMDBitSet& bs) {
    assert(size() == bs.size());
    std::transform(bs.bitset_.begin(), bs.bitset_.end(), bitset_.begin(),
                   bitset_.begin(), std::bit_xor<bool>());
  }

  void OpenMDBitSet::setBits(int fromIndex, int toIndex, bool value) {
    assert(fromIndex <= toIndex);
    assert(fromIndex >= 0);
    assert(toIndex <= size());
    std::vector<bool>::iterator first = bitset_.begin() + fromIndex;
    std::vector<bool>::iterator last  = bitset_.begin() + toIndex;
    std::fill(first, last, value);
  }

  void OpenMDBitSet::resize(int nbits) {
    int oldSize = size();
    bitset_.resize(nbits);
    if (nbits > oldSize) {
      std::fill(bitset_.begin() + oldSize, bitset_.end(), false);
    }
  }

  OpenMDBitSet operator|(const OpenMDBitSet& bs1, const OpenMDBitSet& bs2) {
    assert(bs1.size() == bs2.size());

    OpenMDBitSet result(bs1);
    result |= bs2;
    return result;
  }

  OpenMDBitSet operator&(const OpenMDBitSet& bs1, const OpenMDBitSet& bs2) {
    assert(bs1.size() == bs2.size());

    OpenMDBitSet result(bs1);
    result &= bs2;
    return result;
  }

  OpenMDBitSet operator^(const OpenMDBitSet& bs1, const OpenMDBitSet& bs2) {
    assert(bs1.size() == bs2.size());

    OpenMDBitSet result(bs1);
    result ^= bs2;
    return result;
  }

  OpenMDBitSet operator-(const OpenMDBitSet& bs1, const OpenMDBitSet& bs2) {
    assert(bs1.size() == bs2.size());

    OpenMDBitSet result(bs1);
    result -= bs2;
    return result;
  }

  bool operator==(const OpenMDBitSet& bs1, const OpenMDBitSet& bs2) {
    assert(bs1.size() == bs2.size());
    return std::equal(bs1.bitset_.begin(), bs1.bitset_.end(),
                      bs2.bitset_.begin());
  }

  OpenMDBitSet OpenMDBitSet::parallelReduce() {
    OpenMDBitSet result;

#ifdef IS_MPI

    // This is necessary because std::vector<bool> isn't really a
    // std::vector, so we can't pass the address of the first element
    // to the MPI call that follows.  We first have to convert to a
    // std::vector<int> to do the logical_or Allreduce call, then back
    // convert it into the vector<bool>.

    std::vector<int> bsInt(bitset_.begin(), bitset_.end());

    MPI_Allreduce(MPI_IN_PLACE, &bsInt[0], bsInt.size(), MPI_INT, MPI_LOR,
                  MPI_COMM_WORLD);

    std::transform(bsInt.begin(), bsInt.end(),
                   std::back_inserter(result.bitset_),
                   [](int val) { return val != 0; });
#else

    // Not in MPI?  Just return a copy of the current bitset:
    std::copy(bitset_.begin(), bitset_.end(),
              std::back_inserter(result.bitset_));
#endif

    return result;
  }

  // std::istream& operator>> ( std::istream& is, const OpenMDBitSet& bs) {
  //
  //    return is;
  //}

  std::ostream& operator<<(std::ostream& os, const OpenMDBitSet& bs) {
    for (unsigned int i = 0; i < bs.bitset_.size(); ++i) {
      std::string val = bs[i] ? "true" : "false";
      os << "OpenMDBitSet[" << i << "] = " << val << std::endl;
    }

    return os;
  }

}  // namespace OpenMD
