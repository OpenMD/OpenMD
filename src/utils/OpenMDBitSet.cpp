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

  void OpenMDBitSet::flip(size_t fromIndex, size_t toIndex) {
    assert(fromIndex <= toIndex);
    assert(fromIndex >= 0);
    assert(toIndex <= size());
    std::vector<bool>::iterator first = bitset_.begin() + fromIndex;
    std::vector<bool>::iterator last  = bitset_.begin() + toIndex;

    std::transform(first, last, first, std::logical_not<bool>());
  }

  OpenMDBitSet OpenMDBitSet::get(size_t fromIndex, size_t toIndex) {
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
    while (fromIndex < static_cast<int>(size())) {
      if (!bitset_[fromIndex]) { return fromIndex; }
      ++fromIndex;
    }

    return -1;
  }

  int OpenMDBitSet::nthOffBit(unsigned long int n) const {
    std::vector<int> indices;
    for (int i = firstOffBit(); i != -1; i = nextOffBit(i)) {
      indices.push_back(i);
    }

    if (n < indices.size()) return indices[n];
    return -1;
  }

  int OpenMDBitSet::nextOnBit(int fromIndex) const {
    if (fromIndex <= -1) {
      // in case -1 or other negative number is passed to this function
      return -1;
    }

    ++fromIndex;
    while (fromIndex < static_cast<int>(size())) {
      if (bitset_[fromIndex]) { return fromIndex; }
      ++fromIndex;
    }

    return -1;
  }

  int OpenMDBitSet::nthOnBit(unsigned long int n) const {
    std::vector<int> indices;
    for (int i = firstOnBit(); i != -1; i = nextOnBit(i)) {
      indices.push_back(i);
    }

    if (n < indices.size()) return indices[n];
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

  void OpenMDBitSet::setBits(size_t fromIndex, size_t toIndex, bool value) {
    assert(fromIndex <= toIndex);
    assert(fromIndex >= 0);
    assert(toIndex <= size());
    std::vector<bool>::iterator first = bitset_.begin() + fromIndex;
    std::vector<bool>::iterator last  = bitset_.begin() + toIndex;
    std::fill(first, last, value);
  }

  void OpenMDBitSet::resize(size_t nbits) {
    size_t oldSize = size();
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
    for (size_t i = 0; i < bs.bitset_.size(); ++i) {
      std::string val = bs[i] ? "true" : "false";
      os << "OpenMDBitSet[" << i << "] = " << val << std::endl;
    }

    return os;
  }

}  // namespace OpenMD
