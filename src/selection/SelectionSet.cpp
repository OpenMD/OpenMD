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

#include "selection/SelectionSet.hpp"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <string>

#ifdef IS_MPI
#include <mpi.h>
#endif

namespace OpenMD {

  SelectionSet::SelectionSet(std::vector<int> nbits) {
    bitsets_.resize(N_SELECTIONTYPES);
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      bitsets_[i] = OpenMDBitSet(nbits[i]);
    clearAll();
  }

  std::vector<int> SelectionSet::countBits() {
    std::vector<int> result(N_SELECTIONTYPES, 0);
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      result[i] = bitsets_[i].countBits();
    return result;
  }

  void SelectionSet::flip(std::vector<int> bitIndex) {
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      bitsets_[i].flip(bitIndex[i]);
  }

  void SelectionSet::flip(std::vector<int> fromIndex,
                          std::vector<int> toIndex) {
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      bitsets_[i].flip(fromIndex[i], toIndex[i]);
  }

  void SelectionSet::flip() {
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      bitsets_[i].flip();
  }

  std::vector<bool> SelectionSet::get(std::vector<int> bitIndex) {
    std::vector<bool> result(N_SELECTIONTYPES);
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      result[i] = bitsets_[i].get(bitIndex[i]);
    return result;
  }

  SelectionSet SelectionSet::get(std::vector<int> fromIndex,
                                 std::vector<int> toIndex) {
    SelectionSet result;

    for (int i = 0; i < N_SELECTIONTYPES; i++)
      result.bitsets_[i] = bitsets_[i].get(fromIndex[i], toIndex[i]);

    return result;
  }

  std::vector<bool> SelectionSet::any() {
    std::vector<bool> result(N_SELECTIONTYPES);
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      result[i] = bitsets_[i].any();
    return result;
  }

  std::vector<bool> SelectionSet::none() {
    std::vector<bool> result(N_SELECTIONTYPES);
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      result[i] = bitsets_[i].none();
    return result;
  }

  std::vector<int> SelectionSet::firstOffBit() const {
    std::vector<int> result(N_SELECTIONTYPES);
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      result[i] = bitsets_[i].firstOffBit();
    return result;
  }

  std::vector<int> SelectionSet::nextOffBit(std::vector<int> fromIndex) const {
    std::vector<int> result(N_SELECTIONTYPES);
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      result[i] = bitsets_[i].nextOffBit(fromIndex[i]);
    return result;
  }

  std::vector<int> SelectionSet::firstOnBit() const {
    std::vector<int> result(N_SELECTIONTYPES);
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      result[i] = bitsets_[i].firstOnBit();
    return result;
  }

  std::vector<int> SelectionSet::nextOnBit(std::vector<int> fromIndex) const {
    std::vector<int> result(N_SELECTIONTYPES);
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      result[i] = bitsets_[i].nextOnBit(fromIndex[i]);
    return result;
  }

  void SelectionSet::andOperator(const SelectionSet& ss) {
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      bitsets_[i] &= ss.bitsets_[i];
  }

  void SelectionSet::orOperator(const SelectionSet& ss) {
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      bitsets_[i] |= ss.bitsets_[i];
  }

  void SelectionSet::xorOperator(const SelectionSet& ss) {
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      bitsets_[i] ^= ss.bitsets_[i];
  }

  // void SelectionSet::setBits(std::vector<int> fromIndex,
  //                           std::vector<int> toIndex, bool value) {
  //  for (int i = 0; i < N_SELECTIONTYPES; i++)
  //    bitsets_[i].setBits(fromIndex[i], toIndex[i], value);
  //}

  void SelectionSet::clearAll() {
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      bitsets_[i].clearAll();
  }

  void SelectionSet::setAll() {
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      bitsets_[i].setAll();
  }

  std::vector<int> SelectionSet::size() const {
    std::vector<int> result(N_SELECTIONTYPES);
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      result[i] = bitsets_[i].size();
    return result;
  }

  void SelectionSet::resize(std::vector<int> nbits) {
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      bitsets_[i].resize(nbits[i]);
  }

  SelectionSet operator|(const SelectionSet& ss1, const SelectionSet& ss2) {
    SelectionSet result(ss1);
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      result.bitsets_[i] |= ss2.bitsets_[i];
    return result;
  }

  SelectionSet operator&(const SelectionSet& ss1, const SelectionSet& ss2) {
    SelectionSet result(ss1);
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      result.bitsets_[i] &= ss2.bitsets_[i];
    return result;
  }

  SelectionSet operator^(const SelectionSet& ss1, const SelectionSet& ss2) {
    SelectionSet result(ss1);
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      result.bitsets_[i] ^= ss2.bitsets_[i];
    return result;
  }

  SelectionSet operator-(const SelectionSet& ss1, const SelectionSet& ss2) {
    SelectionSet result(ss1);
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      result.bitsets_[i] -= ss2.bitsets_[i];
    return result;
  }

  bool operator==(const SelectionSet& ss1, const SelectionSet& ss2) {
    for (int i = 0; i < N_SELECTIONTYPES; i++) {
      assert(ss1.bitsets_[i].size() == ss2.bitsets_[i].size());
      if (!(ss1.bitsets_[i] == ss2.bitsets_[i])) return false;
    }

    return true;
  }

  SelectionSet SelectionSet::parallelReduce() {
    SelectionSet result;
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      result.bitsets_[i] = bitsets_[i].parallelReduce();

    return result;
  }

  // std::istream& operator>> ( std::istream& is, const OpenMDBitSet& bs) {
  //
  //    return is;
  //}

  std::ostream& operator<<(std::ostream& os, const SelectionSet& ss) {
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      os << "SelectionSet[" << i << "] = " << ss.bitsets_[i] << std::endl;
    return os;
  }

  // void SelectionSet::setBit(std::vector<int> bitIndex, bool value) {
  //  for (int i = 0; i < N_SELECTIONTYPES; i++)
  //    bitsets_[i].setBit(bitIndex[i], value);
  //}

}  // namespace OpenMD
