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

#ifndef UTILS_OPENMDBITSET_HPP
#define UTILS_OPENMDBITSET_HPP

#include <iostream>
#include <vector>

namespace OpenMD {

  /**
   * @class OpenMDBitSet OpenMDBitSet.hpp "OpenMDBitSet.hpp"
   * @brief OpenMDBitSet is a wrapper class of std::vector<bool> to act as a
   * growable std::bitset
   */
  class OpenMDBitSet {
  public:
    /** */
    OpenMDBitSet() {}
    /** */
    OpenMDBitSet(size_t nbits) : bitset_(nbits) { clearAll(); }

    /** Returns the number of bits set to true in this OpenMDBitSet.  */
    int countBits();

    /** Sets the bit at the specified index to to the complement of its current
     * value. */
    void flip(size_t bitIndex) { bitset_[bitIndex] = !bitset_[bitIndex]; }

    /** Sets each bit from the specified fromIndex(inclusive) to the specified
     * toIndex(exclusive) to the complement of its current value. */
    void flip(size_t fromIndex, size_t toIndex);

    /** Sets each bit to the complement of its current value. */
    void flip() { flip(0, size()); }

    /** Returns the value of the bit with the specified index. */
    bool get(size_t bitIndex) { return bitset_[bitIndex]; }

    /** Returns a new OpenMDBitSet composed of bits from this OpenMDBitSet from
     * fromIndex(inclusive) to toIndex(exclusive). */
    OpenMDBitSet get(size_t fromIndex, size_t toIndex);

    /** Returns true if any bits are set to true */
    bool any() { return !none(); }

    /** Returns true if no bits are set to true */
    bool none();

    int firstOffBit() const { return !bitset_[0] ? 0 : nextOffBit(0); }

    /** Returns the index of the first bit that is set to false that occurs on
     * or after the specified starting index.*/
    int nextOffBit(int fromIndex) const;

    /** Returns the index of the n^th bit that is set to false. */
    int nthOffBit(unsigned long int n) const;

    int firstOnBit() const { return bitset_[0] ? 0 : nextOnBit(0); }

    /** Returns the index of the first bit that is set to true that occurs on or
     * after the specified starting index. */
    int nextOnBit(int fromIndex) const;

    /** Returns the index of the n^th bit that is set to true. */
    int nthOnBit(unsigned long int n) const;

    /** Performs a logical AND of this target bit set with the argument bit set.
     */
    void andOperator(const OpenMDBitSet& bs);

    /** Performs a logical OR of this bit set with the bit set argument. */
    void orOperator(const OpenMDBitSet& bs);

    /** Performs a logical XOR of this bit set with the bit set argument. */
    void xorOperator(const OpenMDBitSet& bs);

    void setBitOn(size_t bitIndex) { setBit(bitIndex, true); }

    void setBitOff(size_t bitIndex) { setBit(bitIndex, false); }

    void setRangeOn(size_t fromIndex, size_t toIndex) {
      setBits(fromIndex, toIndex, true);
    }

    void setRangeOff(size_t fromIndex, size_t toIndex) {
      setBits(fromIndex, toIndex, false);
    }

    /** Sets all of the bits in this OpenMDBitSet to false. */
    void clearAll() { setRangeOff(0, size()); }

    void setAll() { setRangeOn(0, size()); }

    /** Returns the number of bits of space actually in use by this OpenMDBitSet
     * to represent bit values. */
    size_t size() const { return bitset_.size(); }

    /** Changes the size of OpenMDBitSet*/
    void resize(size_t nbits);

    OpenMDBitSet& operator&=(const OpenMDBitSet& bs) {
      andOperator(bs);
      return *this;
    }
    OpenMDBitSet& operator|=(const OpenMDBitSet& bs) {
      orOperator(bs);
      return *this;
    }
    OpenMDBitSet& operator^=(const OpenMDBitSet& bs) {
      xorOperator(bs);
      return *this;
    }
    OpenMDBitSet& operator-=(const OpenMDBitSet& bs) {
      OpenMDBitSet tmp = *this ^ bs;
      *this &= tmp;
      return *this;
    }

    OpenMDBitSet parallelReduce();

    bool operator[](int bitIndex) const { return bitset_[bitIndex]; }
    friend OpenMDBitSet operator|(const OpenMDBitSet& bs1,
                                  const OpenMDBitSet& bs2);
    friend OpenMDBitSet operator&(const OpenMDBitSet& bs1,
                                  const OpenMDBitSet& bs2);
    friend OpenMDBitSet operator^(const OpenMDBitSet& bs1,
                                  const OpenMDBitSet& bs2);
    friend OpenMDBitSet operator-(const OpenMDBitSet& bs1,
                                  const OpenMDBitSet& bs2);

    friend bool operator==(const OpenMDBitSet& bs1, const OpenMDBitSet& bs2);

    // friend std::istream& operator>> ( std::istream&, const OpenMDBitSet& bs);
    friend std::ostream& operator<<(std::ostream&, const OpenMDBitSet& bs);

  private:
    /** Sets the bit at the specified index to the specified value. */
    void setBit(size_t bitIndex, bool value) { bitset_[bitIndex] = value; }

    /** Sets the bits from the specified fromIndex(inclusive) to the specified
     * toIndex(exclusive) to the specified value. */
    void setBits(size_t fromIndex, size_t toIndex, bool value);

    std::vector<bool> bitset_;
  };
}  // namespace OpenMD

#endif
