/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#ifndef UTILS_OPENMDBITSET_HPP
#define UTILS_OPENMDBITSET_HPP

#include <iostream>
#include <vector>
namespace OpenMD {

  /**
   * @class OpenMDBitSet OpenMDBitSet.hpp "OpenMDBitSet.hpp"
   * @brief OpenMDBitSet is a wrapper class of std::vector<bool> to act as a growable std::bitset 
   */
  class OpenMDBitSet {
  public:
    /** */
    OpenMDBitSet() {}
    /** */
    OpenMDBitSet(int nbits) : bitset_(nbits) {clearAll(); }

    /** Returns the number of bits set to true in this OpenMDBitSet.  */
    int countBits();

    /** Sets the bit at the specified index to to the complement of its current value. */
    void flip(int bitIndex) {  bitset_[bitIndex] = !bitset_[bitIndex];  }
 
    /** Sets each bit from the specified fromIndex(inclusive) to the specified toIndex(exclusive) to the complement of its current value. */
    void flip(int fromIndex, int toIndex); 

    /** Sets each bit to the complement of its current value. */
    void flip() { flip(0, size()); }
        
    /** Returns the value of the bit with the specified index. */
    bool get(int bitIndex) {  return bitset_[bitIndex];  }
        
    /** Returns a new OpenMDBitSet composed of bits from this OpenMDBitSet from fromIndex(inclusive) to toIndex(exclusive). */
    OpenMDBitSet get(int fromIndex, int toIndex); 
        
    /** Returns true if any bits are set to true */
    bool any() {return !none(); }

    /** Returns true if no bits are set to true */
    bool none();

    int firstOffBit() const { return !bitset_[0] ? 0 : nextOffBit(0); }
        
    /** Returns the index of the first bit that is set to false that occurs on or after the specified starting index.*/
    int nextOffBit(int fromIndex) const; 

    int firstOnBit() const { return bitset_[0] ? 0 : nextOnBit(0); }
        
    /** Returns the index of the first bit that is set to true that occurs on or after the specified starting index. */
    int nextOnBit(int fromIndex) const; 
        
    /** Performs a logical AND of this target bit set with the argument bit set. */
    void andOperator (const OpenMDBitSet& bs);
       
    /** Performs a logical OR of this bit set with the bit set argument. */
    void orOperator (const OpenMDBitSet& bs); 
        
    /** Performs a logical XOR of this bit set with the bit set argument. */
    void xorOperator (const OpenMDBitSet& bs);        
               
    void setBitOn(int bitIndex) {  setBit(bitIndex, true);  }

    void setBitOff(int bitIndex) {  setBit(bitIndex, false);  }

    void setRangeOn(int fromIndex, int toIndex) {  setBits(fromIndex, toIndex, true);  }

    void setRangeOff(int fromIndex, int toIndex) {  setBits(fromIndex, toIndex, false);  }        

    /** Sets all of the bits in this OpenMDBitSet to false. */
    void clearAll() {  setRangeOff(0, size());  }         

    void setAll() {  setRangeOn(0, size());  }        
        
    /** Returns the number of bits of space actually in use by this OpenMDBitSet to represent bit values. */
    int size() const {  return bitset_.size();  }

    /** Changes the size of OpenMDBitSet*/
    void resize(int nbits);
        
    OpenMDBitSet& operator&= (const OpenMDBitSet &bs) {  andOperator (bs); return *this; }
    OpenMDBitSet& operator|= (const OpenMDBitSet &bs) { orOperator (bs); return *this; }
    OpenMDBitSet& operator^= (const OpenMDBitSet &bs) { xorOperator (bs); return *this; }
    OpenMDBitSet& operator-= (const OpenMDBitSet &bs) { 
      OpenMDBitSet tmp = *this ^ bs;
      *this &= tmp;
      return *this;
    }
        
    bool operator[] (int bitIndex)  const {  return bitset_[bitIndex];  }
    friend OpenMDBitSet operator| (const OpenMDBitSet& bs1, const OpenMDBitSet& bs2);
    friend OpenMDBitSet operator& (const OpenMDBitSet& bs1, const OpenMDBitSet& bs2);
    friend OpenMDBitSet operator^ (const OpenMDBitSet& bs1, const OpenMDBitSet& bs2);
    friend OpenMDBitSet operator- (const OpenMDBitSet& bs1, const OpenMDBitSet& bs2);
        
    friend bool operator== (const OpenMDBitSet & bs1, const OpenMDBitSet &bs2);

    //friend std::istream& operator>> ( std::istream&, const OpenMDBitSet& bs);
    friend std::ostream& operator<< ( std::ostream&, const OpenMDBitSet& bs) ;

  private:

    /** Sets the bit at the specified index to the specified value. */
    void setBit(int bitIndex, bool value) { bitset_[bitIndex] = value; }
        
    /** Sets the bits from the specified fromIndex(inclusive) to the specified toIndex(exclusive) to the specified value. */
    void setBits(int fromIndex, int toIndex, bool value);
        
    std::vector<bool> bitset_;
  };


}
#endif
