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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#ifndef SELECTION_SELECTIONSET_HPP
#define SELECTION_SELECTIONSET_HPP

#include <iostream>
#include <vector>
#include "utils/OpenMDBitSet.hpp"

namespace OpenMD {

  /** 
   * The SelectionType enum.
   *
   * This is used to sort different types of selections by object type
   */
  enum SelectionType {
    STUNTDOUBLE = 0,     /**< StuntDoubles (Atoms & RigidBodies) */
    BOND = 1,            /**< Bonds */
    BEND = 2,            /**< Bends */
    TORSION = 3,         /**< Torsions */
    INVERSION = 4,       /**< Inversions */
    N_SELECTIONTYPES = 5
  };

  class SelectionSet {
  public:
    /** */
    SelectionSet() { bitsets_.resize(N_SELECTIONTYPES); }
    /** */
    SelectionSet(std::vector<int> nbits);

    /** Returns the number of bits set to true in this SelectionSet.  */
    std::vector<int> countBits();

    /** Sets the bit at the specified index to to the complement of its current value. */
    void flip(std::vector<int> bitIndex);
 
    /** Sets each bit from the specified fromIndex(inclusive) to the specified toIndex(exclusive) to the complement of its current value. */
    void flip(std::vector<int> fromIndex, std::vector<int> toIndex); 

    /** Sets each bit to the complement of its current value. */
    void flip();
        
    /** Returns the value of the bit with the specified index. */
    std::vector<bool> get(std::vector<int> bitIndex);
        
    /** Returns a new SelectionSet composed of bits from this SelectionSet from fromIndex(inclusive) to toIndex(exclusive). */
    SelectionSet get(std::vector<int> fromIndex, std::vector<int> toIndex); 
        
    /** Returns true if any bits are set to true */
    std::vector<bool> any();

    /** Returns true if no bits are set to true */
    std::vector<bool> none();

    std::vector<int> firstOffBit() const;
        
    /** Returns the index of the first bit that is set to false that occurs on or after the specified starting index.*/
    std::vector<int> nextOffBit(std::vector<int> fromIndex) const; 

    std::vector<int> firstOnBit() const;
        
    /** Returns the index of the first bit that is set to true that occurs on or after the specified starting index. */
    std::vector<int> nextOnBit(std::vector<int> fromIndex) const; 
        
    /** Performs a logical AND of this target bit set with the argument bit set. */
    void andOperator (const SelectionSet& bs);
       
    /** Performs a logical OR of this bit set with the bit set argument. */
    void orOperator (const SelectionSet& bs); 
        
    /** Performs a logical XOR of this bit set with the bit set argument. */
    void xorOperator (const SelectionSet& bs);        
               
    void setBitOn(std::vector<int> bitIndex);

    void setBitOff(std::vector<int> bitIndex);

    //void setRangeOn(std::vector<int> fromIndex, std::vector<int> toIndex) {  setBits(fromIndex, toIndex, true);  }

    //void setRangeOff(std::vector<int> fromIndex, std::vector<int> toIndex) {  setBits(fromIndex, toIndex, false);  }        

    /** Sets all of the bits in this SelectionSet to false. */
    void clearAll();

    void setAll();
        
    /** Returns the number of bits of space actually in use by this SelectionSet to represent bit values. */
    std::vector<int> size() const;

    /** Changes the size of SelectionSet*/
    void resize(std::vector<int> nbits);
        
    SelectionSet& operator&= (const SelectionSet &ss) {  andOperator (ss); return *this; }
    SelectionSet& operator|= (const SelectionSet &ss) { orOperator (ss); return *this; }
    SelectionSet& operator^= (const SelectionSet &ss) { xorOperator (ss); return *this; }
    SelectionSet& operator-= (const SelectionSet &ss) { 
      SelectionSet tmp = *this ^ ss;
      *this &= tmp;
      return *this;
    }

    SelectionSet parallelReduce();
        
    std::vector<bool>  operator[] (std::vector<int> bitIndex)  const {  
      std::vector<bool> result(N_SELECTIONTYPES);
      for (int i = 0; i < N_SELECTIONTYPES; i++) 
        result[i] = bitsets_[i][bitIndex[i]];
      return result;
    }

    friend SelectionSet operator| (const SelectionSet& bs1, const SelectionSet& bs2);
    friend SelectionSet operator& (const SelectionSet& bs1, const SelectionSet& bs2);
    friend SelectionSet operator^ (const SelectionSet& bs1, const SelectionSet& bs2);
    friend SelectionSet operator- (const SelectionSet& bs1, const SelectionSet& bs2);
        
    friend bool operator== (const SelectionSet & bs1, const SelectionSet &bs2);

    //friend std::istream& operator>> ( std::istream&, const SelectionSet& bs);
    friend std::ostream& operator<< ( std::ostream&, const SelectionSet& bs) ;

    std::vector<OpenMDBitSet> bitsets_;

    

  private:

    /** Sets the bit at the specified index to the specified value. */
    //void setBit(std::vector<int> bitIndex, bool value);
        
    /** Sets the bits from the specified fromIndex(inclusive) to the specified toIndex(exclusive) to the specified value. */
    //void setBits(std::vector<int> fromIndex, std::vector<int> toIndex, bool value);
        

  };


}
#endif
