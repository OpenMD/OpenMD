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

#ifdef IS_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <cassert>
#include <string>
#include <iterator>

#include "selection/SelectionSet.hpp"
#include "utils/Algorithm.hpp"

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

  void SelectionSet::andOperator (const SelectionSet& ss) {
    for (int i = 0; i < N_SELECTIONTYPES; i++) 
      bitsets_[i] &= ss.bitsets_[i];
  }

  void SelectionSet::orOperator (const SelectionSet& ss) {
    for (int i = 0; i < N_SELECTIONTYPES; i++) 
      bitsets_[i] |= ss.bitsets_[i];
  }

  void SelectionSet::xorOperator (const SelectionSet& ss) {
    for (int i = 0; i < N_SELECTIONTYPES; i++) 
      bitsets_[i] ^= ss.bitsets_[i];
  }
   
  //void SelectionSet::setBits(std::vector<int> fromIndex, 
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

  SelectionSet operator| (const SelectionSet& ss1, const SelectionSet& ss2) {
    SelectionSet result(ss1);
    for (int i = 0; i < N_SELECTIONTYPES; i++)           
      result.bitsets_[i] |= ss2.bitsets_[i];
    return result;
  }

  SelectionSet operator& (const SelectionSet& ss1, const SelectionSet& ss2) {
    SelectionSet result(ss1);
    for (int i = 0; i < N_SELECTIONTYPES; i++)           
      result.bitsets_[i] &= ss2.bitsets_[i];
    return result;
  }

  SelectionSet operator^ (const SelectionSet& ss1, const SelectionSet& ss2) {
    SelectionSet result(ss1);
    for (int i = 0; i < N_SELECTIONTYPES; i++)           
      result.bitsets_[i] ^= ss2.bitsets_[i];
    return result;
  }

  SelectionSet operator- (const SelectionSet& ss1, const SelectionSet& ss2) {
    SelectionSet result(ss1);
    for (int i = 0; i < N_SELECTIONTYPES; i++)           
      result.bitsets_[i] -= ss2.bitsets_[i];
    return result;
  }

  bool operator== (const SelectionSet & ss1, const SelectionSet &ss2) {

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



  //std::istream& operator>> ( std::istream& is, const OpenMDBitSet& bs) {
  //
  //    return is;
  //}

  std::ostream& operator<< ( std::ostream& os, const SelectionSet& ss) {
    for (int i = 0; i < N_SELECTIONTYPES; i++)
      os << "SelectionSet[" << i << "] = " << ss.bitsets_[i] << std::endl;
    return os;
  }

  //void SelectionSet::setBit(std::vector<int> bitIndex, bool value) { 
  //  for (int i = 0; i < N_SELECTIONTYPES; i++)      
  //    bitsets_[i].setBit(bitIndex[i], value);
  //}
          
}
