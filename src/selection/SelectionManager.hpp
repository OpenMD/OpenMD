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

#ifndef SELECTION_SELECTIONMANAGER_HPP
#define SELECTION_SELECTIONMANAGER_HPP

#include "selection/SelectionSet.hpp"
#include "primitives/StuntDouble.hpp"
#include "primitives/Bond.hpp"
#include "primitives/Bend.hpp"
#include "primitives/Torsion.hpp"
#include "primitives/Inversion.hpp"

namespace OpenMD {

  class SimInfo;
  class SelectionManager {
  public:
    SelectionManager(SimInfo* info);

    void addSelection(StuntDouble* sd) {
      
      ss_.bitsets_[STUNTDOUBLE].setBitOn(sd->getGlobalIndex());
    }
    void addSelection(Bond* b) {
      ss_.bitsets_[BOND].setBitOn(b->getGlobalIndex());
    }
    void addSelection(Bend* b) {
      ss_.bitsets_[BEND].setBitOn(b->getGlobalIndex());
    }
    void addSelection(Torsion* t) {
      ss_.bitsets_[TORSION].setBitOn(t->getGlobalIndex());
    }
    void addSelection(Inversion* i) {
      ss_.bitsets_[INVERSION].setBitOn(i->getGlobalIndex());
    }
        
    void addSelectionSet(const SelectionSet& bs) {
      ss_.bitsets_[STUNTDOUBLE] |= bs.bitsets_[STUNTDOUBLE];
    }
    void addBondSelectionSet(const SelectionSet& bs) {
      ss_.bitsets_[BOND] |= bs.bitsets_[BOND];
    }
    void addBendSelectionSet(const SelectionSet& bs) {
      ss_.bitsets_[BEND] |= bs.bitsets_[BEND];
    }
    void addTorsionSelectionSet(const SelectionSet& bs) {
      ss_.bitsets_[TORSION] |= bs.bitsets_[TORSION];
    }
    void addInversionSelectionSet(const SelectionSet& bs) {
      ss_.bitsets_[INVERSION] |= bs.bitsets_[INVERSION];
    }

    bool isEmpty() {
      return ss_.bitsets_[STUNTDOUBLE].none() && ss_.bitsets_[BOND].none() 
        && ss_.bitsets_[BEND].none()  && ss_.bitsets_[TORSION].none() 
        && ss_.bitsets_[INVERSION].none();
    }

    void setSelectionSet(const SelectionSet& bs) {
      for (int i = 0; i < N_SELECTIONTYPES; i++) 
        ss_.bitsets_[i] = bs.bitsets_[i];
    }

    //void setSelectionSet(const SelectionSet& bs) {
    //  ss_.bitsets_[STUNTDOUBLE] = bs.bitsets_[];           
    //}
    void setBondSelectionSet(const SelectionSet& bs) {
      ss_.bitsets_[BOND] = bs.bitsets_[BOND];           
    }
    void setBendSelectionSet(const SelectionSet& bs) {
      ss_.bitsets_[BEND] = bs.bitsets_[BEND];           
    }
    void setTorsionSelectionSet(const SelectionSet& bs) {
      ss_.bitsets_[TORSION] = bs.bitsets_[TORSION];           
    }
    void setInversionSelectionSet(const SelectionSet& bs) {
      ss_.bitsets_[INVERSION] = bs.bitsets_[INVERSION];           
    }

    std::vector<int> getSelectionCounts() {
      std::vector<int> counts(N_SELECTIONTYPES,0);
      for (int i = 0; i < N_SELECTIONTYPES; i++) {
        counts[i] = ss_.bitsets_[i].countBits();
      }
      return counts;
    }
     
    int getSelectionCount() {
      return ss_.bitsets_[STUNTDOUBLE].countBits();
    }
    int getBondSelectionCount() {
      return ss_.bitsets_[BOND].countBits();
    }
    int getBendSelectionCount() {
      return ss_.bitsets_[BEND].countBits();
    }
    int getTorsionSelectionCount() {
      return ss_.bitsets_[TORSION].countBits();
    }
    int getInversionSelectionCount() {
      return ss_.bitsets_[INVERSION].countBits();
    }
    
    SelectionSet getSelectionSet() {
      return ss_;
    }
    /*    SelectionSet getBondSelectionSet() {
      return ss_.bitsets_[BOND];
    }
    SelectionSet getBendSelectionSet() {
      return ss_.bitsets_[BEND];
    }
    SelectionSet getTorsionSelectionSet() {
      return ss_.bitsets_[TORSION];
    }
    SelectionSet getInversionSelectionSet() {
      return ss_.bitsets_[INVERSION];
    }
    */

    void setSelection(StuntDouble* sd) {
      ss_.bitsets_[STUNTDOUBLE].clearAll();
      ss_.bitsets_[STUNTDOUBLE].setBitOn(sd->getGlobalIndex());
    }
    void setSelection(Bond* b) {
      ss_.bitsets_[BOND].clearAll();
      ss_.bitsets_[BOND].setBitOn(b->getGlobalIndex());
    }
    void setSelection(Bend* b) {
      ss_.bitsets_[BEND].clearAll();
      ss_.bitsets_[BEND].setBitOn(b->getGlobalIndex());
    }
    void setSelection(Torsion* t) {
      ss_.bitsets_[TORSION].clearAll();
      ss_.bitsets_[TORSION].setBitOn(t->getGlobalIndex());
    }
    void setSelection(Inversion* i) {
      ss_.bitsets_[INVERSION].clearAll();
      ss_.bitsets_[INVERSION].setBitOn(i->getGlobalIndex());
    }

    void toggleSelection(StuntDouble* sd) {
      ss_.bitsets_[STUNTDOUBLE].flip(sd->getGlobalIndex());
    }
    void toggleSelection(Bond* b) {
      ss_.bitsets_[BOND].flip(b->getGlobalIndex());
    }
    void toggleSelection(Bend* b) {
      ss_.bitsets_[BEND].flip(b->getGlobalIndex());
    }
    void toggleSelection(Torsion* t) {
      ss_.bitsets_[TORSION].flip(t->getGlobalIndex());
    }
    void toggleSelection(Inversion* i) {
      ss_.bitsets_[INVERSION].flip(i->getGlobalIndex());
    }

    void toggleSelection() {
      for (int i = 0; i < N_SELECTIONTYPES; i++) 
        ss_.bitsets_[i].flip();
    }
        
    void selectAll() {
      for (int i = 0; i < N_SELECTIONTYPES; i++) 
        ss_.bitsets_[i].setAll();
    }

    void clearSelection() {
      for (int i = 0; i < N_SELECTIONTYPES; i++) 
        ss_.bitsets_[i].clearAll();
    }

    void clearSelection(StuntDouble* sd) {
      ss_.bitsets_[STUNTDOUBLE].setBitOff(sd->getGlobalIndex());
    }
    void clearSelection(Bond* b) {
      ss_.bitsets_[BOND].setBitOff(b->getGlobalIndex());
    }
    void clearSelection(Bend* b) {
      ss_.bitsets_[BEND].setBitOff(b->getGlobalIndex());
    }
    void clearSelection(Torsion* t) {
      ss_.bitsets_[TORSION].setBitOff(t->getGlobalIndex());
    }
    void clearSelection(Inversion* i) {
      ss_.bitsets_[INVERSION].setBitOff(i->getGlobalIndex());
    }

    bool isSelected(StuntDouble* sd) {
      return ss_.bitsets_[STUNTDOUBLE][sd->getGlobalIndex()];
    }
    bool isSelected(Bond* b) {
      return ss_.bitsets_[BOND][b->getGlobalIndex()];
    }
    bool isSelected(Bend* b) {
      return ss_.bitsets_[BEND][b->getGlobalIndex()];
    }
    bool isSelected(Torsion* t) {
      return ss_.bitsets_[TORSION][t->getGlobalIndex()];
    }
    bool isSelected(Inversion* i) {
      return ss_.bitsets_[INVERSION][i->getGlobalIndex()];
    }

    StuntDouble* beginSelected(int& i);
    StuntDouble* nextSelected(int& i);
    StuntDouble* beginUnselected(int& i);
    StuntDouble* nextUnSelected(int& i);

    Bond* beginSelectedBond(int& i);
    Bond* nextSelectedBond(int& i);
    Bond* beginUnselectedBond(int& i);
    Bond* nextUnSelectedBond(int& i);

    Bend* beginSelectedBend(int& i);
    Bend* nextSelectedBend(int& i);
    Bend* beginUnselectedBend(int& i);
    Bend* nextUnSelectedBend(int& i);

    Torsion* beginSelectedTorsion(int& i);
    Torsion* nextSelectedTorsion(int& i);
    Torsion* beginUnselectedTorsion(int& i);
    Torsion* nextUnSelectedTorsion(int& i);

    Inversion* beginSelectedInversion(int& i);
    Inversion* nextSelectedInversion(int& i);
    Inversion* beginUnselectedInversion(int& i);
    Inversion* nextUnSelectedInversion(int& i);

    SelectionManager& operator&= (const SelectionManager &sman) {
      for (int i = 0; i < N_SELECTIONTYPES; i++) 
        ss_.bitsets_[i] &= sman.ss_.bitsets_[i];
      return *this; 
    }
        
    SelectionManager& operator|= (const SelectionManager &sman) {
      for (int i = 0; i < N_SELECTIONTYPES; i++) 
        ss_.bitsets_[i] |= sman.ss_.bitsets_[i];
      return *this; 
    }
        
    SelectionManager& operator^= (const SelectionManager &sman) {
      for (int i = 0; i < N_SELECTIONTYPES; i++) 
        ss_.bitsets_[i] ^= sman.ss_.bitsets_[i];
      return *this; 
    }

    SelectionManager& operator-= (const SelectionManager &sman) {
      for (int i = 0; i < N_SELECTIONTYPES; i++) 
        ss_.bitsets_[i] -= sman.ss_.bitsets_[i];
      return *this; 
    }
        
    friend SelectionManager operator| (const SelectionManager& sman1, const SelectionManager& sman2);
    friend SelectionManager operator& (const SelectionManager& sman1, const SelectionManager& sman2);
    friend SelectionManager operator^ (const SelectionManager& sman1, const SelectionManager& sman2);
    friend SelectionManager operator-(const SelectionManager& sman1, const SelectionManager& sman2);
        
  private:
    SimInfo* info_;
    SelectionSet ss_;
    std::vector<int> nObjects_;
    std::vector<StuntDouble*> stuntdoubles_;
    std::vector<Bond*> bonds_;
    std::vector<Bend*> bends_;
    std::vector<Torsion*> torsions_;
    std::vector<Inversion*> inversions_;
  };

}
#endif
