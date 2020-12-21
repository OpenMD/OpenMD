/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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

#ifndef SELECTION_SELECTIONMANAGER_HPP
#define SELECTION_SELECTIONMANAGER_HPP

#include "selection/SelectionSet.hpp"
#include "primitives/StuntDouble.hpp"
#include "primitives/Bond.hpp"
#include "primitives/Bend.hpp"
#include "primitives/Torsion.hpp"
#include "primitives/Inversion.hpp"
#include "primitives/Molecule.hpp"

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
    void addSelection(Molecule* m) {
      ss_.bitsets_[MOLECULE].setBitOn(m->getGlobalIndex());
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
    void addMoleculeSelectionSet(const SelectionSet& bs) {
      ss_.bitsets_[MOLECULE] |= bs.bitsets_[MOLECULE];
    }

    bool isEmpty() {
      return ss_.bitsets_[STUNTDOUBLE].none() && ss_.bitsets_[BOND].none() 
        && ss_.bitsets_[BEND].none()  && ss_.bitsets_[TORSION].none() 
        && ss_.bitsets_[INVERSION].none() && ss_.bitsets_[MOLECULE].none();
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
    void setMoleculeSelectionSet(const SelectionSet& bs) {
      ss_.bitsets_[MOLECULE] = bs.bitsets_[MOLECULE];           
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
    int getMoleculeSelectionCount() {
      return ss_.bitsets_[MOLECULE].countBits();
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
    void setSelection(Molecule* m) {
      ss_.bitsets_[MOLECULE].clearAll();
      ss_.bitsets_[MOLECULE].setBitOn(m->getGlobalIndex());
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
    void toggleSelection(Molecule* m) {
      ss_.bitsets_[MOLECULE].flip(m->getGlobalIndex());
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
    void clearSelection(Molecule* m) {
      ss_.bitsets_[MOLECULE].setBitOff(m->getGlobalIndex());
    }


    bool isGlobalIDSelected(int globalIndex) {
      return ss_.bitsets_[STUNTDOUBLE][globalIndex];
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
    bool isSelected(Molecule* m) {
      return ss_.bitsets_[MOLECULE][m->getGlobalIndex()];
    }

    /**
     * Finds the first selected StuntDouble in the selection.  In
     * parallel, this is the first selected StuntDouble that is the
     * responsibility of the local processor, not the first
     * StuntDouble in the global selection.
     * @return a pointer to the StuntDouble object, returns NULL if no StuntDouble was found.
     * @param i iterator used to keep track of the selection
     */
    StuntDouble* beginSelected(int& i);
    /**
     * Finds the next selected StuntDouble in the selection.  In
     * parallel, this is the next selected StuntDouble that is the
     * responsibility of the local processor, not the next
     * StuntDouble in the global selection.
     * @return a pointer to the StuntDouble object, returns NULL if no StuntDouble was found.
     * @param i iterator used to keep track of the selection
     */
    StuntDouble* nextSelected(int& i);
    /**
     * Finds the first unselected StuntDouble.  In
     * parallel, this is the first unselected StuntDouble that is the
     * responsibility of the local processor, not the first
     * StuntDouble in the global unselected pool.
     * @return a pointer to the StuntDouble object, returns NULL if no StuntDouble was found.
     * @param i iterator used to keep track of the selection
     */
    StuntDouble* beginUnselected(int& i);
    /**
     * Finds the next unselected StuntDouble.  In
     * parallel, this is the next unselected StuntDouble that is the
     * responsibility of the local processor, not the next
     * StuntDouble in the global unselected pool.
     * @return a pointer to the StuntDouble object, returns NULL if no StuntDouble was found.
     * @param i iterator used to keep track of the selection
     */
    StuntDouble* nextUnselected(int& i);
    /**
     * Finds the first selected Bond in the selection.  In parallel,
     * this is the first selected Bond that is the responsibility of
     * the local processor, not the first Bond in the global
     * selection.
     * @return a pointer to the Bond object, returns NULL if no Bond was found.
     * @param i iterator used to keep track of the selection
     */
    Bond* beginSelectedBond(int& i);
    /**
     * Finds the next selected Bond in the selection.  In parallel,
     * this is the next selected Bond that is the responsibility of
     * the local processor, not the next Bond in the global selection.
     * @return a pointer to the Bond object, returns NULL if no Bond was found.
     * @param i iterator used to keep track of the selection
     */
    Bond* nextSelectedBond(int& i);
    /**
     * Finds the first unselected Bond.  In parallel, this is the
     * first unselected Bond that is the responsibility of the local
     * processor, not the first Bond in the global unselected pool.
     * @return a pointer to the Bond object, returns NULL if no Bond was found.
     * @param i iterator used to keep track of the selection
     */
    Bond* beginUnselectedBond(int& i);
    /**
     * Finds the next unselected Bond.  In parallel, this is the
     * next unselected Bond that is the responsibility of the local
     * processor, not the next Bond in the global unselected pool.
     * @return a pointer to the Bond object, returns NULL if no Bond was found.
     * @param i iterator used to keep track of the selection
     */
    Bond* nextUnselectedBond(int& i);

    /**
     * Finds the first selected Bend in the selection.  In parallel,
     * this is the first selected Bend that is the responsibility of
     * the local processor, not the first Bend in the global
     * selection.
     * @return a pointer to the Bend object, returns NULL if no Bend was found.
     * @param i iterator used to keep track of the selection
     */
    Bend* beginSelectedBend(int& i);
    /**
     * Finds the next selected Bend in the selection.  In parallel,
     * this is the next selected Bend that is the responsibility of
     * the local processor, not the next Bend in the global selection.
     * @return a pointer to the Bend object, returns NULL if no Bend was found.
     * @param i iterator used to keep track of the selection
     */
    Bend* nextSelectedBend(int& i);
    /**
     * Finds the first unselected Bend.  In parallel, this is the
     * first unselected Bend that is the responsibility of the local
     * processor, not the first Bend in the global unselected pool.
     * @return a pointer to the Bend object, returns NULL if no Bend was found.
     * @param i iterator used to keep track of the selection
     */
    Bend* beginUnselectedBend(int& i);
    /**
     * Finds the next unselected Bend.  In parallel, this is the
     * next unselected Bend that is the responsibility of the local
     * processor, not the next Bend in the global unselected pool.
     * @return a pointer to the Bend object, returns NULL if no Bend was found.
     * @param i iterator used to keep track of the selection
     */
    Bend* nextUnselectedBend(int& i);
    /**
     * Finds the first selected Torsion in the selection.  In parallel,
     * this is the first selected Torsion that is the responsibility of
     * the local processor, not the first Torsion in the global
     * selection.
     * @return a pointer to the Torsion object, returns NULL if no Torsion was found.
     * @param i iterator used to keep track of the selection
     */
    Torsion* beginSelectedTorsion(int& i);
    /**
     * Finds the next selected Torsion in the selection.  In parallel,
     * this is the next selected Torsion that is the responsibility of
     * the local processor, not the next Torsion in the global selection.
     * @return a pointer to the Torsion object, returns NULL if no Torsion was found.
     * @param i iterator used to keep track of the selection
     */
    Torsion* nextSelectedTorsion(int& i);
    /**
     * Finds the first unselected Torsion.  In parallel, this is the
     * first unselected Torsion that is the responsibility of the local
     * processor, not the first Torsion in the global unselected pool.
     * @return a pointer to the Torsion object, returns NULL if no Torsion was found.
     * @param i iterator used to keep track of the selection
     */
    Torsion* beginUnselectedTorsion(int& i);
    /**
     * Finds the next unselected Torsion.  In parallel, this is the
     * next unselected Torsion that is the responsibility of the local
     * processor, not the next Torsion in the global unselected pool.
     * @return a pointer to the Torsion object, returns NULL if no Torsion was found.
     * @param i iterator used to keep track of the selection
     */
    Torsion* nextUnselectedTorsion(int& i);
    /**
     * Finds the first selected Inversion in the selection.  In parallel,
     * this is the first selected Inversion that is the responsibility of
     * the local processor, not the first Inversion in the global
     * selection.
     * @return a pointer to the Inversion object, returns NULL if no Inversion was found.
     * @param i iterator used to keep track of the selection
     */
    Inversion* beginSelectedInversion(int& i);
    /**
     * Finds the next selected Inversion in the selection.  In parallel,
     * this is the next selected Inversion that is the responsibility of
     * the local processor, not the next Inversion in the global selection.
     * @return a pointer to the Inversion object, returns NULL if no Inversion was found.
     * @param i iterator used to keep track of the selection
     */
    Inversion* nextSelectedInversion(int& i);
    /**
     * Finds the first unselected Inversion.  In parallel, this is the
     * first unselected Inversion that is the responsibility of the local
     * processor, not the first Inversion in the global unselected pool.
     * @return a pointer to the Inversion object, returns NULL if no Inversion was found.
     * @param i iterator used to keep track of the selection
     */
    Inversion* beginUnselectedInversion(int& i);
    /**
     * Finds the next unselected Inversion.  In parallel, this is the
     * next unselected Inversion that is the responsibility of the local
     * processor, not the next Inversion in the global unselected pool.
     * @return a pointer to the Inversion object, returns NULL if no Inversion was found.
     * @param i iterator used to keep track of the selection
     */
    Inversion* nextUnselectedInversion(int& i);
    /**
     * Finds the first selected Molecule in the selection.  In parallel,
     * this is the first selected Molecule that is the responsibility of
     * the local processor, not the first Molecule in the global
     * selection.
     * @return a pointer to the Molecule object, returns NULL if no Molecule was found.
     * @param i iterator used to keep track of the selection
     */
    Molecule* beginSelectedMolecule(int& i);
    /**
     * Finds the next selected Molecule in the selection.  In parallel,
     * this is the next selected Molecule that is the responsibility of
     * the local processor, not the next Molecule in the global selection.
     * @return a pointer to the Molecule object, returns NULL if no Molecule was found.
     * @param i iterator used to keep track of the selection
     */
    Molecule* nextSelectedMolecule(int& i);
    /**
     * Finds the first unselected Molecule.  In parallel, this is the
     * first unselected Molecule that is the responsibility of the local
     * processor, not the first Molecule in the global unselected pool.
     * @return a pointer to the Molecule object, returns NULL if no Molecule was found.
     * @param i iterator used to keep track of the selection
     */
    Molecule* beginUnselectedMolecule(int& i);
    /**
     * Finds the next unselected Molecule.  In parallel, this is the
     * next unselected Molecule that is the responsibility of the local
     * processor, not the next Molecule in the global unselected pool.
     * @return a pointer to the Molecule object, returns NULL if no Molecule was found.
     * @param i iterator used to keep track of the selection
     */
    Molecule* nextUnselectedMolecule(int& i);

    std::set<AtomType*> getSelectedAtomTypes();

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
    SimInfo* info_ {nullptr};
    SelectionSet ss_;
    std::vector<int> nObjects_;
    std::vector<StuntDouble*> stuntdoubles_;
    std::vector<Bond*> bonds_;
    std::vector<Bend*> bends_;
    std::vector<Torsion*> torsions_;
    std::vector<Inversion*> inversions_;
    std::vector<Molecule*> molecules_;
  };

}
#endif
