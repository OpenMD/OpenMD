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
 
#ifndef TYPES_MOLECULESTAMP_HPP
#define TYPES_MOLECULESTAMP_HPP
#include <vector>
#include <utility>
#include "types/AtomStamp.hpp"
#include "types/BondStamp.hpp"
#include "types/BendStamp.hpp"
#include "types/TorsionStamp.hpp"
#include "types/InversionStamp.hpp"
#include "types/RigidBodyStamp.hpp"
#include "types/CutoffGroupStamp.hpp"
#include "types/FragmentStamp.hpp"

namespace OpenMD {
  class MoleculeStamp : public DataHolder {
    DeclareParameter(Name, std::string);
    DeclareParameter(ConstrainTotalCharge, bool);

  public:
    MoleculeStamp();
    virtual ~MoleculeStamp();
    
    bool addAtomStamp( AtomStamp* atom);
    bool addBondStamp( BondStamp* bond);
    bool addBendStamp( BendStamp* bend);
    bool addTorsionStamp( TorsionStamp* torsion);  
    bool addInversionStamp( InversionStamp* inversion);  
    bool addRigidBodyStamp( RigidBodyStamp* rigidbody);
    bool addCutoffGroupStamp( CutoffGroupStamp* cutoffgroup);
    bool addFragmentStamp( FragmentStamp* fragment);
    
    int  getNAtoms() { return atomStamps_.size(); }
    int  getNBonds() { return bondStamps_.size(); }
    int  getNBends() { return bendStamps_.size(); }
    int  getNTorsions() { return torsionStamps_.size(); }
    int  getNInversions() { return inversionStamps_.size(); }
    int  getNRigidBodies() { return rigidBodyStamps_.size(); }
    int  getNCutoffGroups() { return cutoffGroupStamps_.size(); }  
    int getNIntegrable() { return nintegrable_;}
    int getNFreeAtoms() { return freeAtoms_.size(); }
    virtual void validate();
    
    AtomStamp* getAtomStamp(int index) { return atomStamps_[index]; }
    BondStamp* getBondStamp(int index) { return bondStamps_[index]; }
    BendStamp* getBendStamp(int index) { return bendStamps_[index]; }
    TorsionStamp* getTorsionStamp(int index) { return torsionStamps_[index]; }
    InversionStamp* getInversionStamp(int index) { return inversionStamps_[index]; }
    RigidBodyStamp* getRigidBodyStamp(int index) { return rigidBodyStamps_[index]; }
    CutoffGroupStamp* getCutoffGroupStamp(int index) { return cutoffGroupStamps_[index]; }
    FragmentStamp* getFragmentStamp(int index) { return fragmentStamps_[index]; }
    
    bool isBondInSameRigidBody(BondStamp*bond);
    bool isAtomInRigidBody(int atomIndex);  
    bool isAtomInRigidBody(int atomIndex, int& whichRigidBody, 
                           int& consAtomIndex);  
    std::vector<std::pair<int, int> > getJointAtoms(int rb1, int rb2);
    
  private:
    
    void fillBondInfo();
    void checkAtoms();
    void checkBonds();
    void checkBends();
    void checkTorsions();
    void checkInversions();
    void checkRigidBodies();
    void checkCutoffGroups();
    void checkFragments();

    template <class Cont, class T>
    bool addIndexSensitiveStamp(Cont& cont, T* stamp) {
      unsigned int index = stamp->getIndex();
      bool ret = false;
      size_t size = cont.size();
      
      if (size >= index +1) {
        if (cont[index]!= NULL) {
          ret = false;
        }else {
          cont[index] = stamp;
          ret = true;
        }
      } else {
        cont.insert(cont.end(), index - cont.size() + 1, NULL);
        cont[index] = stamp;
        ret = true;
      }
      
      return ret;
    }
    
    std::vector<AtomStamp*> atomStamps_;
    std::vector<int> freeAtoms_;
    std::vector<BondStamp*> bondStamps_;
    std::vector<BendStamp*> bendStamps_;
    std::vector<TorsionStamp*> torsionStamps_;
    std::vector<InversionStamp*> inversionStamps_;
    std::vector<RigidBodyStamp*> rigidBodyStamps_;
    std::vector<CutoffGroupStamp*> cutoffGroupStamps_;
    std::vector<FragmentStamp*> fragmentStamps_;
    std::vector<int> atom2Rigidbody;
    int nintegrable_;
  };
  
}
#endif
