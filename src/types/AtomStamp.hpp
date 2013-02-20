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

#ifndef TYPES_ATOMSTAMP_HPP
#define TYPES_ATOMSTAMP_HPP

#include <set>
#include <vector>

#include "types/DataHolder.hpp"
namespace OpenMD {
  
  class AtomStamp  : public DataHolder {
    DeclareParameter(Type, std::string);
  public:
    AtomStamp(int index);
  public:
    
    bool setPosition(const std::vector<RealType>& pos);
    bool setOrientation(const std::vector<RealType>& ort);
    bool havePosition() { return havePos_; }
    bool haveOrientation() { return haveOrt_; }      
    RealType getPosX() { return position_[0]; }
    RealType getPosY() { return position_[1]; }
    RealType getPosZ() { return position_[2]; }
    RealType getEulerPhi()   { return orientation_[0]; }
    RealType getEulerTheta() { return orientation_[1]; }
    RealType getEulerPsi()   { return orientation_[2]; }
    int getIndex() { return index_;}
    virtual void validate();
    typedef std::set<int>::iterator AtomIter;
    typedef std::vector<int>::iterator BondIter;
    int getFirstBondedAtom(AtomIter& ai) {
      ai = bondedAtoms_.begin();
      return ai != bondedAtoms_.end() ? *ai : -1;
    }
    int getNextBondedAtom(AtomIter& ai) {
      ++ai;
      return ai != bondedAtoms_.end() ? *ai : -1;
    }
    int getFirstBond(BondIter& bi) {
      bi = bonds_.begin();
      return bi != bonds_.end()? *bi: -1;
    }
    int getNextBond(BondIter& bi) {
      ++bi; 
      return bi != bonds_.end()? *bi: -1;
    }
    void addBond(int bondIndex) {bonds_.push_back(bondIndex);}
    void addBondedAtom(int atomIndex) {bondedAtoms_.insert(atomIndex);}
    int getBondCount() { return bonds_.size(); }
  private:
    Vector3d position_;
    Vector3d orientation_;
    bool havePos_;
    bool haveOrt_;
    int index_;
    std::vector<int> bonds_;
    std::set<int> bondedAtoms_;
  };
  
}
#endif
