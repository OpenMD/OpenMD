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
 
/**
 * @file Bend.hpp
 * @author    tlin
 * @date  11/01/2004
 * @version 1.0
 */ 


#ifndef PRIMITIVES_BEND_HPP
#define PRIMITIVES_BEND_HPP

#include "primitives/ShortRangeInteraction.hpp"
#include "primitives/Atom.hpp"
#include "types/BendType.hpp"

namespace OpenMD {
  struct BendData {
    RealType angle;
    RealType potential;
  };
  
  struct BendDataSet {
    RealType deltaV;
    BendData prev;
    BendData curr;
  };
  
  class Bend : public ShortRangeInteraction {
  public:
    using ShortRangeInteraction::getValue;
    using ShortRangeInteraction::getPrevValue;
    Bend(Atom* atom1, Atom* atom2, Atom* atom3, BendType* bt)
      : ShortRangeInteraction(), bendType_(bt) {
      atoms_.resize(3);
      atoms_[0] = atom1;
      atoms_[1] = atom2;
      atoms_[2] = atom3;
    }
    
    virtual ~Bend() {}
    virtual void calcForce(RealType& angle, bool doParticlePot);
    
    RealType getValue(int snapshotNo) {
      Vector3d pos1 = atoms_[0]->getPos(snapshotNo);
      Vector3d pos2 = atoms_[1]->getPos(snapshotNo);
      Vector3d pos3 = atoms_[2]->getPos(snapshotNo);
      
      Vector3d r21 = pos1 - pos2;
      RealType d21 = r21.length();
      
      Vector3d r23 = pos3 - pos2;
      RealType d23 = r23.length();
      
      RealType cosTheta = dot(r21, r23) / (d21 * d23);
      
      //check roundoff     
      if (cosTheta > 1.0) {
        cosTheta = 1.0;
      } else if (cosTheta < -1.0) {
        cosTheta = -1.0;
      }
      
      return acos(cosTheta);
    }


    RealType getPotential() {
      return potential_;
    }
    
    Atom* getAtomA() {
      return atoms_[0];
    }
    
    Atom* getAtomB() {
      return atoms_[1];
    }
    
    Atom* getAtomC() {
      return atoms_[2];
    }
    
    BendType * getBendType() {
      return bendType_;
    }
    
    virtual std::string getName() { return name_;}        
    /** Sets the name of this bend for selections */
    virtual void setName(const std::string& name) { name_ = name;}

    void accept(BaseVisitor* v) {
      v->visit(this);
    }    

  protected:
    
    RealType potential_;
    BendType* bendType_; /**< bend type */
    std::string name_;        

  };    
} //end namespace OpenMD
#endif //PRIMITIVES_BEND_HPP

