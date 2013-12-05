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
 * @file Torsion.hpp
 * @author    tlin
 * @date  11/01/2004
 * @version 1.0
 */

#ifndef PRIMITIVES_TORSION_HPP
#define PRIMITIVES_TORSION_HPP

#include "primitives/ShortRangeInteraction.hpp"
#include "primitives/Atom.hpp"
#include "types/TorsionType.hpp"
#include <limits>

namespace OpenMD {
struct TorsionData {
    RealType angle;
    RealType potential;
};

struct TorsionDataSet {
    RealType deltaV;
    TorsionData prev;
    TorsionData curr;
};


  /**
   * @class Torsion Torsion.hpp "types/Torsion.hpp"
   */
  class Torsion : public ShortRangeInteraction {
  public:
    Torsion(Atom* atom1, Atom* atom2, Atom* atom3, Atom* atom4, TorsionType* tt);
    virtual ~Torsion() {}
    virtual void calcForce(RealType& angle, bool doParticlePot);
                
    RealType getValue(int snapshotNo) {
      Vector3d pos1 = atoms_[0]->getPos(snapshotNo);
      Vector3d pos2 = atoms_[1]->getPos(snapshotNo);
      Vector3d pos3 = atoms_[2]->getPos(snapshotNo);
      Vector3d pos4 = atoms_[3]->getPos(snapshotNo);
      
      Vector3d r21 = pos1 - pos2;
      Vector3d r32 = pos2 - pos3;
      Vector3d r43 = pos3 - pos4;
      
      //  Calculate the cross products and distances
      Vector3d A = cross(r21, r32);
      RealType rA = A.length();
      Vector3d B = cross(r32, r43);
      RealType rB = B.length();
      
      /* 
         If either of the two cross product vectors is tiny, that means
         the three atoms involved are colinear, and the torsion angle is
         going to be undefined.  The easiest check for this problem is
         to use the product of the two lengths.
      */
      if (rA * rB < OpenMD::epsilon) return numeric_limits<double>::quiet_NaN();
      
      A.normalize();
      B.normalize();  
      
      //  Calculate the sin and cos
      RealType cos_phi = dot(A, B) ;
      if (cos_phi > 1.0) cos_phi = 1.0;
      if (cos_phi < -1.0) cos_phi = -1.0;            
      return acos(cos_phi);
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

    Atom* getAtomD() {
      return atoms_[3];
    }

    TorsionType * getTorsionType() {
      return torsionType_;
    }
        
    virtual std::string getName() { return name_;}        
    /** Sets the name of this torsion for selections */
    virtual void setName(const std::string& name) { name_ = name;}

    void accept(BaseVisitor* v) {
      v->visit(this);
    }    

  protected:

    TorsionType* torsionType_;
    std::string name_;        

    RealType potential_;
  };    

}
#endif //PRIMITIVES_TORSION_HPP
