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
 * @file Inversion.hpp
 * @author    tlin
 * @date  11/01/2004
 * @version 1.0
 */

#ifndef PRIMITIVES_INVERSION_HPP
#define PRIMITIVES_INVERSION_HPP

#include "primitives/ShortRangeInteraction.hpp"
#include "primitives/Atom.hpp"
#include "types/InversionType.hpp"

namespace OpenMD {
  struct InversionData {
    RealType angle;
    RealType potential;
  };

  struct InversionDataSet {
    RealType deltaV;
    InversionData prev;
    InversionData curr;
  };

  /**
   * @class Inversion Inversion.hpp "primitives/Inversion.hpp"
   */
  class Inversion : public ShortRangeInteraction {
  public:
    Inversion(Atom* atom1, Atom* atom2, Atom* atom3, Atom* atom4, InversionType* it);
    virtual ~Inversion() {}
    virtual void calcForce(RealType& angle, bool doParticlePot);
        
    RealType getValue(int snapshotNo) {
      // In OpenMD's version of an inversion, the central atom
      // comes first.  However, to get the planarity in a typical cosine
      // version of this potential (i.e. Amber-style), the central atom
      // is treated as atom *3* in a standard torsion form:
      
      Vector3d pos1 = atoms_[1]->getPos(snapshotNo);
      Vector3d pos2 = atoms_[2]->getPos(snapshotNo);
      Vector3d pos3 = atoms_[0]->getPos(snapshotNo);
      Vector3d pos4 = atoms_[3]->getPos(snapshotNo);
      
      Vector3d r31 = pos1 - pos3;
      Vector3d r23 = pos3 - pos2;
      Vector3d r43 = pos3 - pos4;
      
      //  Calculate the cross products and distances
      Vector3d A = cross(r31, r43);
      RealType rA = A.length();
      Vector3d B = cross(r43, r23);
      RealType rB = B.length();
      //Vector3d C = cross(r23, A);
      //RealType rC = C.length();
      
      A.normalize();
      B.normalize();
      //C.normalize();
      
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

    InversionType * getInversionType() {
      return inversionType_;
    }
    virtual std::string getName() { return name_;}        
    /** Sets the name of this inversion for selections */
    virtual void setName(const std::string& name) { name_ = name;}

    void accept(BaseVisitor* v) {
      v->visit(this);
    }    

  protected:
    InversionType* inversionType_;
    InversionKey inversionKey_;
    std::string name_;        

    RealType potential_;
  };    

}
#endif //PRIMITIVES_INVERSION_HPP
