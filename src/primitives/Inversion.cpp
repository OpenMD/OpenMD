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
 
#include "config.h"
#include <cmath>

#include "primitives/Inversion.hpp"

namespace OpenMD {

  Inversion::Inversion(Atom *atom1, Atom *atom2, Atom *atom3, 
		       Atom *atom4, InversionType *it) :
    atom1_(atom1), atom2_(atom2), atom3_(atom3), atom4_(atom4), 
    inversionType_(it) { }
  
  void Inversion::calcForce(RealType& angle, bool doParticlePot) {
    
    // In OpenMD's version of an inversion, the central atom
    // comes first.  However, to get the planarity in a typical cosine
    // version of this potential (i.e. Amber-style), the central atom
    // is treated as atom *3* in a standard torsion form:

    Vector3d pos1 = atom2_->getPos();
    Vector3d pos2 = atom3_->getPos();
    Vector3d pos3 = atom1_->getPos();
    Vector3d pos4 = atom4_->getPos();

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

    RealType dVdcosPhi;
    inversionType_->calcForce(cos_phi, potential_, dVdcosPhi);
    Vector3d f1 ;
    Vector3d f2 ;
    Vector3d f3 ;

    Vector3d dcosdA = (cos_phi * A - B) /rA;
    Vector3d dcosdB = (cos_phi * B - A) /rB;

    f1 = dVdcosPhi * cross(r43, dcosdA);
    f2 = dVdcosPhi * ( cross(r23, dcosdB) - cross(r31, dcosdA));
    f3 = dVdcosPhi * cross(dcosdB, r43);

    // In OpenMD's version of an improper torsion, the central atom
    // comes first.  However, to get the planarity in a typical cosine
    // version of this potential (i.e. Amber-style), the central atom
    // is treated as atom *3* in a standard torsion form:

    //  AMBER:   I - J - K - L   (e.g. K is sp2 hybridized carbon)
    //  OpenMD:  I - (J - K - L)  (e.g. I is sp2 hybridized carbon)

    // Confusing enough?  Good.

    atom2_->addFrc(f1);
    atom1_->addFrc(f2 - f1 + f3);
    atom4_->addFrc(-f2);
    atom3_->addFrc(-f3);

    if (doParticlePot) { 
      atom1_->addParticlePot(potential_);
      atom2_->addParticlePot(potential_);
      atom3_->addParticlePot(potential_);
      atom4_->addParticlePot(potential_);
    }
    
    angle = acos(cos_phi) /M_PI * 180.0;
  }

}
