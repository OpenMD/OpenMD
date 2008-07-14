/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
 
#include "primitives/Inversion.hpp"

namespace oopse {

  Inversion::Inversion(Atom *atom1, Atom *atom2, Atom *atom3, 
		       Atom *atom4, InversionType *it) :
    atom1_(atom1), atom2_(atom2), atom3_(atom3), atom4_(atom4), 
    inversionType_(it) { }
  
  void Inversion::calcForce(RealType& angle) {
    
    // In OOPSE's version of an inversion, the central atom
    // comes first.  However, to get the planarity in a typical cosine
    // version of this potential (i.e. Amber-style), the central atom
    // is treated as atom *3* in a standard torsion form:

    Vector3d pos1 = atom2_->getPos();
    Vector3d pos2 = atom3_->getPos();
    Vector3d pos3 = atom1_->getPos();
    Vector3d pos4 = atom4_->getPos();

    Vector3d r21 = pos1 - pos2;
    Vector3d r32 = pos2 - pos3;
    Vector3d r43 = pos3 - pos4;

    //  Calculate the cross products and distances
    Vector3d A = cross(r21, r32);
    RealType rA = A.length();
    Vector3d B = cross(r32, r43);
    RealType rB = B.length();
    Vector3d C = cross(r32, A);
    RealType rC = C.length();

    A.normalize();
    B.normalize();
    C.normalize();
    
    //  Calculate the sin and cos
    RealType cos_phi = dot(A, B) ;
    if (cos_phi > 1.0) cos_phi = 1.0;
    if (cos_phi < -1.0) cos_phi = -1.0; 


    RealType dVdcosPhi;
    inversionType_->calcForce(cos_phi, potential_, dVdcosPhi);
    Vector3d f1;
    Vector3d f2;
    Vector3d f3;

    Vector3d dcosdA = (cos_phi * A - B) /rA;
    Vector3d dcosdB = (cos_phi * B - A) /rB;

    f1 = dVdcosPhi * cross(r32, dcosdA);
    f2 = dVdcosPhi * ( cross(r43, dcosdB) - cross(r21, dcosdA));
    f3 = dVdcosPhi * cross(dcosdB, r32);

    // In OOPSE's version of an improper torsion, the central atom
    // comes first.  However, to get the planarity in a typical cosine
    // version of this potential (i.e. Amber-style), the central atom
    // is treated as atom *3* in a standard torsion form:

    //  AMBER:   I - J - K - L   (e.g. K is sp2 hybridized carbon)
    //  OOPSE:   I - (J - K - L)  (e.g. I is sp2 hybridized carbon)

    // Confusing enough?  Good.

    atom3_->addFrc(f1);
    atom1_->addFrc(f2 - f1);
    atom2_->addFrc(f3 - f2);
    atom4_->addFrc(-f3);
    angle = acos(cos_phi) /M_PI * 180.0;
  }

}
