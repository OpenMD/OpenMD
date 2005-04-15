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
 
#include "primitives/Torsion.hpp"

namespace oopse {

  Torsion::Torsion(Atom *atom1, Atom *atom2, Atom *atom3, Atom *atom4,
		   TorsionType *tt) :
    atom1_(atom1), atom2_(atom2), atom3_(atom3), atom4_(atom4), torsionType_(tt) { }

  void Torsion::calcForce() {
    Vector3d pos1 = atom1_->getPos();
    Vector3d pos2 = atom2_->getPos();
    Vector3d pos3 = atom3_->getPos();
    Vector3d pos4 = atom4_->getPos();

    Vector3d r21 = pos1 - pos2;
    Vector3d r32 = pos2 - pos3;
    Vector3d r43 = pos3 - pos4;

    //  Calculate the cross products and distances
    Vector3d A = cross(r21, r32);
    double rA = A.length();
    Vector3d B = cross(r32, r43);
    double rB = B.length();
    Vector3d C = cross(r32, A);
    double rC = C.length();

    A.normalize();
    B.normalize();
    C.normalize();
    
    //  Calculate the sin and cos
    double cos_phi = dot(A, B) ;
    double sin_phi = dot(C, B);

    double dVdPhi;
    torsionType_->calcForce(cos_phi, sin_phi, potential_, dVdPhi);

    Vector3d f1;
    Vector3d f2;
    Vector3d f3;

    //  Next, we want to calculate the forces.  In order
    //  to do that, we first need to figure out whether the
    //  sin or cos form will be more stable.  For this,
    //  just look at the value of phi
    //if (fabs(sin_phi) > 0.1) {
    //  use the sin version to avoid 1/cos terms

    Vector3d dcosdA = (cos_phi * A - B) /rA;
    Vector3d dcosdB = (cos_phi * B - A) /rB;

    double dVdcosPhi = -dVdPhi / sin_phi;

    f1 = dVdcosPhi * cross(r32, dcosdA);
    f2 = dVdcosPhi * ( cross(r43, dcosdB) - cross(r21, dcosdA));
    f3 = dVdcosPhi * cross(dcosdB, r32);

    /** @todo fix below block, must be something wrong with the sign somewhere */
    //} else {
    //  This angle is closer to 0 or 180 than it is to
    //  90, so use the cos version to avoid 1/sin terms

    //double dVdsinPhi = dVdPhi /cos_phi;
    //Vector3d dsindB = (sin_phi * B - C) /rB;
    //Vector3d dsindC = (sin_phi * C - B) /rC;

    //f1.x() = dVdsinPhi*((r32.y()*r32.y() + r32.z()*r32.z())*dsindC.x() - r32.x()*r32.y()*dsindC.y() - r32.x()*r32.z()*dsindC.z());

    //f1.y() = dVdsinPhi*((r32.z()*r32.z() + r32.x()*r32.x())*dsindC.y() - r32.y()*r32.z()*dsindC.z() - r32.y()*r32.x()*dsindC.x());

    //f1.z() = dVdsinPhi*((r32.x()*r32.x() + r32.y()*r32.y())*dsindC.z() - r32.z()*r32.x()*dsindC.x() - r32.z()*r32.y()*dsindC.y());

    //f2.x() = dVdsinPhi*(-(r32.y()*r21.y() + r32.z()*r21.z())*dsindC.x() + (2.0*r32.x()*r21.y() - r21.x()*r32.y())*dsindC.y()
    //+ (2.0*r32.x()*r21.z() - r21.x()*r32.z())*dsindC.z() + dsindB.z()*r43.y() - dsindB.y()*r43.z());

    //f2.y() = dVdsinPhi*(-(r32.z()*r21.z() + r32.x()*r21.x())*dsindC.y() + (2.0*r32.y()*r21.z() - r21.y()*r32.z())*dsindC.z()
    //+ (2.0*r32.y()*r21.x() - r21.y()*r32.x())*dsindC.x() + dsindB.x()*r43.z() - dsindB.z()*r43.x());

    //f2.z() = dVdsinPhi*(-(r32.x()*r21.x() + r32.y()*r21.y())*dsindC.z() + (2.0*r32.z()*r21.x() - r21.z()*r32.x())*dsindC.x()
    //+(2.0*r32.z()*r21.y() - r21.z()*r32.y())*dsindC.y() + dsindB.y()*r43.x() - dsindB.x()*r43.y());

    //f3 = dVdsinPhi * cross(r32, dsindB);

    //}

    atom1_->addFrc(f1);
    atom2_->addFrc(f2 - f1);
    atom3_->addFrc(f3 - f2);
    atom4_->addFrc(-f3);
  }

}
