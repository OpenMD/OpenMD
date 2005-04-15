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
 
#include "primitives/GhostTorsion.hpp"

namespace oopse {

  GhostTorsion::GhostTorsion(Atom *atom1, Atom *atom2,  DirectionalAtom* ghostAtom,
			     TorsionType *tt) : Torsion(atom1, atom2, ghostAtom, ghostAtom, tt) {}

  void GhostTorsion::calcForce() {
    DirectionalAtom* ghostAtom = static_cast<DirectionalAtom*>(atom3_);    

    Vector3d pos1 = atom1_->getPos();
    Vector3d pos2 = atom2_->getPos();
    Vector3d pos3 = ghostAtom->getPos();

    Vector3d r21 = pos1 - pos2;
    Vector3d r32 = pos2 - pos3;
    Vector3d r43 = ghostAtom->getElectroFrame().getColumn(2);

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

    Vector3d dcosdA = (cos_phi * A - B) /rA;
    Vector3d dcosdB = (cos_phi * B - A) /rB;

    double dVdcosPhi = -dVdPhi / sin_phi;

    Vector3d f1 = dVdcosPhi * cross(r32, dcosdA);
    Vector3d f2 = dVdcosPhi * ( cross(r43, dcosdB) - cross(r21, dcosdA));
    Vector3d f3 = dVdcosPhi * cross(dcosdB, r32);

    atom1_->addFrc(f1);
    atom2_->addFrc(f2 - f1);

    ghostAtom->addFrc(-f2);

    f3.negate();
    ghostAtom->addTrq(cross(r43, f3));    
  }

}

