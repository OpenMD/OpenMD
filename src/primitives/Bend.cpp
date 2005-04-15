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
 
#include "primitives/Bend.hpp"

namespace oopse {

  /**@todo still a lot left to improve*/
  void Bend::calcForce() {
    Vector3d pos1 = atom1_->getPos();
    Vector3d pos2 = atom2_->getPos();
    Vector3d pos3 = atom3_->getPos();

    Vector3d r21 = pos1 - pos2;
    double d21 = r21.length();

    double d21inv = 1.0 / d21;

    Vector3d r23 = pos3 - pos2;
    double d23 = r23.length();

    double d23inv = 1.0 / d23;

    double cosTheta = dot(r21, r23) / (d21 * d23);

    //check roundoff     
    if (cosTheta > 1.0) {
      cosTheta = 1.0;
    } else if (cosTheta < -1.0) {
      cosTheta = -1.0;
    }

    double theta = acos(cosTheta);

    double dVdTheta;

    bendType_->calcForce(theta, potential_, dVdTheta);

    double sinTheta = sqrt(1.0 - cosTheta * cosTheta);

    if (fabs(sinTheta) < 1.0E-6) {
      sinTheta = 1.0E-6;
    }

    double commonFactor1 = dVdTheta / sinTheta * d21inv;
    double commonFactor2 = dVdTheta / sinTheta * d23inv;

    Vector3d force1 = commonFactor1 * (r23 * d23inv - r21*d21inv*cosTheta);
    Vector3d force3 = commonFactor2 * (r21 * d21inv - r23*d23inv*cosTheta);

    //total force in current bend is zero
    Vector3d force2 = force1 + force3;
    force2 *= -1.0;

    atom1_->addFrc(force1);
    atom2_->addFrc(force2);
    atom3_->addFrc(force3);
  }

} //end namespace oopse
