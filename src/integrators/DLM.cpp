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
 
#include "DLM.hpp"

namespace oopse {

void DLM::doRotate(StuntDouble* sd, Vector3d& ji, double dt) {
    double dt2 = 0.5 * dt;    
    double angle;

    RotMat3x3d A = sd->getA();
    Mat3x3d I = sd->getI();

    // use the angular velocities to propagate the rotation matrix a full time step
    if (sd->isLinear()) {

        int i = sd->linearAxis();
        int j = (i+1)%3;
        int k = (i+2)%3;

        angle = dt2 * ji[j] / I(j, j);
        rotateStep( k, i, angle, ji, A );

        angle = dt * ji[k] / I(k, k);
        rotateStep( i, j, angle, ji, A);

        angle = dt2 * ji[j] / I(j, j);
        rotateStep( k, i, angle, ji, A );

    } else {
        // rotate about the x-axis
        angle = dt2 * ji[0] / I(0, 0);
        rotateStep( 1, 2, angle, ji, A );

        // rotate about the y-axis
        angle = dt2 * ji[1] / I(1, 1);
        rotateStep( 2, 0, angle, ji, A );

        // rotate about the z-axis
        angle = dt * ji[2] / I(2, 2);
        sd->addZangle(angle);
        rotateStep( 0, 1, angle, ji, A);

        // rotate about the y-axis
        angle = dt2 * ji[1] / I(1, 1);
        rotateStep( 2, 0, angle, ji, A );

        // rotate about the x-axis
        angle = dt2 * ji[0] / I(0, 0);
        rotateStep( 1, 2, angle, ji, A );

    }

    sd->setA( A  );
}


void DLM::rotateStep(int axes1, int axes2, double angle, Vector3d& ji, RotMat3x3d& A) {

    double sinAngle;
    double cosAngle;
    double angleSqr;
    double angleSqrOver4;
    double top, bottom;

    RotMat3x3d tempA(A);  // initialize the tempA
    Vector3d tempJ(0.0);

    RotMat3x3d rot = RotMat3x3d::identity(); // initalize rot as a unit matrix

    // use a small angle aproximation for sin and cosine

    //angleSqr = angle * angle;
    //angleSqrOver4 = angleSqr / 4.0;
    //top = 1.0 - angleSqrOver4;
    //bottom = 1.0 + angleSqrOver4;

    //cosAngle = top / bottom;
    //sinAngle = angle / bottom;
    cosAngle = cos(angle);
    sinAngle = sin(angle);
    rot(axes1, axes1) = cosAngle;
    rot(axes2, axes2) = cosAngle;

    rot(axes1, axes2) = sinAngle;
    rot(axes2, axes1) = -sinAngle;

    // rotate the momentum acoording to: ji[] = rot[][] * ji[]
    ji = rot * ji;

    // rotate the Rotation matrix acording to:
    // A[][] = A[][] * transpose(rot[][])
    // transpose(A[][]) = transpose(A[][]) * transpose(rot[][])

    A = rot * A; //? A = A* rot.transpose();
  
}


}
