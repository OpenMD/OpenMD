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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include "DLM.hpp"

namespace OpenMD {

  void DLM::doRotate(StuntDouble* sd, Vector3d& ji, RealType dt) {
    RealType dt2 = 0.5 * dt;    
    RealType angle;

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


  void DLM::rotateStep(int axes1, int axes2, RealType angle, Vector3d& ji, RotMat3x3d& A) {

    RealType sinAngle;
    RealType cosAngle;
    RealType angleSqr;
    RealType angleSqrOver4;
    RealType top, bottom;

    RotMat3x3d tempA(A);  // initialize the tempA
    Vector3d tempJ(0.0);

    RotMat3x3d rot = RotMat3x3d::identity(); // initalize rot as a unit matrix

    // use a small angle aproximation for sin and cosine
    
    angleSqr = angle * angle;
    angleSqrOver4 = angleSqr / 4.0;
    top = 1.0 - angleSqrOver4;
    bottom = 1.0 + angleSqrOver4;

    cosAngle = top / bottom;
    sinAngle = angle / bottom;
    
    // or don't use the small angle approximation:
    //cosAngle = cos(angle);
    //sinAngle = sin(angle);

    rot(axes1, axes1) = cosAngle;
    rot(axes2, axes2) = cosAngle;

    rot(axes1, axes2) = sinAngle;
    rot(axes2, axes1) = -sinAngle;

    // rotate the momentum acoording to: ji[] = rot[][] * ji[]
    ji = rot * ji;

    // This code comes from converting an algorithm detailed in 
    // J. Chem. Phys. 107 (15), pp. 5840-5851 by Dullweber, 
    // Leimkuhler and McLachlan (DLM) for use in our code.
    // In Appendix A, the DLM paper has the change to the rotation 
    // matrix as: Q = Q * rot.transpose(), but our rotation matrix 
    // A is actually equivalent to Q.transpose(). This fact can be 
    // seen on page 5849 of the DLM paper where a lab frame 
    // dipole \mu_i(t) is expressed in terms of a body-fixed
    // reference orientation, \bar{\mu_i} and the rotation matrix, Q:
    //  \mu_i(t) = Q * \bar{\mu_i}
    // Our code computes lab frame vectors from body-fixed reference
    // vectors using:
    //   v_{lab} = A.transpose() * v_{body}
    //  (See StuntDouble.hpp for confirmation of this fact).
    //
    // So, using the identity:
    //  (A * B).transpose() = B.transpose() * A.transpose(),  we
    // get the equivalent of Q = Q * rot.transpose() for our code to be:

    A = rot * A;
  
  }


}
