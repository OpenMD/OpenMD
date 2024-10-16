/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#include "DLM.hpp"

namespace OpenMD {

  void DLM::doRotate(StuntDouble* sd, Vector3d& ji, RealType dt) {
    RealType dt2 = 0.5 * dt;
    RealType angle;

    RotMat3x3d A = sd->getA();
    Mat3x3d I    = sd->getI();

    // use the angular velocities to propagate the rotation matrix a full time
    // step
    if (sd->isLinear()) {
      int i = sd->linearAxis();
      int j = (i + 1) % 3;
      int k = (i + 2) % 3;

      angle = dt2 * ji[j] / I(j, j);
      rotateStep(k, i, angle, ji, A);

      angle = dt * ji[k] / I(k, k);
      rotateStep(i, j, angle, ji, A);

      angle = dt2 * ji[j] / I(j, j);
      rotateStep(k, i, angle, ji, A);

    } else {
      // rotate about the x-axis
      angle = dt2 * ji[0] / I(0, 0);
      rotateStep(1, 2, angle, ji, A);

      // rotate about the y-axis
      angle = dt2 * ji[1] / I(1, 1);
      rotateStep(2, 0, angle, ji, A);

      // rotate about the z-axis
      angle = dt * ji[2] / I(2, 2);
      rotateStep(0, 1, angle, ji, A);

      // rotate about the y-axis
      angle = dt2 * ji[1] / I(1, 1);
      rotateStep(2, 0, angle, ji, A);

      // rotate about the x-axis
      angle = dt2 * ji[0] / I(0, 0);
      rotateStep(1, 2, angle, ji, A);
    }

    sd->setA(A);
  }

  void DLM::rotateStep(int axes1, int axes2, RealType angle, Vector3d& ji,
                       RotMat3x3d& A) {
    RealType sinAngle;
    RealType cosAngle;
    // RealType angleSqr;
    // RealType angleSqrOver4;
    // RealType top, bottom;

    RotMat3x3d tempA(A);  // initialize the tempA
    Vector3d tempJ(0.0);

    RotMat3x3d rot = RotMat3x3d::identity();  // initalize rot as a unit matrix

    // use a small angle aproximation for sin and cosine

    // angleSqr = angle * angle;
    // angleSqrOver4 = angleSqr / 4.0;
    // top = 1.0 - angleSqrOver4;
    // bottom = 1.0 + angleSqrOver4;

    // cosAngle = top / bottom;
    // sinAngle = angle / bottom;

    // or don't use the small angle approximation:
    cosAngle          = cos(angle);
    sinAngle          = sin(angle);
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

}  // namespace OpenMD
