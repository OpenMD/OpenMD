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

/**
 * @file DirectionalAtom.hpp
 * @author    tlin
 * @date  10/23/2004
 * @version 1.0
 */

#ifndef PRIMITIVES_DIRECTIONALATOM_HPP
#define PRIMITIVES_DIRECTIONALATOM_HPP

#include "primitives/Atom.hpp"
#include "types/DirectionalAtomType.hpp"
namespace OpenMD {
  class DirectionalAtom : public Atom {
  public:
    DirectionalAtom(DirectionalAtomType* dAtomType);
    /**
     * Returns the inertia tensor of this stuntdouble
     * @return the inertia tensor of this stuntdouble
     */
    virtual Mat3x3d getI();

    /**
     * Sets  the previous rotation matrix of this stuntdouble
     * @param a  new rotation matrix
     */
    virtual void setPrevA(const RotMat3x3d& a);

    /**
     * Sets  the current rotation matrix of this stuntdouble
     * @param a  new rotation matrix
     */
    virtual void setA(const RotMat3x3d& a);

    /**
     * Sets  the rotation matrix of this stuntdouble in specified snapshot
     * @param a rotation matrix to be set
     * @param snapshotNo
     * @see #getA
     */
    virtual void setA(const RotMat3x3d& a, int snapshotNo);

    /**
     * Left multiple rotation matrix by another rotation matrix
     * @param m a rotation matrix
     */
    void rotateBy(const RotMat3x3d& m);

    /** Sets the internal unit frame of this stuntdouble by three euler angles
     */
    void setUnitFrameFromEuler(double phi, double theta, double psi);

    /**
     * Returns the gradient of this stuntdouble
     * @return the gradient of this stuntdouble
     */
    virtual std::vector<double> getGrad();

    virtual void accept(BaseVisitor* v);

  protected:
    Mat3x3d inertiaTensor_; /**< inertial tensor */
    RotMat3x3d sU_;         /**< body fixed standard unit vector */
  };

}  // namespace OpenMD

#endif  // PRIMITIVES_DIRECTIONALATOM_HPP
