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

#ifndef HYDRODYNAMICS_HYDROPROP_HPP
#define HYDRODYNAMICS_HYDROPROP_HPP

#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"

namespace OpenMD {

  /**
   * @class HydroProp HydroProp.hpp "hydrodynamics/HydroProp.hpp"
   * Container for information about the hydrodynamic behavior of objects
   * interacting with surroundings.
   * @note the units for the center of resistance (a location) are in Angstroms
   *
   * @note the units of the resistance tensor are:
   *    Xitt (Translation-translation) : kcal fs mol^-1 Angstroms^-2
   *    Xirt (Rotation-translation) : kcal fs mol^-1 Angstroms^-1 radians^-1
   *    Xitr (Translation-rotation) : kcal fs mol^-1 Angstroms^-1 radians^-1
   *    Xirr (Rotation-rotation) : kcal fs mol^-1 radians^-2
   *
   * @note the units of D, the diffusion tensor are:
   *    Dtt (Translation-translation) : Angstroms^2 fs^-1
   *    Drt (Rotation-translation) :    Angstroms fs^-1
   *    Dtr (Translation-rotation) :    Angstroms fs^-1
   *    Drr (Rotation-rotation) :       fs^-1
   *
   * @note after setting the value of Xi manually, the complete() function
   * should be called to perform the Cholesky Decomposition.
   */
  class HydroProp {
  public:
    HydroProp();
    HydroProp(Vector3d cor, Mat6x6d Xi);

    void setName(std::string name) { name_ = name; }
    std::string getName() { return name_; }

    void setCenterOfResistance(Vector3d cor) {
      cor_    = cor;
      hasCOR_ = true;
    }
    Vector3d getCenterOfResistance() { return cor_; }

    void setResistanceTensor(Mat6x6d Xi) {
      Xi_    = Xi;
      hasXi_ = true;
      complete();
    }
    Mat6x6d getResistanceTensor() { return Xi_; }
    Mat3x3d getXitt();
    Mat3x3d getXitr();
    Mat3x3d getXirt();
    Mat3x3d getXirr();

    void complete();
    /*
     * Returns the result of a Cholesky decomposition to obtain the
     * square root matrix of the resistance tensor,
     *  \f[ \Xi = S S^T \f]
     * where S is a lower triangular matrix.
     */
    Mat6x6d getS();

    Mat6x6d getDiffusionTensor(RealType temperature);

    /*
     * Recomputes the Resistance Tensor at a new location
     */
    Mat6x6d getResistanceTensorAtPos(Vector3d pos);
    /*
     * Recomputes the Diffusion Tensor at a new location (and temperature)
     */
    Mat6x6d getDiffusionTensorAtPos(Vector3d pos, RealType temperature);
    /*
     * Computes the Center Of Diffusion at a particular temperature
     */
    Vector3d getCenterOfDiffusion(RealType temperature);

    Mat3x3d getPitchMatrix();
    RealType getScalarPitch();
    void pitchAxes(Mat3x3d& pitchAxes, Vector3d& pitches,
                   RealType& pitchScalar);

    /*
     * Computes the Center of Pitch
     */
    Vector3d getCenterOfPitch();

  private:
    std::string name_;
    Vector3d cor_;  // Center of Resistance
    Mat6x6d Xi_;    // Resistance Tensor
    Mat6x6d S_;     // Cholesky Decomposition of Xi_
    bool hasCOR_;
    bool hasXi_;
    bool hasS_;
  };
}  // namespace OpenMD

#endif
