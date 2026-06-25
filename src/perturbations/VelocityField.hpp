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

/*! \file perturbations/VelocityField.hpp
    \brief Linear (homogeneous) ambient velocity field, queried at a location
*/

#ifndef PERTURBATIONS_VELOCITYFIELD_HPP
#define PERTURBATIONS_VELOCITYFIELD_HPP

#include "brains/SimInfo.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "perturbations/VelocityFieldParameters.hpp"

namespace OpenMD {

  //! A linear (homogeneous) ambient velocity field, v(r) = v0 + K . r
  /*! Built from a velocityField{ ... } block.  The velocity gradient
      \f$ K_{ij} = \partial v_i / \partial r_j \f$ is assembled from an
      optional constant-stress (rate-of-strain) part and an optional
      constant-vorticity part:

      \f$ \mathsf{K} = \mathsf{E} + \mathsf{W}, \qquad
          \mathsf{E} = \mathrm{sym}(\mathsf{K}), \quad
          \mathsf{W} = \mathrm{antisym}(\mathsf{K}). \f$

      Unlike the electric-field gradient in UniformGradient (which Laplace's
      equation forces to be symmetric and traceless), K carries an
      independent antisymmetric part.  Because E is built traceless and W is
      antisymmetric, this decomposed field is divergence-free by
      construction (incompressible, \f$\nabla\cdot\mathbf{v}=0\f$).

      This is a passive object: it is queried by the (Langevin /
      hydrodynamic) integrator rather than registered as a ForceModifier.
      Accessors expose the quantities the mobility coupling consumes: the
      ambient velocity (leading drag), the rate of strain (stresslet), and
      the vorticity / co-rotation rate (rotational coupling).  Since the
      flow is linear, \f$\nabla^2\mathbf{v}=0\f$ and the a^2/6 Faxen term
      vanishes identically.

      The velocity field can be applied by specifying one of these blocks in the
      omd file:
      
      \code{.unparsed}
      
      // pure planar/uniaxial extension along x
      velocityField {
         useVelocityField = true;
         strainRate       = 1.0e-3;        // fs^-1
         strainDirection1 = (1, 0, 0);
         strainDirection2 = (1, 0, 0);
      }

      // rigid rotation about z (vorticity only)
      velocityField {
         useVelocityField = true;
         vorticity        = (0, 0, 2.0e-3);  // fs^-1
      }

      // background velocity field
      velocityField {
         useVelocityField   = true;
         backgroundVelocity = (0, 0, 1.0e-4);  // Angstrom fs^-1
      }

      // simple shear  v_x = gammaDot * y  (equal strain + vorticity)
      velocityField {
         useVelocityField = true;
         strainRate       = 1.0e-3;
         strainDirection1 = (1, 0, 0);
         strainDirection2 = (0, 1, 0);
         vorticity        = (0, 0, -1.0e-3);
      }
      \endcode

      This last block is important: simple shear is not a pure strain
      or pure vorticity state — it's \f$ K_{xy} = \f$ split as \f$
      E_{xy} = \dot{\gamma} / 2 \f$ plus \f$ \omega = -\dot{\gamma}
      \f$, which is why both strain and vorticity appear.
  */
  class VelocityField {
  public:
    explicit VelocityField(SimInfo* info);

    //! true when a valid velocityField{} block was supplied
    bool isActive() const { return doVelocityField_; }

    //! ambient velocity at position r:  v0 + K . r
    Vector3d getVelocity(const Vector3d& r) const { return v0_ + K_ * r; }

    //! constant velocity gradient, K_ij = d v_i / d r_j
    const Mat3x3d& getVelocityGradient() const { return K_; }

    //! rate-of-strain tensor E = sym(K)  (the constant-stress part)
    const Mat3x3d& getRateOfStrain() const { return E_; }

    //! spin tensor W = antisym(K)
    const Mat3x3d& getSpin() const { return W_; }

    //! vorticity vector, omega = curl v  (the constant-vorticity part)
    const Vector3d& getVorticity() const { return omega_; }

    //! angular velocity an immersed sphere co-rotates with, omega / 2
    Vector3d getAngularVelocity() const { return 0.5 * omega_; }

    //! uniform background velocity v0
    const Vector3d& getBackgroundVelocity() const { return v0_; }

  private:
    void initialize();

    SimInfo* info_;
    Perturbations::VelocityFieldParameters* vfParams_;

    bool initialized_ {false};
    bool doVelocityField_ {false};

    Vector3d v0_ {V3Zero};     // uniform background velocity
    Mat3x3d K_ {};             // velocity gradient
    Mat3x3d E_ {};             // sym(K), rate of strain
    Mat3x3d W_ {};             // antisym(K), spin
    Vector3d omega_ {V3Zero};  // vorticity = curl v
  };
}  // namespace OpenMD

#endif  // PERTURBATIONS_VELOCITYFIELD_HPP
