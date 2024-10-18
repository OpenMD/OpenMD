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

/*! \file perturbations/Light.hpp
    \brief Rotating Electric Field perturbation
*/

#ifndef PERTURBATIONS_LIGHT_HPP
#define PERTURBATIONS_LIGHT_HPP

#include <complex>

#include "brains/ForceModifier.hpp"
#include "brains/SimInfo.hpp"
#include "perturbations/LightParameters.hpp"

namespace OpenMD::Perturbations {

  //! Applies an EM field representing light to the system
  /*! The field is applied as an external perturbation.  The user specifies

  \code{.unparsed}
    light{
      useLight = true;
      intensity = c;
      wavelength = l;
      propagationDirection = (0, 0, 1);
      polarization = "+";
    }
  \endcode

    in the .omd file where the propagationDirection vector, \f$
    \mathbf{k} \f$ is a unit vector, the value of \f$ c \f$ is in
    units of \f$ W cm^{-2} \f$, and wavelength, \f$ l \f$ is in nm.

    Alternatively, the user can specify

  \code{.unparsed}
    light{
      useLight = true;
      intensity = c;
      frequency = w;
      propagationDirection = (0, 0, 1);
      polarization = "+";
    }
  \endcode

   with frequency, \f$ w \f$ measured in Hz, or alternatively:

  \code{.unparsed}
    light{
      useLight = true;
      intensity = c;
      waveVector = (kx, ky, kz);
      polarization = "+";
    }
  \endcode

  with the components of the wavevector specified in \f$ \AA^{-1} \f$.
  In all cases, options for polarization include "x", "y", "+", or "-".

  */
  class Light : public ForceModifier {
  public:
    Light(SimInfo* info);

  private:
    void modifyForces() override;
    void initialize();

    bool initialized;
    bool doLight;
    bool doParticlePot;

    SimInfo* info_;
    Perturbations::LightParameters* lightParams;
    Vector3d khat_;
    Vector3d k_;
    RealType omega_;
    RealType lambda_;
    RealType kmag_;
    RealType E0_;
    std::vector<std::complex<RealType>>
        jones_;  // Jones vector for polarization
    RotMat3x3d A_;
    RotMat3x3d Ainv_;

    enum LightPolarization {
      lightX,      // x-plane polarized light
      lightY,      // y-plane polarized light
      lightPlus,   // circularly polarized light (clockwise)
      lightMinus,  // circularly polarized light (counterclockwise)
      lightUnknownPolarization
    };
  };
}  // namespace OpenMD::Perturbations

#endif
