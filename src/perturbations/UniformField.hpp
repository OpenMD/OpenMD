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
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

/*! \file perturbations/UniformField.hpp
    \brief Uniform Electric Field perturbation
*/

#ifndef PERTURBATIONS_UNIFORMFIELD_HPP
#define PERTURBATIONS_UNIFORMFIELD_HPP

#include "brains/ForceModifier.hpp"
#include "brains/SimInfo.hpp"

namespace OpenMD {

  //! Applies a uniform (vector) electric field to the system
  /*! The field is applied as an external perturbation.  The user specifies

  \code{.unparsed}
    uniformField = (a, b, c);
  \endcode

    in the .omd file where the values of a, b, and c are in units of
    \f$ V / \AA \f$

    The electrostatic potential corresponding to this uniform field is

    \f$ \phi(\mathbf{r})  = - a x - b y - c z \f$

    which grows unbounded and is not periodic.  For these reasons,
    care should be taken in using a Uniform field with point charges.

    The field itself is

    \f$ \mathbf{E} = \left( \begin{array}{c} a \\ b \\ c \end{array} \right) \f$

   The external field applies a force on charged atoms, \f$ \mathbf{F}
   = C \mathbf{E} \f$.  For dipolar atoms, the field applies both a
   potential, \f$ U = - \mathbf{D} \cdot \mathbf{E} \f$ and a torque,
   \f$ \mathbf{\tau} = \mathbf{D} \times \mathbf{E} \f$.
  */
  class UniformField : public ForceModifier {
  public:
    UniformField(SimInfo* info);

  private:
    void modifyForces() override;
    void initialize();

    bool initialized;
    bool doUniformField;
    bool doParticlePot;

    Globals* simParams;
    Vector3d EF;
  };
}  // namespace OpenMD

#endif
