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

#ifndef OPENMD_PERTURBATIONS_VELOCITYFIELDPARAMETERS_HPP
#define OPENMD_PERTURBATIONS_VELOCITYFIELDPARAMETERS_HPP

#include <vector>

#include "types/DataHolder.hpp"
#include "utils/ParameterManager.hpp"

namespace OpenMD::Perturbations {

  /*! \class VelocityFieldParameters
      \brief Parsed contents of the velocityField{ ... } block.

      Specifies a linear (homogeneous) ambient velocity field,
      \f$ \mathbf{v}(\mathbf{r}) = \mathbf{v}_0 + \mathsf{K}\cdot\mathbf{r}\f$,
      via its physically meaningful parts.  The user may supply a
      constant-stress (rate-of-strain) block, a constant-vorticity vector,
      or both; at least one is required.

      The rate-of-strain block follows the UniformGradient idiom: two
      direction vectors and a scalar rate build the symmetric, traceless
      tensor
      \f$ E_{ij} = \dot\varepsilon\left[\tfrac12(a_i b_j + a_j b_i)
                    - \tfrac13(\mathbf{a}\cdot\mathbf{b})\,\delta_{ij}\right]. \f$

      The vorticity vector \f$\boldsymbol\omega\f$ builds the antisymmetric
      spin tensor \f$ W = \tfrac12[\boldsymbol\omega]_\times \f$ so that the
      rotational part of the flow is \f$\tfrac12\,\boldsymbol\omega\times
      \mathbf{r}\f$ and \f$\nabla\times\mathbf{v}=\boldsymbol\omega\f$.

      Units (OpenMD internal): strainRate and vorticity in fs^-1,
      backgroundVelocity in Angstrom fs^-1.
  */
  class VelocityFieldParameters : public DataHolder {
    DeclareParameter(UseVelocityField, bool);
    DeclareParameter(StrainRate, RealType);
    DeclareParameter(StrainDirection1, std::vector<RealType>);
    DeclareParameter(StrainDirection2, std::vector<RealType>);
    DeclareParameter(Vorticity, std::vector<RealType>);
    DeclareParameter(BackgroundVelocity, std::vector<RealType>);

  public:
    VelocityFieldParameters();
    virtual ~VelocityFieldParameters() = default;
    virtual void validate();
  };
}  // namespace OpenMD::Perturbations

#endif  // OPENMD_PERTURBATIONS_VELOCITYFIELDPARAMETERS_HPP
