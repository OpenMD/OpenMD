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

#ifndef APPLICATIONS_DYNAMICPROPS_COLLECTIVEDIPOLEDISPLACEMENT_HPP
#define APPLICATIONS_DYNAMICPROPS_COLLECTIVEDIPOLEDISPLACEMENT_HPP

#include "applications/dynamicProps/TimeCorrFunc.hpp"
#include "brains/Thermo.hpp"

namespace OpenMD {
  //! Calculates the collective dipole displacement function
  /*! This time correlation function is the Helfand moment conjugate
      to the current density. Helfand moments are used to calculate
      the Einstein-Helfand relations for transport that are formally
      equivalent to Green-Kubo expressions using a related flux.  In
      this case, the flux,

      \f[ \mathbf{J}(t) = \sum_{i=1}^{N} q_i \mathbf{v}_{\mathrm{cm},i}(t) \f]

      is normally used to calculate an ionic conductivity,

      \f[ \sigma = \frac{1}{3V k_b T} \int_0^\infty \left< \mathbf{J}(0) \cdot
     \mathbf{J}(t) \right> dt \f]

      The cm subscript denotes center of mass locations for all molecules.

      This class computes the collective translational dipole moment,

      \f[ \mathbf{M}_\mathrm{trans}(t) = \sum_{i=1}^{N} q_i
     \mathbf{r}_{\mathrm{cm},i}(t) \f]

      as well as total contributions to the system's net dipole moment

      \f[ \mathbf{M}_\mathrm{tot}(t) =  \sum_{i=1}^{N} \sum_{a} q_{ia}
     \mathbf{r}_{ia}(t) = \sum_{i=1}^{N} q_i \mathbf{r}_{\mathrm{cq},i}(t) \f]

      where cq denotes the molecular center of charge.  It also
      calculates the rotational contribution,

      \f[ \mathbf{M}_\mathrm{rot}(t) = \sum_{i=1}^{N} q_i \left[
     \mathbf{r}_{\mathrm{cq},i}(t) - \mathbf{r}_{\mathrm{cm},i}(t) \right] \f]

      The correlation functions are the displacements of these terms
      from their values at an earlier time,

      \f[ \left< \left| \mathbf{M}_\mathrm{trans}(t) -
     \mathbf{M}_\mathrm{trans}(0) \right|^2 \right> \f]

      and identical quantities for the total and rotational contributions.
  */
  class CollectiveDipoleDisplacement : public SystemACF<Vector3d> {
  public:
    CollectiveDipoleDisplacement(SimInfo* info, const std::string& filename,
                                 const std::string& sele1,
                                 const std::string& sele2);

  private:
    virtual void computeProperty1(int frame);
    virtual Vector3d calcCorrVal(int frame1, int frame2);

    Thermo* thermo_;

    std::vector<Vector3d> CRcm_;
    std::vector<Vector3d> CRtot_;
    std::vector<Vector3d> CRrot_;
  };
}  // namespace OpenMD

#endif
