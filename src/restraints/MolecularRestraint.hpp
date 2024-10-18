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

#ifndef RESTRAINTS_MOLECULARRESTRAINT_HPP
#define RESTRAINTS_MOLECULARRESTRAINT_HPP

#include <vector>

#include "math/Vector3.hpp"
#include "restraints/Restraint.hpp"

namespace OpenMD {
  /**
   * @class MolecularRestraint
   *
   * MolecularRestraint is the restraint (both positional and
   * orientational) for the configuration of a flexible Molecule
   * relative to some reference structure for the same Molecule.  The
   * angles that define the deflection away from the reference
   * structure are the Euler angles taken from the rotation matrix
   * that gives the lowest root mean square deviation (RMSD), while
   * the displacement of the molecule is simply the displacement of
   * the center of mass relative to the reference structure.
   */
  class MolecularRestraint : public Restraint {
  public:
    MolecularRestraint() : Restraint() {}

    void setReferenceStructure(std::vector<Vector3d> ref, Vector3d refCom) {
      ref_    = ref;
      refCom_ = refCom;

      std::vector<Vector3d>::iterator i;

      for (i = ref_.begin(); i != ref_.end(); ++i) {
        (*i) = (*i) - refCom_;
      }

      forces_.clear();
      forces_.resize(ref_.size());
    }

    void calcForce(std::vector<Vector3d> struc, Vector3d molCom);

    std::vector<Vector3d> getRestraintForces() { return forces_; }

  private:
    std::vector<Vector3d> ref_;
    std::vector<Vector3d> forces_;
    Vector3d refCom_;
  };
}  // namespace OpenMD

#endif
