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

#ifndef INTEGRATORS_LAGNEVINPISTON_HPP
#define INTEGRATORS_LAGNEVINPISTON_HPP

#include <memory>
#include <random>

#include "integrators/NPT.hpp"
#include "utils/RandNumGen.hpp"

namespace OpenMD {

  /**
   * @class LangevinPiston
   * Constant pressure and temperature integrator
   *
   * The Langevin Piston Nosé-Hoover method in OpenMD combines the
   * Nosé-Hoover constant pressure method as described in
   * G.J. Martyna, D.J. Tobias and M.L. Klein, "Constant pressure
   * molecular dynamics algorithms", J. Chem. Phys. 101, 4177 (1994);
   * https://doi.org/10.1063/1.467468 , with piston fluctuation
   * control implemented using Langevin dynamics as in S.E. Feller,
   * Y. Zhang, R.W. Pastor and B.R. Brooks, "Constant pressure
   * molecular dynamics simulation: The Langevin piston method",
   * J. Chem. Phys. 103, 4613 (1995); https://doi.org/10.1063/1.470648
   */
  class LangevinPiston : public NPT {
  public:
    LangevinPiston(SimInfo* info);

  private:
    /* We need to implement moveA and moveB separately to leave out
     * the extended system thermostatting. Temperature control is
     * supplied in Langevin methods by connecting friction and random
     * forces via the second fluctuation dissipation theorem.
     */

    virtual void moveA();
    virtual void moveB();

    virtual void evolveEtaA();
    virtual void evolveEtaB();

    virtual bool etaConverged();

    virtual void getVelScaleA(Vector3d& sc, const Vector3d& vel);
    virtual void getVelScaleB(Vector3d& sc, int index);
    virtual void getPosScale(const Vector3d& pos, const Vector3d& COM,
                             int index, Vector3d& sc);

    virtual void calcVelScale();

    virtual void scaleSimBox();
    virtual RealType calcConservedQuantity() { return 0.0; }

    virtual void loadEta();
    virtual void saveEta();
    void genRandomForce(RealType& randomForce);

    RealType eta;
    RealType oldEta;
    RealType prevEta;
    RealType vScale;

    Utils::RandNumGenPtr randNumGen_;
    std::normal_distribution<RealType> forceDistribution_;

    RealType W_;
    RealType gamma_;
    RealType randomForce_;
  };
}  // namespace OpenMD

#endif  // INTEGRATORS_LANGEVINPISTON_HPP
