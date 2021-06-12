/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
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
   * J. Chem. Phys. 103, 4613 (1995); https://doi.org/10.1063/1.47064
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
