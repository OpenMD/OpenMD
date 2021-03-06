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

#ifndef INTEGRATOR_LDFORCEMANAGER_HPP
#define INTEGRATOR_LDFORCEMANAGER_HPP

#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "brains/ForceManager.hpp"
#include "brains/Velocitizer.hpp"
#include "hydrodynamics/Shape.hpp"
#include "primitives/Molecule.hpp"
#include "utils/RandNumGen.hpp"

namespace OpenMD {

  struct SDShape {
    StuntDouble* sd;
    Shape* shape;
  };

  struct MomentData {
    Vector3d rcr;   /**< Distance between CoM and Center of Resistance */
    Mat3x3d Icr;    /**< Moment of Inertia at Center of Resistance */
    Mat3x3d IcrInv; /**< Icr^{-1}  */
  };

  /**
   * @class LDForceManager
   * Force manager for Lagevin Dynamics applying friction and random
   * forces as well as torques.
   */
  class LDForceManager : public ForceManager {
  public:
    LDForceManager(SimInfo* info);

    int getMaxIterationNumber() { return maxIterNum_; }

    void setMaxIterationNumber(int maxIter) { maxIterNum_ = maxIter; }

    RealType getForceTolerance() { return forceTolerance_; }

    void setForceTolerance(RealType tol) { forceTolerance_ = tol; }

    RealType getDt2() { return dt2_; }

    void setDt2(RealType dt2) { dt2_ = dt2; }

  protected:
    virtual void postCalculation();

  private:
    std::map<std::string, HydroProp*> parseFrictionFile(
        const std::string& filename);
    MomentData* getMomentData(StuntDouble* sd);

    void genRandomForceAndTorque(Vector3d& force, Vector3d& torque,
                                 unsigned int index);

    std::map<std::string, HydroProp*> hydroPropMap_;
    std::vector<HydroProp*> hydroProps_;

    std::map<std::string, MomentData*> momentsMap_;
    std::vector<MomentData*> moments_;

    // convergence parameters:
    int maxIterNum_;
    RealType forceTolerance_;
    RealType dt2_;

    // random number generation:
    Utils::RandNumGenPtr randNumGen_;
    std::normal_distribution<RealType> forceDistribution_;

    RealType langevinBufferRadius_;
    RealType frozenBufferRadius_;
    bool sphericalBoundaryConditions_;
    Globals* simParams;
    std::unique_ptr<Velocitizer> veloMunge {nullptr};
  };
}  // namespace OpenMD
#endif  // BRAINS_FORCEMANAGER_HPP
