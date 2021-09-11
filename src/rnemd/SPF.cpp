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

#include "rnemd/SPF.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "brains/ForceManager.hpp"
#include "brains/Thermo.hpp"
#include "io/Globals.hpp"
#include "math/ConvexHull.hpp"
#include "math/Polynomial.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector.hpp"
#include "math/Vector3.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/StuntDouble.hpp"
#include "rnemd/RNEMD.hpp"
#include "rnemd/RNEMDParameters.hpp"
#include "rnemd/SPFForceManager.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/Accumulator.hpp"
#include "utils/Constants.hpp"

namespace OpenMD::RNEMD {

  SPFMethod::SPFMethod(SimInfo* info, ForceManager* forceMan) :
      RNEMD {info, forceMan}, sourceSman_ {info} {
    rnemdMethodLabel_ = "SPF";

    if (SPFForceManager* spfForceManager =
            dynamic_cast<SPFForceManager*>(forceMan)) {
      forceManager_ = spfForceManager;
    } else {
      sprintf(painCave.errMsg,
              "SPF-RNEMD cannot be used with the default ForceManager.\n");
      painCave.isFatal  = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    RNEMDParameters* rnemdParams = info->getSimParams()->getRNEMDParameters();

    bool hasParticleFlux = rnemdParams->haveParticleFlux();

    bool methodFluxMismatch = false;
    bool hasCorrectFlux     = false;

    switch (rnemdFluxType_) {
    case rnemdParticle:
      hasCorrectFlux = hasParticleFlux;
      break;
    default:
      methodFluxMismatch = true;
      break;
    }

    if (methodFluxMismatch) {
      sprintf(painCave.errMsg,
              "RNEMD: The current method, SPF\n"
              "\tcannot be used with the current flux type, %s\n",
              rnemdFluxTypeLabel_.c_str());
      painCave.isFatal  = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    if (!hasCorrectFlux) {
      sprintf(painCave.errMsg,
              "RNEMD: The current method, SPF, and flux type, %s,\n"
              "\tdid not have the correct flux value specified. Options\n"
              "\tinclude: particleFlux.\n",
              rnemdFluxTypeLabel_.c_str());
      painCave.isFatal  = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    if (hasParticleFlux) {
      setParticleFlux(rnemdParams->getParticleFlux());
    } else {
      setParticleFlux(0.0);
    }

    selectMolecule();
  }

  void SPFMethod::doRNEMDImpl(SelectionManager& smanA,
                              SelectionManager& smanB) {
    // First kick will succeed
    if (!hasMovingMolecule_) { selectMolecule(); }

    if (!doRNEMD_) return;
    int selei;
    int selej;

    StuntDouble* sd;

    std::vector<StuntDouble*> binA, binB;

    Vector3d P_a {V3Zero};
    Vector3d P_b {V3Zero};
    RealType M_a {};
    RealType M_b {};
    RealType K_a {};
    RealType K_b {};

    for (sd = smanA.beginSelected(selei); sd != NULL;
         sd = smanA.nextSelected(selei)) {
      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:
      if (usePeriodicBoundaryConditions_) currentSnap_->wrapVector(pos);

      RealType mass = sd->getMass();
      Vector3d vel  = sd->getVel();

      binA.push_back(sd);

      P_a += mass * vel;
      M_a += mass;
      K_a += mass * vel.lengthSquare();
    }

    for (sd = smanB.beginSelected(selej); sd != NULL;
         sd = smanB.nextSelected(selej)) {
      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:
      if (usePeriodicBoundaryConditions_) currentSnap_->wrapVector(pos);

      RealType mass = sd->getMass();
      Vector3d vel  = sd->getVel();

      binB.push_back(sd);

      P_b += mass * vel;
      M_b += mass;
      K_b += mass * vel.lengthSquare();
    }

    K_a *= 0.5;
    K_b *= 0.5;

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &P_a[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &P_b[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &M_a, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &M_b, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &K_a, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &K_b, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
#endif

    bool successfulExchange = false;

    if ((M_a > 0.0) && (M_b > 0.0) &&
        hasMovingMolecule_) {  // both slabs are not empty
      Vector3d v_a = P_a / M_a;
      Vector3d v_b = P_b / M_b;

      RealType numerator = forceManager_->getDeltaU() * particleTarget_ *
                           Constants::energyConvert;

      RealType denominator = K_a + K_b;
      denominator -= 0.5 * M_a * v_a.lengthSquare();
      denominator -= 0.5 * M_b * v_b.lengthSquare();

      RealType a2 = (numerator / denominator) + 1.0;

      if (a2 > 0.0) {
        RealType a = std::sqrt(a2);

        if ((a > 0.9) && (a < 1.1)) {  // restrict scaling coefficients
          std::vector<StuntDouble*>::iterator sdi;
          Vector3d vel;

          for (sdi = binA.begin(); sdi != binA.end(); ++sdi) {
            vel = ((*sdi)->getVel() - v_a) * a;

            (*sdi)->setVel(vel);
          }

          for (sdi = binB.begin(); sdi != binB.end(); ++sdi) {
            vel = ((*sdi)->getVel() - v_b) * a;

            (*sdi)->setVel(vel);
          }

          successfulExchange = true;
          particleExchange_ += particleTarget_;
          updateLambda();
        }
      }
    }

    if (successfulExchange != true) {
      sprintf(painCave.errMsg,
              "SPF exchange NOT performed - roots that solve\n"
              "\tthe constraint equations may not exist or there may be\n"
              "\tno selected objects in one or both slabs.\n");
      painCave.isFatal  = 0;
      painCave.severity = OPENMD_INFO;
      simError();
      failTrialCount_++;
    }

    // First kick will fail
    // if (!hasMovingMolecule_) { selectMolecule(); }
  }

  void SPFMethod::updateLambda() {
    lambda_ += std::fabs(particleTarget_);

    if (lambda_ > 1.0) {
      lambda_ -= 1.0;
      forceManager_->updateLambda(lambda_);
      forceManager_->acceptGhost();

      selectMolecule();
    } else if (std::fabs(lambda_ - 1.0) < 1e-6) {
      lambda_ = 0.0;
      forceManager_->updateLambda(lambda_);
      forceManager_->acceptGhost();

      selectMolecule();
    } else {
      forceManager_->updateLambda(lambda_);
    }
  }

  void SPFMethod::selectMolecule() {
    // The sign of our flux determines which slab is the source and which is
    // the sink
    if (particleTarget_ > 0.0) {
      sourceSman_       = commonA_;
      targetSlabCenter_ = slabBCenter_;
    } else {
      sourceSman_       = commonB_;
      targetSlabCenter_ = slabACenter_;
    }

    if (sourceSman_.getMoleculeSelectionCount() == 0) {
      hasMovingMolecule_ = false;
      forceManager_->setSelectedMolecule(NULL, V3Zero);
      return;
    }

    int whichSelectedID {};
    Molecule* selectedMolecule;

    Utils::RandNumGenPtr randNumGen = info_->getRandomNumberGenerator();

#ifdef IS_MPI
    int worldRank {};
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    if (worldRank == 0) {
#endif
      std::uniform_int_distribution<> selectedMoleculeDistribution {
          0, sourceSman_.getMoleculeSelectionCount() - 1};

      whichSelectedID = selectedMoleculeDistribution(*randNumGen);
#ifdef IS_MPI
    }

    MPI_Bcast(&whichSelectedID, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    selectedMolecule = sourceSman_.nthSelectedMolecule(whichSelectedID);

    if (selectedMolecule) {
      int axis0 = (rnemdPrivilegedAxis_ + 1) % 3;
      int axis1 = (rnemdPrivilegedAxis_ + 2) % 3;
      int axis2 = rnemdPrivilegedAxis_;

      std::uniform_real_distribution<RealType> distr0 {0, hmat_(axis0, axis0)};
      std::uniform_real_distribution<RealType> distr1 {0, hmat_(axis1, axis1)};

      std::normal_distribution<RealType> distr2 {targetSlabCenter_,
                                                 0.25 * slabWidth_};

      Vector3d newPosition {V3Zero};

      newPosition[axis0] = distr0(*randNumGen);
      newPosition[axis1] = distr1(*randNumGen);
      newPosition[axis2] = distr2(*randNumGen);

      forceManager_->setSelectedMolecule(selectedMolecule, newPosition);
    }

    hasMovingMolecule_ = true;
  }
}  // namespace OpenMD::RNEMD