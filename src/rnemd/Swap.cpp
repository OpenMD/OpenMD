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

#include "rnemd/Swap.hpp"

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
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/Accumulator.hpp"
#include "utils/Constants.hpp"
#include "utils/Tuple.hpp"

#define HONKING_LARGE_VALUE 1.0e10

namespace OpenMD {
  namespace RNEMD {

    SwapMethod::SwapMethod(SimInfo* info, ForceManager* forceMan) :
        RNEMD {info, forceMan} {
      rnemdMethodLabel_ = "Swap";

      RNEMDParameters* rnemdParams = info->getSimParams()->getRNEMDParameters();

      bool hasKineticFlux  = rnemdParams->haveKineticFlux();
      bool hasMomentumFlux = rnemdParams->haveMomentumFlux();

      bool methodFluxMismatch = false;
      bool hasCorrectFlux     = false;

      switch (rnemdFluxType_) {
      case rnemdKE:
        hasCorrectFlux = hasKineticFlux;
        break;
      case rnemdPx:
      case rnemdPy:
      case rnemdPz:
        hasCorrectFlux = hasMomentumFlux;
        break;
      default:
        methodFluxMismatch = true;
        break;
      }

      if (methodFluxMismatch) {
        sprintf(painCave.errMsg,
                "RNEMD: The current method,\n"
                "\t\tSwap\n"
                "\tcannot be used with the current flux type, %s\n",
                rnemdFluxTypeLabel_.c_str());
        painCave.isFatal  = 1;
        painCave.severity = OPENMD_ERROR;
        simError();
      }

      if (!hasCorrectFlux) {
        sprintf(painCave.errMsg,
                "RNEMD: The current method, Swap, and flux type, %s,\n"
                "\tdid not have the correct flux value specified. Options\n"
                "\tinclude: kineticFlux and momentumFlux.\n",
                rnemdFluxTypeLabel_.c_str());
        painCave.isFatal  = 1;
        painCave.severity = OPENMD_ERROR;
        simError();
      }

      if (hasKineticFlux) {
        setKineticFlux(rnemdParams->getKineticFlux());
      } else {
        setKineticFlux(0.0);
      }

      if (hasMomentumFlux) {
        RealType momentumFlux = rnemdParams->getMomentumFlux();
        std::vector<RealType> momentumFluxVector(3);

        switch (rnemdFluxType_) {
        case rnemdPx:
          momentumFluxVector[0] = momentumFlux;
          break;
        case rnemdPy:
          momentumFluxVector[1] = momentumFlux;
          break;
        case rnemdPz:
          momentumFluxVector[2] = momentumFlux;
          break;
        default:
          break;
        }

        setMomentumFluxVector(momentumFluxVector);
      }
    }

    void SwapMethod::doRNEMDImpl(SelectionManager& smanA,
                                 SelectionManager& smanB) {
      if (!doRNEMD_) return;
      int selei;
      int selej;

      StuntDouble* sd;

      RealType min_val(0.0);
      int min_found       = 0;
      StuntDouble* min_sd = NULL;

      RealType max_val(0.0);
      int max_found       = 0;
      StuntDouble* max_sd = NULL;

      for (sd = smanA.beginSelected(selei); sd != NULL;
           sd = smanA.nextSelected(selei)) {
        Vector3d pos = sd->getPos();

        // wrap the stuntdouble's position back into the box:

        if (usePeriodicBoundaryConditions_) currentSnap_->wrapVector(pos);

        RealType mass = sd->getMass();
        Vector3d vel  = sd->getVel();
        RealType value(0.0);

        switch (rnemdFluxType_) {
        case rnemdKE:

          value = mass * vel.lengthSquare();

          if (sd->isDirectional()) {
            Vector3d angMom = sd->getJ();
            Mat3x3d I       = sd->getI();

            if (sd->isLinear()) {
              int i = sd->linearAxis();
              int j = (i + 1) % 3;
              int k = (i + 2) % 3;
              value += angMom[j] * angMom[j] / I(j, j) +
                       angMom[k] * angMom[k] / I(k, k);
            } else {
              value += angMom[0] * angMom[0] / I(0, 0) +
                       angMom[1] * angMom[1] / I(1, 1) +
                       angMom[2] * angMom[2] / I(2, 2);
            }
          }  // angular momenta exchange enabled
          value *= 0.5;
          break;
        case rnemdPx:
          value = mass * vel[0];
          break;
        case rnemdPy:
          value = mass * vel[1];
          break;
        case rnemdPz:
          value = mass * vel[2];
          break;
        default:
          break;
        }
        if (!max_found) {
          max_val   = value;
          max_sd    = sd;
          max_found = 1;
        } else {
          if (max_val < value) {
            max_val = value;
            max_sd  = sd;
          }
        }
      }

      for (sd = smanB.beginSelected(selej); sd != NULL;
           sd = smanB.nextSelected(selej)) {
        Vector3d pos = sd->getPos();

        // wrap the stuntdouble's position back into the box:

        if (usePeriodicBoundaryConditions_) currentSnap_->wrapVector(pos);

        RealType mass = sd->getMass();
        Vector3d vel  = sd->getVel();
        RealType value(0.0);

        switch (rnemdFluxType_) {
        case rnemdKE:

          value = mass * vel.lengthSquare();

          if (sd->isDirectional()) {
            Vector3d angMom = sd->getJ();
            Mat3x3d I       = sd->getI();

            if (sd->isLinear()) {
              int i = sd->linearAxis();
              int j = (i + 1) % 3;
              int k = (i + 2) % 3;
              value += angMom[j] * angMom[j] / I(j, j) +
                       angMom[k] * angMom[k] / I(k, k);
            } else {
              value += angMom[0] * angMom[0] / I(0, 0) +
                       angMom[1] * angMom[1] / I(1, 1) +
                       angMom[2] * angMom[2] / I(2, 2);
            }
          }  // angular momenta exchange enabled
          value *= 0.5;
          break;
        case rnemdPx:
          value = mass * vel[0];
          break;
        case rnemdPy:
          value = mass * vel[1];
          break;
        case rnemdPz:
          value = mass * vel[2];
          break;
        default:
          break;
        }

        if (!min_found) {
          min_val   = value;
          min_sd    = sd;
          min_found = 1;
        } else {
          if (min_val > value) {
            min_val = value;
            min_sd  = sd;
          }
        }
      }

#ifdef IS_MPI
      int worldRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

      int my_min_found = min_found;
      int my_max_found = max_found;

      // Even if we didn't find a minimum, did someone else?
      MPI_Allreduce(&my_min_found, &min_found, 1, MPI_INT, MPI_LOR,
                    MPI_COMM_WORLD);
      // Even if we didn't find a maximum, did someone else?
      MPI_Allreduce(&my_max_found, &max_found, 1, MPI_INT, MPI_LOR,
                    MPI_COMM_WORLD);
#endif

      if (max_found && min_found) {
#ifdef IS_MPI
        struct {
          RealType val;
          int rank;
        } max_vals, min_vals;

        if (my_min_found) {
          min_vals.val = min_val;
        } else {
          min_vals.val = HONKING_LARGE_VALUE;
        }
        min_vals.rank = worldRank;

        // Who had the minimum?
        MPI_Allreduce(&min_vals, &min_vals, 1, MPI_REALTYPE_INT, MPI_MINLOC,
                      MPI_COMM_WORLD);
        min_val = min_vals.val;

        if (my_max_found) {
          max_vals.val = max_val;
        } else {
          max_vals.val = -HONKING_LARGE_VALUE;
        }
        max_vals.rank = worldRank;

        // Who had the maximum?
        MPI_Allreduce(&max_vals, &max_vals, 1, MPI_REALTYPE_INT, MPI_MAXLOC,
                      MPI_COMM_WORLD);
        max_val = max_vals.val;
#endif

        if (min_val < max_val) {
#ifdef IS_MPI
          if (max_vals.rank == worldRank && min_vals.rank == worldRank) {
            // I have both maximum and minimum, so proceed like a single
            // processor version:
#endif

            Vector3d min_vel = min_sd->getVel();
            Vector3d max_vel = max_sd->getVel();
            RealType temp_vel;

            switch (rnemdFluxType_) {
            case rnemdKE:
              min_sd->setVel(max_vel);
              max_sd->setVel(min_vel);
              if (min_sd->isDirectional() && max_sd->isDirectional()) {
                Vector3d min_angMom = min_sd->getJ();
                Vector3d max_angMom = max_sd->getJ();
                min_sd->setJ(max_angMom);
                max_sd->setJ(min_angMom);
              }  // angular momenta exchange enabled
                 // assumes same rigid body identity
              break;
            case rnemdPx:
              temp_vel    = min_vel.x();
              min_vel.x() = max_vel.x();
              max_vel.x() = temp_vel;
              min_sd->setVel(min_vel);
              max_sd->setVel(max_vel);
              break;
            case rnemdPy:
              temp_vel    = min_vel.y();
              min_vel.y() = max_vel.y();
              max_vel.y() = temp_vel;
              min_sd->setVel(min_vel);
              max_sd->setVel(max_vel);
              break;
            case rnemdPz:
              temp_vel    = min_vel.z();
              min_vel.z() = max_vel.z();
              max_vel.z() = temp_vel;
              min_sd->setVel(min_vel);
              max_sd->setVel(max_vel);
              break;
            default:
              break;
            }

#ifdef IS_MPI
            // the rest of the cases only apply in parallel simulations:
          } else if (max_vals.rank == worldRank) {
            // I had the max, but not the minimum

            Vector3d min_vel;
            Vector3d max_vel = max_sd->getVel();
            MPI_Status status;

            // point-to-point swap of the velocity vector
            MPI_Sendrecv(max_vel.getArrayPointer(), 3, MPI_REALTYPE,
                         min_vals.rank, 0, min_vel.getArrayPointer(), 3,
                         MPI_REALTYPE, min_vals.rank, 0, MPI_COMM_WORLD,
                         &status);

            switch (rnemdFluxType_) {
            case rnemdKE:
              max_sd->setVel(min_vel);
              // angular momenta exchange enabled
              if (max_sd->isDirectional()) {
                Vector3d min_angMom;
                Vector3d max_angMom = max_sd->getJ();

                // point-to-point swap of the angular momentum vector
                MPI_Sendrecv(max_angMom.getArrayPointer(), 3, MPI_REALTYPE,
                             min_vals.rank, 1, min_angMom.getArrayPointer(), 3,
                             MPI_REALTYPE, min_vals.rank, 1, MPI_COMM_WORLD,
                             &status);

                max_sd->setJ(min_angMom);
              }
              break;
            case rnemdPx:
              max_vel.x() = min_vel.x();
              max_sd->setVel(max_vel);
              break;
            case rnemdPy:
              max_vel.y() = min_vel.y();
              max_sd->setVel(max_vel);
              break;
            case rnemdPz:
              max_vel.z() = min_vel.z();
              max_sd->setVel(max_vel);
              break;
            default:
              break;
            }
          } else if (min_vals.rank == worldRank) {
            // I had the minimum but not the maximum:

            Vector3d max_vel;
            Vector3d min_vel = min_sd->getVel();
            MPI_Status status;

            // point-to-point swap of the velocity vector
            MPI_Sendrecv(min_vel.getArrayPointer(), 3, MPI_REALTYPE,
                         max_vals.rank, 0, max_vel.getArrayPointer(), 3,
                         MPI_REALTYPE, max_vals.rank, 0, MPI_COMM_WORLD,
                         &status);

            switch (rnemdFluxType_) {
            case rnemdKE:
              min_sd->setVel(max_vel);
              // angular momenta exchange enabled
              if (min_sd->isDirectional()) {
                Vector3d min_angMom = min_sd->getJ();
                Vector3d max_angMom;

                // point-to-point swap of the angular momentum vector
                MPI_Sendrecv(min_angMom.getArrayPointer(), 3, MPI_REALTYPE,
                             max_vals.rank, 1, max_angMom.getArrayPointer(), 3,
                             MPI_REALTYPE, max_vals.rank, 1, MPI_COMM_WORLD,
                             &status);

                min_sd->setJ(max_angMom);
              }
              break;
            case rnemdPx:
              min_vel.x() = max_vel.x();
              min_sd->setVel(min_vel);
              break;
            case rnemdPy:
              min_vel.y() = max_vel.y();
              min_sd->setVel(min_vel);
              break;
            case rnemdPz:
              min_vel.z() = max_vel.z();
              min_sd->setVel(min_vel);
              break;
            default:
              break;
            }
          }
#endif

          switch (rnemdFluxType_) {
          case rnemdKE:
            kineticExchange_ += max_val - min_val;
            break;
          case rnemdPx:
            momentumExchange_.x() += max_val - min_val;
            break;
          case rnemdPy:
            momentumExchange_.y() += max_val - min_val;
            break;
          case rnemdPz:
            momentumExchange_.z() += max_val - min_val;
            break;
          default:
            break;
          }
        } else {
          sprintf(painCave.errMsg, "RNEMD::doSwap exchange NOT performed "
                                   "because min_val > max_val\n");
          painCave.isFatal  = 0;
          painCave.severity = OPENMD_INFO;
          simError();
          failTrialCount_++;
        }
      } else {
        sprintf(painCave.errMsg,
                "Swap exchange NOT performed because selected object\n"
                "\twas not present in at least one of the two slabs.\n");
        painCave.isFatal  = 0;
        painCave.severity = OPENMD_INFO;
        simError();
        failTrialCount_++;
      }
    }
  }  // namespace RNEMD
}  // namespace OpenMD
