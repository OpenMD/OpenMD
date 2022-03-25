/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

#include "rnemd/VSS.hpp"

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

namespace OpenMD::RNEMD {

  VSSMethod::VSSMethod(SimInfo* info, ForceManager* forceMan) :
      RNEMD {info, forceMan} {
    rnemdMethodLabel_ = "VSS";

    RNEMDParameters* rnemdParams = info->getSimParams()->getRNEMDParameters();

    bool hasKineticFlux         = rnemdParams->haveKineticFlux();
    bool hasMomentumFlux        = rnemdParams->haveMomentumFlux();
    bool hasMomentumFluxVector  = rnemdParams->haveMomentumFluxVector();
    bool hasAngularMomentumFlux = rnemdParams->haveAngularMomentumFlux();
    bool hasAngularMomentumFluxVector =
        rnemdParams->haveAngularMomentumFluxVector();

    bool methodFluxMismatch = false;
    bool hasCorrectFlux     = false;

    switch (rnemdFluxType_) {
    case rnemdKE:
    case rnemdRotKE:
    case rnemdFullKE:
      hasCorrectFlux = hasKineticFlux;
      break;
    case rnemdPx:
    case rnemdPy:
    case rnemdPz:
      hasCorrectFlux = hasMomentumFlux;
      break;
    case rnemdLx:
    case rnemdLy:
    case rnemdLz:
      hasCorrectFlux = hasAngularMomentumFlux;
      break;
    case rnemdPvector:
      hasCorrectFlux = hasMomentumFluxVector;
      break;
    case rnemdLvector:
      hasCorrectFlux = hasAngularMomentumFluxVector;
      break;
    case rnemdKePx:
    case rnemdKePy:
      hasCorrectFlux = hasMomentumFlux && hasKineticFlux;
      break;
    case rnemdKeLx:
    case rnemdKeLy:
    case rnemdKeLz:
      hasCorrectFlux = hasAngularMomentumFlux && hasKineticFlux;
      break;
    case rnemdKePvector:
      hasCorrectFlux = hasMomentumFluxVector && hasKineticFlux;
      break;
    case rnemdKeLvector:
      hasCorrectFlux = hasAngularMomentumFluxVector && hasKineticFlux;
      break;
    default:
      methodFluxMismatch = true;
      break;
    }

    if (methodFluxMismatch) {
      sprintf(painCave.errMsg,
              "RNEMD: The current method,\n"
              "\t\tVSS\n"
              "\tcannot be used with the current flux type, %s\n",
              rnemdFluxTypeLabel_.c_str());
      painCave.isFatal  = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    if (!hasCorrectFlux) {
      sprintf(painCave.errMsg,
              "RNEMD: The current method, VSS, and flux type, %s,\n"
              "\tdid not have the correct flux value specified. Options\n"
              "\tinclude: kineticFlux, momentumFlux, angularMomentumFlux,\n"
              "\tmomentumFluxVector, and angularMomentumFluxVector.\n",
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

    if (hasMomentumFluxVector) {
      setMomentumFluxVector(rnemdParams->getMomentumFluxVector());
    } else {
      std::vector<RealType> momentumFluxVector(3);

      if (hasMomentumFlux) {
        RealType momentumFlux = rnemdParams->getMomentumFlux();

        switch (rnemdFluxType_) {
        case rnemdPx:
        case rnemdKePx:
          momentumFluxVector[0] = momentumFlux;
          break;
        case rnemdPy:
        case rnemdKePy:
          momentumFluxVector[1] = momentumFlux;
          break;
        case rnemdPz:
          momentumFluxVector[2] = momentumFlux;
          break;
        default:
          break;
        }
      }

      setMomentumFluxVector(momentumFluxVector);
    }

    if (hasAngularMomentumFluxVector) {
      setAngularMomentumFluxVector(rnemdParams->getAngularMomentumFluxVector());
    } else {
      std::vector<RealType> angularMomentumFluxVector(3);

      if (hasAngularMomentumFlux) {
        RealType angularMomentumFlux = rnemdParams->getAngularMomentumFlux();

        switch (rnemdFluxType_) {
        case rnemdLx:
        case rnemdKeLx:
          angularMomentumFluxVector[0] = angularMomentumFlux;
          break;
        case rnemdLy:
        case rnemdKeLy:
          angularMomentumFluxVector[1] = angularMomentumFlux;
          break;
        case rnemdLz:
        case rnemdKeLz:
          angularMomentumFluxVector[2] = angularMomentumFlux;
        default:
          break;
        }
      }

      setAngularMomentumFluxVector(angularMomentumFluxVector);
    }
  }

  void VSSMethod::doRNEMDImpl(SelectionManager& smanA,
                              SelectionManager& smanB) {
    if (!doRNEMD_) return;
    int selei;
    int selej;

    StuntDouble* sd;

    vector<StuntDouble*> hotBin, coldBin;

    Vector3d Ph(V3Zero);
    Vector3d Lh(V3Zero);
    RealType Mh = 0.0;
    Mat3x3d Ih(0.0);
    RealType Kh = 0.0;
    Vector3d Pc(V3Zero);
    Vector3d Lc(V3Zero);
    RealType Mc = 0.0;
    Mat3x3d Ic(0.0);
    RealType Kc = 0.0;

    // Constraints can be on only the linear or angular momentum, but
    // not both.  Usually, the user will specify which they want, but
    // in case they don't, the use of periodic boundaries should make
    // the choice for us.
    bool doLinearPart  = false;
    bool doAngularPart = false;

    switch (rnemdFluxType_) {
    case rnemdPx:
    case rnemdPy:
    case rnemdPz:
    case rnemdPvector:
    case rnemdKePx:
    case rnemdKePy:
    case rnemdKePvector:
      doLinearPart = true;
      break;
    case rnemdLx:
    case rnemdLy:
    case rnemdLz:
    case rnemdLvector:
    case rnemdKeLx:
    case rnemdKeLy:
    case rnemdKeLz:
    case rnemdKeLvector:
      doAngularPart = true;
      break;
    case rnemdKE:
    case rnemdRotKE:
    case rnemdFullKE:
    default:
      if (usePeriodicBoundaryConditions_)
        doLinearPart = true;
      else
        doAngularPart = true;
      break;
    }

    for (sd = smanA.beginSelected(selei); sd != NULL;
         sd = smanA.nextSelected(selei)) {
      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:
      if (usePeriodicBoundaryConditions_) currentSnap_->wrapVector(pos);

      RealType mass = sd->getMass();
      Vector3d vel  = sd->getVel();
      Vector3d rPos = sd->getPos() - coordinateOrigin_;
      RealType r2;

      hotBin.push_back(sd);
      Ph += mass * vel;
      Mh += mass;
      Kh += mass * vel.lengthSquare();
      Lh += mass * cross(rPos, vel);
      Ih -= outProduct(rPos, rPos) * mass;
      r2 = rPos.lengthSquare();
      Ih(0, 0) += mass * r2;
      Ih(1, 1) += mass * r2;
      Ih(2, 2) += mass * r2;

      if (rnemdFluxType_ == rnemdFullKE) {
        if (sd->isDirectional()) {
          Vector3d angMom = sd->getJ();
          Mat3x3d I       = sd->getI();
          if (sd->isLinear()) {
            int i = sd->linearAxis();
            int j = (i + 1) % 3;
            int k = (i + 2) % 3;
            Kh += angMom[j] * angMom[j] / I(j, j) +
                  angMom[k] * angMom[k] / I(k, k);
          } else {
            Kh += angMom[0] * angMom[0] / I(0, 0) +
                  angMom[1] * angMom[1] / I(1, 1) +
                  angMom[2] * angMom[2] / I(2, 2);
          }
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
      Vector3d rPos = sd->getPos() - coordinateOrigin_;
      RealType r2;

      coldBin.push_back(sd);
      Pc += mass * vel;
      Mc += mass;
      Kc += mass * vel.lengthSquare();
      Lc += mass * cross(rPos, vel);
      Ic -= outProduct(rPos, rPos) * mass;
      r2 = rPos.lengthSquare();
      Ic(0, 0) += mass * r2;
      Ic(1, 1) += mass * r2;
      Ic(2, 2) += mass * r2;

      if (rnemdFluxType_ == rnemdFullKE) {
        if (sd->isDirectional()) {
          Vector3d angMom = sd->getJ();
          Mat3x3d I       = sd->getI();
          if (sd->isLinear()) {
            int i = sd->linearAxis();
            int j = (i + 1) % 3;
            int k = (i + 2) % 3;
            Kc += angMom[j] * angMom[j] / I(j, j) +
                  angMom[k] * angMom[k] / I(k, k);
          } else {
            Kc += angMom[0] * angMom[0] / I(0, 0) +
                  angMom[1] * angMom[1] / I(1, 1) +
                  angMom[2] * angMom[2] / I(2, 2);
          }
        }
      }
    }

    Kh *= 0.5;
    Kc *= 0.5;

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &Ph[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Pc[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Lh[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Lc[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Mh, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Kh, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Mc, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Kc, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, Ih.getArrayPointer(), 9, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, Ic.getArrayPointer(), 9, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
#endif

    Vector3d ac, acrec, bc, bcrec;
    Vector3d ah, ahrec, bh, bhrec;

    bool successfulExchange = false;
    if ((Mh > 0.0) && (Mc > 0.0)) {  // both slabs are not empty

      Vector3d vc = Pc / Mc;
      ac          = -momentumTarget_ / Mc + vc;
      acrec       = -momentumTarget_ / Mc;

      // We now need the inverse of the inertia tensor to calculate the
      // angular velocity of the cold slab;
      Mat3x3d Ici     = Ic.inverse();
      Vector3d omegac = Ici * Lc;
      bc              = -(Ici * angularMomentumTarget_) + omegac;
      bcrec           = bc - omegac;

      RealType cNumerator = Kc - kineticTarget_;
      if (doLinearPart) cNumerator -= 0.5 * Mc * ac.lengthSquare();

      if (doAngularPart) cNumerator -= 0.5 * (dot(bc, Ic * bc));

      RealType cDenominator = Kc;

      if (doLinearPart) cDenominator -= 0.5 * Mc * vc.lengthSquare();

      if (doAngularPart) cDenominator -= 0.5 * (dot(omegac, Ic * omegac));

      if (cNumerator / cDenominator > 0.0) {
        RealType c = sqrt(cNumerator / cDenominator);

        if ((c > 0.9) && (c < 1.1)) {  // restrict scaling coefficients

          Vector3d vh = Ph / Mh;
          ah          = momentumTarget_ / Mh + vh;
          ahrec       = momentumTarget_ / Mh;

          // We now need the inverse of the inertia tensor to
          // calculate the angular velocity of the hot slab;
          Mat3x3d Ihi     = Ih.inverse();
          Vector3d omegah = Ihi * Lh;
          bh              = (Ihi * angularMomentumTarget_) + omegah;
          bhrec           = bh - omegah;

          RealType hNumerator = Kh + kineticTarget_;
          if (doLinearPart) hNumerator -= 0.5 * Mh * ah.lengthSquare();

          if (doAngularPart) hNumerator -= 0.5 * (dot(bh, Ih * bh));

          RealType hDenominator = Kh;
          if (doLinearPart) hDenominator -= 0.5 * Mh * vh.lengthSquare();
          if (doAngularPart) hDenominator -= 0.5 * (dot(omegah, Ih * omegah));

          if (hNumerator / hDenominator > 0.0) {
            RealType h = sqrt(hNumerator / hDenominator);

            if ((h > 0.9) && (h < 1.1)) {
              vector<StuntDouble*>::iterator sdi;
              Vector3d vel;
              Vector3d rPos;

              for (sdi = coldBin.begin(); sdi != coldBin.end(); ++sdi) {
                if (doLinearPart) vel = ((*sdi)->getVel() - vc) * c + ac;
                if (doAngularPart) {
                  rPos = (*sdi)->getPos() - coordinateOrigin_;
                  vel  = ((*sdi)->getVel() - cross(omegac, rPos)) * c +
                        cross(bc, rPos);
                }

                (*sdi)->setVel(vel);

                if (rnemdFluxType_ == rnemdFullKE) {
                  if ((*sdi)->isDirectional()) {
                    Vector3d angMom = (*sdi)->getJ() * c;
                    (*sdi)->setJ(angMom);
                  }
                }
              }

              for (sdi = hotBin.begin(); sdi != hotBin.end(); ++sdi) {
                if (doLinearPart) vel = ((*sdi)->getVel() - vh) * h + ah;
                if (doAngularPart) {
                  rPos = (*sdi)->getPos() - coordinateOrigin_;
                  vel  = ((*sdi)->getVel() - cross(omegah, rPos)) * h +
                        cross(bh, rPos);
                }

                (*sdi)->setVel(vel);

                if (rnemdFluxType_ == rnemdFullKE) {
                  if ((*sdi)->isDirectional()) {
                    Vector3d angMom = (*sdi)->getJ() * h;
                    (*sdi)->setJ(angMom);
                  }
                }
              }

              successfulExchange = true;
              kineticExchange_ += kineticTarget_;
              momentumExchange_ += momentumTarget_;
              angularMomentumExchange_ += angularMomentumTarget_;
            }
          }
        }
      }
    }

    if (successfulExchange != true) {
      sprintf(painCave.errMsg,
              "VSS exchange NOT performed - roots that solve\n"
              "\tthe constraint equations may not exist or there may be\n"
              "\tno selected objects in one or both slabs.\n");
      painCave.isFatal  = 0;
      painCave.severity = OPENMD_INFO;
      simError();
      failTrialCount_++;
    }
  }
}  // namespace OpenMD::RNEMD
