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

#include "rnemd/NIVS.hpp"

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
#include "utils/Constants.hpp"

#define HONKING_LARGE_VALUE 1.0e10

namespace OpenMD::RNEMD {

  NIVSMethod::NIVSMethod(SimInfo* info, ForceManager* forceMan) :
      RNEMD {info, forceMan} {
    rnemdMethodLabel_ = "NIVS";

    RNEMDParameters* rnemdParams = info->getSimParams()->getRNEMDParameters();

    bool hasKineticFlux  = rnemdParams->haveKineticFlux();
    bool hasMomentumFlux = rnemdParams->haveMomentumFlux();

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
    case rnemdKePx:
    case rnemdKePy:
      hasCorrectFlux = hasMomentumFlux && hasKineticFlux;
      break;
    default:
      methodFluxMismatch = true;
      break;
    }

    if (methodFluxMismatch) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "RNEMD: The current method,\n"
               "\t\tNIVS\n"
               "\tcannot be used with the current flux type, %s\n",
               rnemdFluxTypeLabel_.c_str());
      painCave.isFatal  = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    if (!hasCorrectFlux) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "RNEMD: The current method, NIVS, and flux type, %s,\n"
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

  void NIVSMethod::doRNEMDImpl(SelectionManager& smanA,
                               SelectionManager& smanB) {
    if (!doRNEMD_) return;
    int selei;
    int selej;

    StuntDouble* sd;

    std::vector<StuntDouble*> hotBin, coldBin;

    RealType Phx = 0.0;
    RealType Phy = 0.0;
    RealType Phz = 0.0;
    RealType Khx = 0.0;
    RealType Khy = 0.0;
    RealType Khz = 0.0;
    RealType Khw = 0.0;
    RealType Pcx = 0.0;
    RealType Pcy = 0.0;
    RealType Pcz = 0.0;
    RealType Kcx = 0.0;
    RealType Kcy = 0.0;
    RealType Kcz = 0.0;
    RealType Kcw = 0.0;

    for (sd = smanA.beginSelected(selei); sd != NULL;
         sd = smanA.nextSelected(selei)) {
      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:

      if (usePeriodicBoundaryConditions_) currentSnap_->wrapVector(pos);

      RealType mass = sd->getMass();
      Vector3d vel  = sd->getVel();

      hotBin.push_back(sd);
      Phx += mass * vel.x();
      Phy += mass * vel.y();
      Phz += mass * vel.z();
      Khx += mass * vel.x() * vel.x();
      Khy += mass * vel.y() * vel.y();
      Khz += mass * vel.z() * vel.z();
      if (sd->isDirectional()) {
        Vector3d angMom = sd->getJ();
        Mat3x3d I       = sd->getI();
        if (sd->isLinear()) {
          int i = sd->linearAxis();
          int j = (i + 1) % 3;
          int k = (i + 2) % 3;
          Khw +=
              angMom[j] * angMom[j] / I(j, j) + angMom[k] * angMom[k] / I(k, k);
        } else {
          Khw += angMom[0] * angMom[0] / I(0, 0) +
                 angMom[1] * angMom[1] / I(1, 1) +
                 angMom[2] * angMom[2] / I(2, 2);
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

      coldBin.push_back(sd);
      Pcx += mass * vel.x();
      Pcy += mass * vel.y();
      Pcz += mass * vel.z();
      Kcx += mass * vel.x() * vel.x();
      Kcy += mass * vel.y() * vel.y();
      Kcz += mass * vel.z() * vel.z();
      if (sd->isDirectional()) {
        Vector3d angMom = sd->getJ();
        Mat3x3d I       = sd->getI();
        if (sd->isLinear()) {
          int i = sd->linearAxis();
          int j = (i + 1) % 3;
          int k = (i + 2) % 3;
          Kcw +=
              angMom[j] * angMom[j] / I(j, j) + angMom[k] * angMom[k] / I(k, k);
        } else {
          Kcw += angMom[0] * angMom[0] / I(0, 0) +
                 angMom[1] * angMom[1] / I(1, 1) +
                 angMom[2] * angMom[2] / I(2, 2);
        }
      }
    }

    Khx *= 0.5;
    Khy *= 0.5;
    Khz *= 0.5;
    Khw *= 0.5;
    Kcx *= 0.5;
    Kcy *= 0.5;
    Kcz *= 0.5;
    Kcw *= 0.5;

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &Phx, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Phy, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Phz, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Pcx, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Pcy, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Pcz, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(MPI_IN_PLACE, &Khx, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Khy, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Khz, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Khw, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(MPI_IN_PLACE, &Kcx, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Kcy, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Kcz, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Kcw, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
#endif

    // solve coldBin coeff's first
    RealType px = Pcx / Phx;
    RealType py = Pcy / Phy;
    RealType pz = Pcz / Phz;
    RealType c(0.0), x(0.0), y(0.0), z(0.0);
    bool successfulScale = false;
    if ((rnemdFluxType_ == rnemdFullKE) || (rnemdFluxType_ == rnemdRotKE)) {
      // may need sanity check Khw & Kcw > 0

      if (rnemdFluxType_ == rnemdFullKE) {
        c = 1.0 - kineticTarget_ / (Kcx + Kcy + Kcz + Kcw);
      } else {
        c = 1.0 - kineticTarget_ / Kcw;
      }

      if ((c > 0.81) && (c < 1.21)) {  // restrict scaling coefficients
        c = sqrt(c);

        RealType w = 0.0;
        if (rnemdFluxType_ == rnemdFullKE) {
          x = 1.0 + px * (1.0 - c);
          y = 1.0 + py * (1.0 - c);
          z = 1.0 + pz * (1.0 - c);
          /* more complicated way
        w = 1.0 + (Kcw - Kcw * c * c - (c * c * (Kcx + Kcy + Kcz
        + Khx * px * px + Khy * py * py + Khz * pz * pz)
        - 2.0 * c * (Khx * px * (1.0 + px) + Khy * py * (1.0 + py)
        + Khz * pz * (1.0 + pz)) + Khx * px * (2.0 + px)
        + Khy * py * (2.0 + py) + Khz * pz * (2.0 + pz)
        - Kcx - Kcy - Kcz)) / Khw; the following is simpler
      */
          if ((std::fabs(x - 1.0) < 0.1) && (std::fabs(y - 1.0) < 0.1) &&
              (std::fabs(z - 1.0) < 0.1)) {
            w = 1.0 + (kineticTarget_ + Khx * (1.0 - x * x) +
                       Khy * (1.0 - y * y) + Khz * (1.0 - z * z)) /
                          Khw;
          }  // no need to calculate w if x, y or z is out of range
        } else {
          w = 1.0 + kineticTarget_ / Khw;
        }
        if ((w > 0.81) && (w < 1.21)) {  // restrict scaling coefficients
          // if w is in the right range, so should be x, y, z.
          std::vector<StuntDouble*>::iterator sdi;
          Vector3d vel;
          for (sdi = coldBin.begin(); sdi != coldBin.end(); ++sdi) {
            if (rnemdFluxType_ == rnemdFullKE) {
              vel = (*sdi)->getVel() * c;
              (*sdi)->setVel(vel);
            }
            if ((*sdi)->isDirectional()) {
              Vector3d angMom = (*sdi)->getJ() * c;
              (*sdi)->setJ(angMom);
            }
          }
          w = sqrt(w);
          for (sdi = hotBin.begin(); sdi != hotBin.end(); ++sdi) {
            if (rnemdFluxType_ == rnemdFullKE) {
              vel = (*sdi)->getVel();
              vel.x() *= x;
              vel.y() *= y;
              vel.z() *= z;
              (*sdi)->setVel(vel);
            }
            if ((*sdi)->isDirectional()) {
              Vector3d angMom = (*sdi)->getJ() * w;
              (*sdi)->setJ(angMom);
            }
          }
          successfulScale = true;
          kineticExchange_ += kineticTarget_;
        }
      }
    } else {
      RealType a000(0.0), a110(0.0), c0(0.0);
      RealType a001(0.0), a111(0.0), b01(0.0), b11(0.0), c1(0.0);
      switch (rnemdFluxType_) {
      case rnemdKE:
        /* used hotBin coeff's & only scale x & y dimensions
        RealType px = Phx / Pcx;
        RealType py = Phy / Pcy;
        a110 = Khy;
        c0 = - Khx - Khy - kineticTarget_;
        a000 = Khx;
        a111 = Kcy * py * py;
        b11 = -2.0 * Kcy * py * (1.0 + py);
        c1 = Kcy * py * (2.0 + py) + Kcx * px * ( 2.0 + px) + kineticTarget_;
        b01 = -2.0 * Kcx * px * (1.0 + px);
        a001 = Kcx * px * px;
      */
        // scale all three dimensions, let c_x = c_y
        a000 = Kcx + Kcy;
        a110 = Kcz;
        c0   = kineticTarget_ - Kcx - Kcy - Kcz;
        a001 = Khx * px * px + Khy * py * py;
        a111 = Khz * pz * pz;
        b01  = -2.0 * (Khx * px * (1.0 + px) + Khy * py * (1.0 + py));
        b11  = -2.0 * Khz * pz * (1.0 + pz);
        c1   = Khx * px * (2.0 + px) + Khy * py * (2.0 + py) +
             Khz * pz * (2.0 + pz) - kineticTarget_;
        break;
      case rnemdPx:
        c    = 1 - momentumTarget_.x() / Pcx;
        a000 = Kcy;
        a110 = Kcz;
        c0   = Kcx * c * c - Kcx - Kcy - Kcz;
        a001 = py * py * Khy;
        a111 = pz * pz * Khz;
        b01  = -2.0 * Khy * py * (1.0 + py);
        b11  = -2.0 * Khz * pz * (1.0 + pz);
        c1   = Khy * py * (2.0 + py) + Khz * pz * (2.0 + pz) +
             Khx * (fastpow(c * px - px - 1.0, 2) - 1.0);
        break;
      case rnemdPy:
        c    = 1 - momentumTarget_.y() / Pcy;
        a000 = Kcx;
        a110 = Kcz;
        c0   = Kcy * c * c - Kcx - Kcy - Kcz;
        a001 = px * px * Khx;
        a111 = pz * pz * Khz;
        b01  = -2.0 * Khx * px * (1.0 + px);
        b11  = -2.0 * Khz * pz * (1.0 + pz);
        c1   = Khx * px * (2.0 + px) + Khz * pz * (2.0 + pz) +
             Khy * (fastpow(c * py - py - 1.0, 2) - 1.0);
        break;
      case rnemdPz:  // we don't really do this, do we?
        c    = 1 - momentumTarget_.z() / Pcz;
        a000 = Kcx;
        a110 = Kcy;
        c0   = Kcz * c * c - Kcx - Kcy - Kcz;
        a001 = px * px * Khx;
        a111 = py * py * Khy;
        b01  = -2.0 * Khx * px * (1.0 + px);
        b11  = -2.0 * Khy * py * (1.0 + py);
        c1   = Khx * px * (2.0 + px) + Khy * py * (2.0 + py) +
             Khz * (fastpow(c * pz - pz - 1.0, 2) - 1.0);
        break;
      default:
        break;
      }

      RealType v1  = a000 * a111 - a001 * a110;
      RealType v2  = a000 * b01;
      RealType v3  = a000 * b11;
      RealType v4  = a000 * c1 - a001 * c0;
      RealType v8  = a110 * b01;
      RealType v10 = -b01 * c0;

      RealType u0 = v2 * v10 - v4 * v4;
      RealType u1 = -2.0 * v3 * v4;
      RealType u2 = -v2 * v8 - v3 * v3 - 2.0 * v1 * v4;
      RealType u3 = -2.0 * v1 * v3;
      RealType u4 = -v1 * v1;
      // rescale coefficients
      RealType maxAbs = fabs(u0);
      if (maxAbs < fabs(u1)) maxAbs = fabs(u1);
      if (maxAbs < fabs(u2)) maxAbs = fabs(u2);
      if (maxAbs < fabs(u3)) maxAbs = fabs(u3);
      if (maxAbs < fabs(u4)) maxAbs = fabs(u4);
      u0 /= maxAbs;
      u1 /= maxAbs;
      u2 /= maxAbs;
      u3 /= maxAbs;
      u4 /= maxAbs;
      // max_element(start, end) is also available.
      Polynomial<RealType> poly;  // same as DoublePolynomial poly;
      poly.setCoefficient(4, u4);
      poly.setCoefficient(3, u3);
      poly.setCoefficient(2, u2);
      poly.setCoefficient(1, u1);
      poly.setCoefficient(0, u0);
      vector<RealType> realRoots = poly.FindRealRoots();

      vector<RealType>::iterator ri;
      RealType r1, r2, alpha0;
      vector<pair<RealType, RealType>> rps;
      for (ri = realRoots.begin(); ri != realRoots.end(); ++ri) {
        r2 = *ri;
        // Check to see if FindRealRoots() gave the right answer:
        if (fabs(u0 + r2 * (u1 + r2 * (u2 + r2 * (u3 + r2 * u4)))) > 1e-6) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "RNEMD Warning: polynomial solve seems to have an error!");
          painCave.isFatal = 0;
          simError();
          failRootCount_++;
        }
        // Might not be useful w/o rescaling coefficients
        alpha0 = -c0 - a110 * r2 * r2;
        if (alpha0 >= 0.0) {
          r1 = sqrt(alpha0 / a000);
          if (fabs(c1 + r1 * (b01 + r1 * a001) + r2 * (b11 + r2 * a111)) <
              1e-6) {
            rps.push_back(make_pair(r1, r2));
          }
          if (r1 > 1e-6) {  // r1 non-negative
            r1 = -r1;
            if (fabs(c1 + r1 * (b01 + r1 * a001) + r2 * (b11 + r2 * a111)) <
                1e-6) {
              rps.push_back(make_pair(r1, r2));
            }
          }
        }
      }
      // Consider combining together the part for solving for the pair
      // w/ the searching for the best solution part so that we don't
      // need the pairs vector:
      if (!rps.empty()) {
        RealType smallestDiff = HONKING_LARGE_VALUE;
        RealType diff(0.0);
        std::pair<RealType, RealType> bestPair = std::make_pair(1.0, 1.0);
        std::vector<std::pair<RealType, RealType>>::iterator rpi;
        for (rpi = rps.begin(); rpi != rps.end(); ++rpi) {
          r1 = (*rpi).first;
          r2 = (*rpi).second;
          switch (rnemdFluxType_) {
          case rnemdKE:
            diff = fastpow(1.0 - r1, 2) + fastpow(1.0 - r2, 2) +
                   fastpow(r1 * r1 / r2 / r2 - Kcz / Kcx, 2) +
                   fastpow(r1 * r1 / r2 / r2 - Kcz / Kcy, 2);
            break;
          case rnemdPx:
            diff = fastpow(1.0 - r1, 2) + fastpow(1.0 - r2, 2) +
                   fastpow(r1 * r1 / r2 / r2 - Kcz / Kcy, 2);
            break;
          case rnemdPy:
            diff = fastpow(1.0 - r1, 2) + fastpow(1.0 - r2, 2) +
                   fastpow(r1 * r1 / r2 / r2 - Kcz / Kcx, 2);
            break;
          case rnemdPz:
            diff = fastpow(1.0 - r1, 2) + fastpow(1.0 - r2, 2) +
                   fastpow(r1 * r1 / r2 / r2 - Kcy / Kcx, 2);
          default:
            break;
          }
          if (diff < smallestDiff) {
            smallestDiff = diff;
            bestPair     = *rpi;
          }
        }
#ifdef IS_MPI
        if (worldRank == 0) {
#endif
          // snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          //         "RNEMD: roots r1= %lf\tr2 = %lf\n",
          //         bestPair.first, bestPair.second);
          // painCave.isFatal = 0;
          // painCave.severity = OPENMD_INFO;
          // simError();
#ifdef IS_MPI
        }
#endif

        switch (rnemdFluxType_) {
        case rnemdKE:
          x = bestPair.first;
          y = bestPair.first;
          z = bestPair.second;
          break;
        case rnemdPx:
          x = c;
          y = bestPair.first;
          z = bestPair.second;
          break;
        case rnemdPy:
          x = bestPair.first;
          y = c;
          z = bestPair.second;
          break;
        case rnemdPz:
          x = bestPair.first;
          y = bestPair.second;
          z = c;
          break;
        default:
          break;
        }
        vector<StuntDouble*>::iterator sdi;
        Vector3d vel;
        for (sdi = coldBin.begin(); sdi != coldBin.end(); ++sdi) {
          vel = (*sdi)->getVel();
          vel.x() *= x;
          vel.y() *= y;
          vel.z() *= z;
          (*sdi)->setVel(vel);
        }
        // convert to hotBin coefficient
        x = 1.0 + px * (1.0 - x);
        y = 1.0 + py * (1.0 - y);
        z = 1.0 + pz * (1.0 - z);
        for (sdi = hotBin.begin(); sdi != hotBin.end(); ++sdi) {
          vel = (*sdi)->getVel();
          vel.x() *= x;
          vel.y() *= y;
          vel.z() *= z;
          (*sdi)->setVel(vel);
        }
        successfulScale = true;
        switch (rnemdFluxType_) {
        case rnemdKE:
          kineticExchange_ += kineticTarget_;
          break;
        case rnemdPx:
        case rnemdPy:
        case rnemdPz:
          momentumExchange_ += momentumTarget_;
          break;
        default:
          break;
        }
      }
    }
    if (successfulScale != true) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "NIVS exchange NOT performed - roots that solve\n"
               "\tthe constraint equations may not exist or there may be\n"
               "\tno selected objects in one or both slabs.\n");
      painCave.isFatal  = 0;
      painCave.severity = OPENMD_INFO;
      simError();
      failTrialCount_++;
    }
  }
}  // namespace OpenMD::RNEMD
