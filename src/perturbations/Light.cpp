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

#include "perturbations/Light.hpp"

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "brains/ForceModifier.hpp"
#include "nonbonded/NonBondedInteraction.hpp"
#include "primitives/Molecule.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "types/MultipoleAdapter.hpp"
#include "utils/Constants.hpp"

namespace OpenMD::Perturbations {
  Light::Light(SimInfo* info) :
      ForceModifier {info}, initialized {false}, doLight {false},
      doParticlePot {false}, info_(info) {
    lightParams = info_->getSimParams()->getLightParameters();
  }

  void Light::initialize() {
    bool haveE0           = false;
    bool haveDirection    = false;
    bool haveFrequency    = false;
    bool havePolarization = false;

    if (lightParams->haveWaveVector()) {
      std::vector<RealType> k = lightParams->getWaveVector();
      // wave vectors are input in inverse angstroms, so no unit conversion:
      k_.x()        = k[0];
      k_.y()        = k[1];
      k_.z()        = k[2];
      kmag_         = k_.length();
      lambda_       = 2.0 * Constants::PI / kmag_;
      omega_        = 2.0 * Constants::PI * Constants::c / lambda_;
      haveFrequency = true;
      khat_         = k_;
      khat_.normalize();
      haveDirection = true;
    }

    if (lightParams->havePropagationDirection()) {
      if (haveDirection) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "light: please specify either waveVector or "
                 "propagationDirection, but not both.\n");
        painCave.isFatal = 1;
        simError();
      }
      std::vector<RealType> pd = lightParams->getPropagationDirection();
      khat_.x()                = pd[0];
      khat_.y()                = pd[1];
      khat_.z()                = pd[2];
      khat_.normalize();
      haveDirection = true;
    }

    if (lightParams->haveWavelength()) {
      if (haveFrequency) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "light: please specify one of: waveVector, wavelength, or"
                 "frequency (but only one of these).\n");
        painCave.isFatal = 1;
        simError();
      }
      // wavelengths are entered in nm to work with experimentalists.
      // Convert to angstroms:
      lambda_       = lightParams->getWavelength() * 10.0;
      omega_        = 2.0 * Constants::PI * Constants::c / lambda_;
      kmag_         = 2.0 * Constants::PI / lambda_;
      haveFrequency = true;
    }

    if (lightParams->haveFrequency()) {
      if (haveFrequency) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "light: please specify one of: waveVector, wavelength, or"
                 "frequency (but only one of these).\n");
        painCave.isFatal = 1;
        simError();
      }
      // frequencies are entered in Hz to work with experimentalists.
      // Convert to fs^-1
      omega_        = lightParams->getFrequency() * 1.0e-15;
      lambda_       = 2.0 * Constants::PI * Constants::c / omega_;
      kmag_         = 2.0 * Constants::PI / lambda_;
      haveFrequency = true;
    }

    if (haveFrequency && haveDirection) { k_ = khat_ * kmag_; }

    if (lightParams->haveIntensity()) {
      RealType intense = lightParams->getIntensity();
      // intensities are input in W/cm^2
      intense *= 1.439326e-11;
      E0_ = std::sqrt(2.0 * intense / (Constants::epsilon0 * Constants::c));
      // E0 now has units of kcal/mol e^-1 Angstroms^-1
      haveE0 = true;
    }

    // Determine Polarization Type
    jones_.clear();
    jones_.reserve(2);
    std::map<std::string, LightPolarization> stringToPolarization;

    stringToPolarization["X"] = lightX;
    stringToPolarization["Y"] = lightY;
    stringToPolarization["+"] = lightPlus;
    stringToPolarization["-"] = lightMinus;

    if (lightParams->havePolarization()) {
      std::string lpl      = lightParams->getPolarization();
      LightPolarization lp = stringToPolarization.find(lpl)->second;
      switch (lp) {
      case lightX:
        jones_[0]        = {1.0, 0.0};
        jones_[1]        = {0.0, 0.0};
        havePolarization = true;
        break;
      case lightY:
        jones_[0]        = {0.0, 0.0};
        jones_[1]        = {1.0, 0.0};
        havePolarization = true;
        break;
      case lightPlus:
        jones_[0]        = {1.0, 0.0};
        jones_[1]        = {0.0, 1.0};
        havePolarization = true;
        break;
      case lightMinus:
        jones_[0]        = {1.0, 0.0};
        jones_[1]        = {0.0, -1.0};
        havePolarization = true;
        break;
      default:
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Light: Unknown polarization type\n");
        painCave.isFatal  = 1;
        painCave.severity = OPENMD_ERROR;
        simError();
        break;
      }

    } else {
      std::string allowedPolarizations;
      int currentLineLength = 0;

      for (std::map<std::string, LightPolarization>::iterator polStrIter =
               stringToPolarization.begin();
           polStrIter != stringToPolarization.end(); ++polStrIter) {
        allowedPolarizations += polStrIter->first + ", ";
        currentLineLength += polStrIter->first.length() + 2;

        if (currentLineLength >= 50) {
          allowedPolarizations += "\n\t\t";
          currentLineLength = 0;
        }
      }

      allowedPolarizations.erase(allowedPolarizations.length() - 2, 2);

      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "Light: No polarization was set in the omd file. This parameter\n"
          "\tmust be set to use Light, and can take any of these values:\n"
          "\t\t%s.\n",
          allowedPolarizations.c_str());
      painCave.isFatal  = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    if (haveE0 && haveDirection && haveFrequency && havePolarization) {
      doLight = true;
    } else {
      if (!haveDirection) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Light: could not determine direction of propagation.\n");
        painCave.isFatal = 1;
        simError();
      }
      if (!haveE0) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Light: intensity not specified.\n");
        painCave.isFatal = 1;
        simError();
      }
      if (!haveFrequency) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Light: could not determine frequency or wavelength.\n");
        painCave.isFatal = 1;
        simError();
      }
      if (!havePolarization) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Light: polarization  not specifieid.\n");
        painCave.isFatal = 1;
        simError();
      }
    }

    int storageLayout_ = info_->getSnapshotManager()->getAtomStorageLayout();
    if (storageLayout_ & DataStorage::dslParticlePot) doParticlePot = true;

    // Relatively simple Euler angles between khat_ and lab frame:

    RealType psi = 0.0;
    RealType theta =
        acos(std::min((RealType)1.0, std::max((RealType)-1.0, khat_[2])));
    RealType phi = std::atan2(-khat_[1], khat_[0]);

    if (phi < 0) phi += 2.0 * Constants::PI;

    A_.setupRotMat(phi, theta, psi);
    Ainv_ = A_.inverse();

    initialized = true;
  }

  void Light::modifyForces() {
    if (!initialized) initialize();

    SimInfo::MoleculeIterator i;
    Molecule::AtomIterator j;
    Molecule* mol;
    Atom* atom;
    AtomType* atype;
    potVec longRangePotential(0.0);
    int l, m, n;
    RealType C {}, U {}, fPot {};
    Vector3d r {}, rp {}, v {}, f {}, trq {}, D {}, J {}, av {};
    Vector3d EFk {}, EF {}, BF {};
    Mat3x3d I {};

    bool isCharge;

    if (doLight) {
      RealType t = info_->getSnapshotManager()->getCurrentSnapshot()->getTime();
      U          = 0.0;
      fPot       = 0.0;

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (atom = mol->beginAtom(j); atom != NULL; atom = mol->nextAtom(j)) {
          isCharge = false;
          C        = 0.0;

          atype = atom->getAtomType();

          // We are not wrapping coordinates for light interactions:
          r = atom->getPos();
          v = atom->getVel();

          rp = A_ * r;  // atom's position in frame of light propagation

          // e^{ i (k*z - omega * t) } is the main oscillatory component:
          std::complex<RealType> argument(0.0, kmag_ * rp.z() - omega_ * t);
          std::complex<RealType> Ex = E0_ * jones_[0] * std::exp(argument);
          std::complex<RealType> Ey = E0_ * jones_[1] * std::exp(argument);

          EFk.x() = Ex.real();
          EFk.y() = Ey.real();
          EFk.z() = 0.0;

          EF = Ainv_ * EFk;  // electric field rotated back into lab coordinates

          // The magnetic field (BF) is perpendicular to both electric field
          // and light propagation direction:

          BF = cross(EF, khat_) / Constants::c;

          atom->addElectricField(EF);

          FixedChargeAdapter fca = FixedChargeAdapter(atype);
          if (fca.isFixedCharge()) {
            isCharge = true;
            C        = fca.getCharge();
          }

          FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
          if (fqa.isFluctuatingCharge()) {
            isCharge = true;
            C += atom->getFlucQPos();
            atom->addFlucQFrc(dot(r, EF));
          }

          if (isCharge) {
            f = C * (EF + cross(v, BF));
            atom->addFrc(f);
            U = -dot(r, f);
            if (doParticlePot) { atom->addParticlePot(U); }
            fPot += U;
          }

          MultipoleAdapter ma = MultipoleAdapter(atype);
          if (ma.isDipole()) {
            D = atom->getDipole() * Constants::dipoleConvert;

            trq += cross(D, EF) + cross(D, cross(v, BF));

            atom->addTrq(trq);

            J = atom->getJ();
            I = atom->getI();
            if (atom->isLinear()) {
              l     = atom->linearAxis();
              m     = (l + 1) % 3;
              n     = (l + 2) % 3;
              av[l] = 0;
              av[m] = J[m] / I(m, m);
              av[n] = J[n] / I(n, n);
            } else {
              av = I.inverse() * J;
            }

            f = cross(cross(av, D), BF);
            atom->addFrc(f);

            U = -dot(D, EF);
            if (doParticlePot) { atom->addParticlePot(U); }
            fPot += U;
          }
        }
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, &fPot, 1, MPI_REALTYPE, MPI_SUM,
                    MPI_COMM_WORLD);
#endif

      Snapshot* snap     = info_->getSnapshotManager()->getCurrentSnapshot();
      longRangePotential = snap->getLongRangePotentials();
      longRangePotential[ELECTROSTATIC_FAMILY] += fPot;
      snap->setLongRangePotentials(longRangePotential);
    }
  }
}  // namespace OpenMD::Perturbations
