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

#include "brains/Thermo.hpp"

#include <cmath>
#include <iostream>

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "primitives/Molecule.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "types/MultipoleAdapter.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"
#ifdef HAVE_QHULL
#include "math/AlphaHull.hpp"
#include "math/ConvexHull.hpp"
#endif

using namespace std;
namespace OpenMD {

  RealType Thermo::getTranslationalKinetic() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasTranslationalKineticEnergy) {
      SimInfo::MoleculeIterator miter;
      vector<StuntDouble*>::iterator iiter;
      Molecule* mol;
      StuntDouble* sd;
      Vector3d vel;
      RealType mass;
      RealType kinetic(0.0);

      for (mol = info_->beginMolecule(miter); mol != NULL;
           mol = info_->nextMolecule(miter)) {
        for (sd = mol->beginIntegrableObject(iiter); sd != NULL;
             sd = mol->nextIntegrableObject(iiter)) {
          mass = sd->getMass();
          vel  = sd->getVel();

          kinetic +=
              mass * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
        }
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, &kinetic, 1, MPI_REALTYPE, MPI_SUM,
                    MPI_COMM_WORLD);
#endif

      kinetic = kinetic * 0.5 / Constants::energyConvert;

      snap->setTranslationalKineticEnergy(kinetic);
    }
    return snap->getTranslationalKineticEnergy();
  }

  RealType Thermo::getRotationalKinetic() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasRotationalKineticEnergy) {
      SimInfo::MoleculeIterator miter;
      vector<StuntDouble*>::iterator iiter;
      Molecule* mol;
      StuntDouble* sd;
      Vector3d angMom;
      Mat3x3d I;
      int i, j, k;
      RealType kinetic(0.0);

      for (mol = info_->beginMolecule(miter); mol != NULL;
           mol = info_->nextMolecule(miter)) {
        for (sd = mol->beginIntegrableObject(iiter); sd != NULL;
             sd = mol->nextIntegrableObject(iiter)) {
          if (sd->isDirectional()) {
            angMom = sd->getJ();
            I      = sd->getI();

            if (sd->isLinear()) {
              i = sd->linearAxis();
              j = (i + 1) % 3;
              k = (i + 2) % 3;
              kinetic += angMom[j] * angMom[j] / I(j, j) +
                         angMom[k] * angMom[k] / I(k, k);
            } else {
              kinetic += angMom[0] * angMom[0] / I(0, 0) +
                         angMom[1] * angMom[1] / I(1, 1) +
                         angMom[2] * angMom[2] / I(2, 2);
            }
          }
        }
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, &kinetic, 1, MPI_REALTYPE, MPI_SUM,
                    MPI_COMM_WORLD);
#endif

      kinetic = kinetic * 0.5 / Constants::energyConvert;

      snap->setRotationalKineticEnergy(kinetic);
    }
    return snap->getRotationalKineticEnergy();
  }

  RealType Thermo::getKinetic() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasKineticEnergy) {
      RealType ke = getTranslationalKinetic() + getRotationalKinetic() +
                    getElectronicKinetic();

      snap->setKineticEnergy(ke);
    }
    return snap->getKineticEnergy();
  }

  RealType Thermo::getPotential() {
    // ForceManager computes the potential and stores it in the
    // Snapshot.  All we have to do is report it.

    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    return snap->getPotentialEnergy();
  }

  potVec Thermo::getSelectionPotentials() {
    // ForceManager computes the selection potentials and stores them
    // in the Snapshot.  All we have to do is report them.

    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    return snap->getSelectionPotentials();
  }

  RealType Thermo::getTotalEnergy() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasTotalEnergy) {
      snap->setTotalEnergy(this->getKinetic() + this->getPotential());
    }

    return snap->getTotalEnergy();
  }

  /*
   * Returns only the nuclear portion of the temperature - see
   * getElectronicTemperature for the electronic portion */
  RealType Thermo::getTemperature() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasTemperature) {
      RealType nuclearKE =
          this->getTranslationalKinetic() + this->getRotationalKinetic();

      RealType temperature {};

      // With no degrees of freedom, T is not well defined, but we'll
      // set to zero.

      if (info_->getNdf() > 0)
        temperature = (2.0 * nuclearKE) / (info_->getNdf() * Constants::kb);
      else
        temperature = 0.0;

      snap->setTemperature(temperature);
    }

    return snap->getTemperature();
  }

  RealType Thermo::getElectronicKinetic() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasElectronicKineticEnergy) {
      SimInfo::MoleculeIterator miter;
      vector<Atom*>::iterator iiter;
      Molecule* mol;
      Atom* atom;
      RealType cvel;
      RealType cmass;
      RealType kinetic(0.0);

      for (mol = info_->beginMolecule(miter); mol != NULL;
           mol = info_->nextMolecule(miter)) {
        for (atom = mol->beginFluctuatingCharge(iiter); atom != NULL;
             atom = mol->nextFluctuatingCharge(iiter)) {
          cmass = atom->getChargeMass();
          cvel  = atom->getFlucQVel();
          kinetic += cmass * cvel * cvel;
        }
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, &kinetic, 1, MPI_REALTYPE, MPI_SUM,
                    MPI_COMM_WORLD);
#endif

      kinetic *= 0.5;
      snap->setElectronicKineticEnergy(kinetic);
    }

    return snap->getElectronicKineticEnergy();
  }

  RealType Thermo::getElectronicTemperature() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasElectronicTemperature) {
      RealType eTemp = (2.0 * this->getElectronicKinetic()) /
                       (info_->getNFluctuatingCharges() * Constants::kb);

      snap->setElectronicTemperature(eTemp);
    }

    return snap->getElectronicTemperature();
  }

  RealType Thermo::getNetCharge() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasNetCharge) {
      SimInfo::MoleculeIterator miter;
      vector<Atom*>::iterator aiter;
      Molecule* mol;
      Atom* atom;
      RealType charge;
      RealType netCharge(0.0);

      for (mol = info_->beginMolecule(miter); mol != NULL;
           mol = info_->nextMolecule(miter)) {
        for (atom = mol->beginAtom(aiter); atom != NULL;
             atom = mol->nextAtom(aiter)) {
          charge = 0.0;

          FixedChargeAdapter fca = FixedChargeAdapter(atom->getAtomType());
          if (fca.isFixedCharge()) { charge = fca.getCharge(); }

          FluctuatingChargeAdapter fqa =
              FluctuatingChargeAdapter(atom->getAtomType());
          if (fqa.isFluctuatingCharge()) { charge += atom->getFlucQPos(); }

          netCharge += charge;
        }
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, &netCharge, 1, MPI_REALTYPE, MPI_SUM,
                    MPI_COMM_WORLD);
#endif

      snap->setNetCharge(netCharge);
    }

    return snap->getNetCharge();
  }

  RealType Thermo::getChargeMomentum() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasChargeMomentum) {
      SimInfo::MoleculeIterator miter;
      vector<Atom*>::iterator iiter;
      Molecule* mol;
      Atom* atom;
      RealType cvel;
      RealType cmass;
      RealType momentum(0.0);

      for (mol = info_->beginMolecule(miter); mol != NULL;
           mol = info_->nextMolecule(miter)) {
        for (atom = mol->beginFluctuatingCharge(iiter); atom != NULL;
             atom = mol->nextFluctuatingCharge(iiter)) {
          cmass = atom->getChargeMass();
          cvel  = atom->getFlucQVel();

          momentum += cmass * cvel;
        }
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, &momentum, 1, MPI_REALTYPE, MPI_SUM,
                    MPI_COMM_WORLD);
#endif

      snap->setChargeMomentum(momentum);
    }

    return snap->getChargeMomentum();
  }

  std::vector<Vector3d> Thermo::getCurrentDensity() {
    Snapshot* snap       = info_->getSnapshotManager()->getCurrentSnapshot();
    AtomTypeSet simTypes = info_->getSimulatedAtomTypes();

    SimInfo::MoleculeIterator miter;
    std::vector<Atom*>::iterator iiter;
    std::vector<RigidBody*>::iterator ri;
    Molecule* mol;
    RigidBody* rb;
    Atom* atom;
    AtomType* atype;
    AtomTypeSet::iterator at;
    Vector3d Jc(0.0);
    std::vector<Vector3d> typeJc(simTypes.size(), V3Zero);

    for (mol = info_->beginMolecule(miter); mol != NULL;
         mol = info_->nextMolecule(miter)) {
      // change the velocities of atoms which belong to the rigidbodies
      for (rb = mol->beginRigidBody(ri); rb != NULL;
           rb = mol->nextRigidBody(ri)) {
        rb->updateAtomVel();
      }

      for (atom = mol->beginAtom(iiter); atom != NULL;
           atom = mol->nextAtom(iiter)) {
        Vector3d v = atom->getVel();
        RealType q = 0.0;
        int typeIndex(-1);

        atype                  = atom->getAtomType();
        FixedChargeAdapter fca = FixedChargeAdapter(atype);
        if (fca.isFixedCharge()) q = fca.getCharge();
        FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
        if (fqa.isFluctuatingCharge()) q += atom->getFlucQPos();

        typeIndex = -1;
        at        = std::find(simTypes.begin(), simTypes.end(), atype);
        if (at != simTypes.end()) {
          typeIndex = std::distance(simTypes.begin(), at);
        }

        if (typeIndex != -1) { typeJc[typeIndex] += q * v; }
        Jc += q * v;
      }
    }

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &Jc[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    for (unsigned int j = 0; j < simTypes.size(); j++) {
      MPI_Allreduce(MPI_IN_PLACE, &typeJc[j][0], 3, MPI_REALTYPE, MPI_SUM,
                    MPI_COMM_WORLD);
    }
#endif

    RealType vol = snap->getVolume();

    Jc /= (vol * Constants::currentDensityConvert);
    for (unsigned int j = 0; j < simTypes.size(); j++) {
      typeJc[j] /= (vol * Constants::currentDensityConvert);
    }

    std::vector<Vector3d> result;
    result.clear();

    result.push_back(Jc);
    for (unsigned int j = 0; j < simTypes.size(); j++)
      result.push_back(typeJc[j]);

    return result;
  }

  RealType Thermo::getVolume() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    return snap->getVolume();
  }

  RealType Thermo::getPressure() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasPressure) {
      // Relies on the calculation of the full molecular pressure tensor

      Mat3x3d tensor;
      RealType pressure;

      tensor = getPressureTensor();

      pressure = Constants::pressureConvert *
                 (tensor(0, 0) + tensor(1, 1) + tensor(2, 2)) / 3.0;

      snap->setPressure(pressure);
    }

    return snap->getPressure();
  }

  RealType Thermo::getPressure(Snapshot* snap) {
    if (!snap->hasPressure) {
      // Relies on the calculation of the full molecular pressure tensor
      Mat3x3d tensor;
      RealType pressure;

      tensor = getPressureTensor(snap);

      pressure = Constants::pressureConvert *
                 (tensor(0, 0) + tensor(1, 1) + tensor(2, 2)) / 3.0;

      snap->setPressure(pressure);
    }

    return snap->getPressure();
  }

  Mat3x3d Thermo::getPressureTensor() {
    // returns pressure tensor in units amu*fs^-2*Ang^-1
    // routine derived via viral theorem description in:
    // Paci, E. and Marchi, M. J.Phys.Chem. 1996, 100, 4314-4322
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasPressureTensor) {
      Mat3x3d pressureTensor;
      Mat3x3d p_tens(0.0);
      RealType mass;
      Vector3d vcom;

      SimInfo::MoleculeIterator i;
      vector<StuntDouble*>::iterator j;
      Molecule* mol;
      StuntDouble* sd;
      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (sd = mol->beginIntegrableObject(j); sd != NULL;
             sd = mol->nextIntegrableObject(j)) {
          mass = sd->getMass();
          vcom = sd->getVel();
          p_tens += mass * outProduct(vcom, vcom);
        }
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, p_tens.getArrayPointer(), 9, MPI_REALTYPE,
                    MPI_SUM, MPI_COMM_WORLD);
#endif

      RealType volume      = this->getVolume();
      Mat3x3d virialTensor = snap->getVirialTensor();

      pressureTensor =
          (p_tens + Constants::energyConvert * virialTensor) / volume;

      snap->setPressureTensor(pressureTensor);
    }
    return snap->getPressureTensor();
  }

  Mat3x3d Thermo::getPressureTensor(Snapshot* snap) {
    // returns pressure tensor in units amu*fs^-2*Ang^-1
    // routine derived via viral theorem description in:
    // Paci, E. and Marchi, M. J.Phys.Chem. 1996, 100, 4314-4322
    if (!snap->hasPressureTensor) {
      Mat3x3d pressureTensor;
      Mat3x3d p_tens(0.0);
      RealType mass;
      Vector3d vcom;

      SimInfo::MoleculeIterator i;
      vector<StuntDouble*>::iterator j;
      Molecule* mol;
      StuntDouble* sd;
      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (sd = mol->beginIntegrableObject(j); sd != NULL;
             sd = mol->nextIntegrableObject(j)) {
          mass = sd->getMass();
          vcom = sd->getVel(snap);
          p_tens += mass * outProduct(vcom, vcom);
        }
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, p_tens.getArrayPointer(), 9, MPI_REALTYPE,
                    MPI_SUM, MPI_COMM_WORLD);
#endif

      RealType volume      = this->getVolume();
      Mat3x3d virialTensor = snap->getVirialTensor();

      pressureTensor =
          (p_tens + Constants::energyConvert * virialTensor) / volume;

      snap->setPressureTensor(pressureTensor);
    }
    return snap->getPressureTensor();
  }

  Vector3d Thermo::getSystemDipole() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasSystemDipole) {
      SimInfo::MoleculeIterator miter;
      vector<Atom*>::iterator aiter;
      Molecule* mol;
      Atom* atom;
      RealType charge;
      Vector3d ri(0.0);
      Vector3d dipoleVector(0.0);
      Vector3d nPos(0.0);
      Vector3d pPos(0.0);
      RealType nChg(0.0);
      RealType pChg(0.0);
      int nCount = 0;
      int pCount = 0;

      RealType chargeToC   = 1.60217733e-19;
      RealType angstromToM = 1.0e-10;
      RealType debyeToCm   = 3.33564095198e-30;

      for (mol = info_->beginMolecule(miter); mol != NULL;
           mol = info_->nextMolecule(miter)) {
        for (atom = mol->beginAtom(aiter); atom != NULL;
             atom = mol->nextAtom(aiter)) {
          charge = 0.0;

          FixedChargeAdapter fca = FixedChargeAdapter(atom->getAtomType());
          if (fca.isFixedCharge()) { charge = fca.getCharge(); }

          FluctuatingChargeAdapter fqa =
              FluctuatingChargeAdapter(atom->getAtomType());
          if (fqa.isFluctuatingCharge()) { charge += atom->getFlucQPos(); }

          charge *= chargeToC;

          ri = atom->getPos();
          snap->wrapVector(ri);
          ri *= angstromToM;

          if (charge < 0.0) {
            nPos += ri;
            nChg -= charge;
            nCount++;
          } else if (charge > 0.0) {
            pPos += ri;
            pChg += charge;
            pCount++;
          }

          if (atom->isDipole()) {
            dipoleVector += atom->getDipole() * debyeToCm;
          }
        }
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, &pChg, 1, MPI_REALTYPE, MPI_SUM,
                    MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &nChg, 1, MPI_REALTYPE, MPI_SUM,
                    MPI_COMM_WORLD);

      MPI_Allreduce(MPI_IN_PLACE, &pCount, 1, MPI_INTEGER, MPI_SUM,
                    MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &nCount, 1, MPI_INTEGER, MPI_SUM,
                    MPI_COMM_WORLD);

      MPI_Allreduce(MPI_IN_PLACE, pPos.getArrayPointer(), 3, MPI_REALTYPE,
                    MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, nPos.getArrayPointer(), 3, MPI_REALTYPE,
                    MPI_SUM, MPI_COMM_WORLD);

      MPI_Allreduce(MPI_IN_PLACE, dipoleVector.getArrayPointer(), 3,
                    MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
#endif

      // first load the accumulated dipole moment (if dipoles were present)
      Vector3d boxDipole = dipoleVector;
      // now include the dipole moment due to charges
      // use the lesser of the positive and negative charge totals
      RealType chg_value = nChg <= pChg ? nChg : pChg;

      // find the average positions
      if (pCount > 0 && nCount > 0) {
        pPos /= pCount;
        nPos /= nCount;
      }

      // dipole is from the negative to the positive (physics notation)
      boxDipole += (pPos - nPos) * chg_value;
      snap->setSystemDipole(boxDipole);
    }

    return snap->getSystemDipole();
  }

  Mat3x3d Thermo::getSystemQuadrupole() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasSystemQuadrupole) {
      SimInfo::MoleculeIterator miter;
      vector<Atom*>::iterator aiter;
      Molecule* mol;
      Atom* atom;
      RealType charge;
      Vector3d ri(0.0);
      Vector3d dipole(0.0);
      Mat3x3d qpole(0.0);

      RealType chargeToC   = 1.60217733e-19;
      RealType angstromToM = 1.0e-10;
      RealType debyeToCm   = 3.33564095198e-30;

      for (mol = info_->beginMolecule(miter); mol != NULL;
           mol = info_->nextMolecule(miter)) {
        for (atom = mol->beginAtom(aiter); atom != NULL;
             atom = mol->nextAtom(aiter)) {
          ri = atom->getPos();
          snap->wrapVector(ri);
          ri *= angstromToM;

          charge = 0.0;

          FixedChargeAdapter fca = FixedChargeAdapter(atom->getAtomType());
          if (fca.isFixedCharge()) { charge = fca.getCharge(); }

          FluctuatingChargeAdapter fqa =
              FluctuatingChargeAdapter(atom->getAtomType());
          if (fqa.isFluctuatingCharge()) { charge += atom->getFlucQPos(); }

          charge *= chargeToC;

          qpole += 0.5 * charge * outProduct(ri, ri);

          MultipoleAdapter ma = MultipoleAdapter(atom->getAtomType());

          if (ma.isDipole()) {
            dipole = atom->getDipole() * debyeToCm;
            qpole += 0.5 * outProduct(dipole, ri);
            qpole += 0.5 * outProduct(ri, dipole);
          }

          if (ma.isQuadrupole()) {
            qpole += atom->getQuadrupole() * debyeToCm * angstromToM;
          }
        }
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, qpole.getArrayPointer(), 9, MPI_REALTYPE,
                    MPI_SUM, MPI_COMM_WORLD);
#endif

      snap->setSystemQuadrupole(qpole);
    }

    return snap->getSystemQuadrupole();
  }

  // Returns the Heat Flux Vector for the system
  Vector3d Thermo::getHeatFlux() {
    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    SimInfo::MoleculeIterator miter;
    vector<StuntDouble*>::iterator iiter;
    Molecule* mol;
    StuntDouble* sd;
    RigidBody::AtomIterator ai;
    Atom* atom;
    Vector3d vel;
    Vector3d angMom;
    Mat3x3d I;
    int i;
    int j;
    int k;
    RealType mass;

    Vector3d x_a;
    RealType kinetic;
    RealType potential;
    RealType eatom;
    // Convective portion of the heat flux
    Vector3d heatFluxJc = V3Zero;

    /* Calculate convective portion of the heat flux */
    for (mol = info_->beginMolecule(miter); mol != NULL;
         mol = info_->nextMolecule(miter)) {
      for (sd = mol->beginIntegrableObject(iiter); sd != NULL;
           sd = mol->nextIntegrableObject(iiter)) {
        mass = sd->getMass();
        vel  = sd->getVel();

        kinetic = mass * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);

        if (sd->isDirectional()) {
          angMom = sd->getJ();
          I      = sd->getI();

          if (sd->isLinear()) {
            i = sd->linearAxis();
            j = (i + 1) % 3;
            k = (i + 2) % 3;
            kinetic += angMom[j] * angMom[j] / I(j, j) +
                       angMom[k] * angMom[k] / I(k, k);
          } else {
            kinetic += angMom[0] * angMom[0] / I(0, 0) +
                       angMom[1] * angMom[1] / I(1, 1) +
                       angMom[2] * angMom[2] / I(2, 2);
          }
        }

        potential = 0.0;

        if (sd->isRigidBody()) {
          RigidBody* rb = dynamic_cast<RigidBody*>(sd);
          for (atom = rb->beginAtom(ai); atom != NULL;
               atom = rb->nextAtom(ai)) {
            potential += atom->getParticlePot();
          }
        } else {
          potential = sd->getParticlePot();
        }

        potential *= Constants::energyConvert;  // amu A^2/fs^2
        // The potential may not be a 1/2 factor
        eatom = (kinetic + potential) / 2.0;  // amu A^2/fs^2
        heatFluxJc[0] += eatom * vel[0];      // amu A^3/fs^3
        heatFluxJc[1] += eatom * vel[1];      // amu A^3/fs^3
        heatFluxJc[2] += eatom * vel[2];      // amu A^3/fs^3
      }
    }

    /* The J_v vector is reduced in the forceManager so everyone has
     *  the global Jv. Jc is computed over the local atoms and must be
     *  reduced among all processors.
     */
#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &heatFluxJc[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
#endif

    // (kcal/mol * A/fs) * conversion => (amu A^3)/fs^3

    Vector3d heatFluxJv =
        currSnapshot->getConductiveHeatFlux() * Constants::energyConvert;

    // Correct for the fact the flux is 1/V (Jc + Jv)
    return (heatFluxJv + heatFluxJc) / this->getVolume();  // amu / fs^3
  }

  Vector3d Thermo::getComVel() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasCOMvel) {
      SimInfo::MoleculeIterator i;
      Molecule* mol;

      Vector3d comVel(0.0);
      RealType totalMass(0.0);

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        RealType mass = mol->getMass();
        totalMass += mass;
        comVel += mass * mol->getComVel();
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, &totalMass, 1, MPI_REALTYPE, MPI_SUM,
                    MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, comVel.getArrayPointer(), 3, MPI_REALTYPE,
                    MPI_SUM, MPI_COMM_WORLD);
#endif

      comVel /= totalMass;
      snap->setCOMvel(comVel);
    }
    return snap->getCOMvel();
  }

  Vector3d Thermo::getCom() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasCOM) {
      SimInfo::MoleculeIterator i;
      Molecule* mol;

      Vector3d com(0.0);
      RealType totalMass(0.0);

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        RealType mass = mol->getMass();
        totalMass += mass;
        com += mass * mol->getCom();
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, &totalMass, 1, MPI_REALTYPE, MPI_SUM,
                    MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, com.getArrayPointer(), 3, MPI_REALTYPE,
                    MPI_SUM, MPI_COMM_WORLD);
#endif

      com /= totalMass;
      snap->setCOM(com);
    }
    return snap->getCOM();
  }

  /**
   * Returns center of mass and center of mass velocity in one
   * function call.
   */
  void Thermo::getComAll(Vector3d& com, Vector3d& comVel) {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!(snap->hasCOM && snap->hasCOMvel)) {
      SimInfo::MoleculeIterator i;
      Molecule* mol;

      RealType totalMass(0.0);

      com    = 0.0;
      comVel = 0.0;

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        RealType mass = mol->getMass();
        totalMass += mass;
        com += mass * mol->getCom();
        comVel += mass * mol->getComVel();
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, &totalMass, 1, MPI_REALTYPE, MPI_SUM,
                    MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, com.getArrayPointer(), 3, MPI_REALTYPE,
                    MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, comVel.getArrayPointer(), 3, MPI_REALTYPE,
                    MPI_SUM, MPI_COMM_WORLD);
#endif

      com /= totalMass;
      comVel /= totalMass;
      snap->setCOM(com);
      snap->setCOMvel(comVel);
    }
    com    = snap->getCOM();
    comVel = snap->getCOMvel();
    return;
  }

  /**
   * \brief Return inertia tensor for entire system and angular momentum
   *  Vector.
   *
   *
   *
   *    [  Ixx -Ixy  -Ixz ]
   * I =| -Iyx  Iyy  -Iyz |
   *    [ -Izx -Iyz   Izz ]
   */
  void Thermo::getInertiaTensor(Mat3x3d& inertiaTensor,
                                Vector3d& angularMomentum) {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    Molecule::IntegrableObjectIterator j;
    StuntDouble* sd;

    if (!(snap->hasInertiaTensor && snap->hasCOMw)) {
      RealType xx = 0.0;
      RealType yy = 0.0;
      RealType zz = 0.0;
      RealType xy = 0.0;
      RealType xz = 0.0;
      RealType yz = 0.0;
      Vector3d com(0.0);
      Vector3d comVel(0.0);

      getComAll(com, comVel);

      SimInfo::MoleculeIterator i;
      Molecule* mol;

      Vector3d r(0.0);
      Vector3d v(0.0);
      Mat3x3d Itmp {};

      RealType m = 0.0;

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (sd = mol->beginIntegrableObject(j); sd != NULL;
             sd = mol->nextIntegrableObject(j)) {
          r = sd->getPos() - com;
          v = sd->getVel() - comVel;

          m = sd->getMass();

          // Compute moment of intertia coefficients.
          xx += r[0] * r[0] * m;
          yy += r[1] * r[1] * m;
          zz += r[2] * r[2] * m;

          // compute products of intertia
          xy += r[0] * r[1] * m;
          xz += r[0] * r[2] * m;
          yz += r[1] * r[2] * m;

          angularMomentum += cross(r, v) * m;

          if (sd->isDirectional()) {
            RotMat3x3d A = sd->getA();
            Itmp += A.transpose() * sd->getI() * A;
            angularMomentum += sd->getJ();
          }
        }
      }

      inertiaTensor(0, 0) = yy + zz;
      inertiaTensor(0, 1) = -xy;
      inertiaTensor(0, 2) = -xz;
      inertiaTensor(1, 0) = -xy;
      inertiaTensor(1, 1) = xx + zz;
      inertiaTensor(1, 2) = -yz;
      inertiaTensor(2, 0) = -xz;
      inertiaTensor(2, 1) = -yz;
      inertiaTensor(2, 2) = xx + yy;

      inertiaTensor += Itmp;

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, inertiaTensor.getArrayPointer(), 9,
                    MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, angularMomentum.getArrayPointer(), 3,
                    MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
#endif

      snap->setCOMw(angularMomentum);
      snap->setInertiaTensor(inertiaTensor);
    }

    angularMomentum = snap->getCOMw();
    inertiaTensor   = snap->getInertiaTensor();

    return;
  }

  Mat3x3d Thermo::getBoundingBox() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!(snap->hasBoundingBox)) {
      SimInfo::MoleculeIterator i;
      Molecule::RigidBodyIterator ri;
      Molecule::AtomIterator ai;
      Molecule* mol;
      RigidBody* rb;
      Atom* atom;
      Vector3d pos, bMax, bMin;
      int index = 0;

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        // change the positions of atoms which belong to the rigidbodies
        for (rb = mol->beginRigidBody(ri); rb != NULL;
             rb = mol->nextRigidBody(ri)) {
          rb->updateAtoms();
        }

        for (atom = mol->beginAtom(ai); atom != NULL;
             atom = mol->nextAtom(ai)) {
          pos = atom->getPos();

          if (index == 0) {
            bMax = pos;
            bMin = pos;
          } else {
            for (int i = 0; i < 3; i++) {
              bMax[i] = max(bMax[i], pos[i]);
              bMin[i] = min(bMin[i], pos[i]);
            }
          }
          index++;
        }
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, &bMax[0], 3, MPI_REALTYPE, MPI_MAX,
                    MPI_COMM_WORLD);

      MPI_Allreduce(MPI_IN_PLACE, &bMin[0], 3, MPI_REALTYPE, MPI_MIN,
                    MPI_COMM_WORLD);
#endif
      Mat3x3d bBox = Mat3x3d(0.0);
      for (int i = 0; i < 3; i++) {
        bBox(i, i) = bMax[i] - bMin[i];
      }
      snap->setBoundingBox(bBox);
    }

    return snap->getBoundingBox();
  }

  // Returns the angular momentum of the system
  Vector3d Thermo::getAngularMomentum() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasCOMw) {
      Vector3d com(0.0);
      Vector3d comVel(0.0);
      Vector3d angularMomentum(0.0);

      getComAll(com, comVel);

      SimInfo::MoleculeIterator i;
      Molecule* mol;

      Vector3d thisr(0.0);
      Vector3d thisp(0.0);

      RealType thisMass;

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        thisMass = mol->getMass();
        thisr    = mol->getCom() - com;
        thisp    = (mol->getComVel() - comVel) * thisMass;

        angularMomentum += cross(thisr, thisp);
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, angularMomentum.getArrayPointer(), 3,
                    MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
#endif

      snap->setCOMw(angularMomentum);
    }

    return snap->getCOMw();
  }

  /**
   * Returns the Volume of the system based on a ellipsoid with
   * semi-axes based on the radius of gyration V=4/3*Pi*R_1*R_2*R_3
   * where R_i are related to the principle inertia moments
   *  R_i = sqrt(C*I_i/N), this reduces to
   *  V = 4/3*Pi*(C/N)^3/2*sqrt(det(I)).
   * See S.E. Baltazar et. al. Comp. Mat. Sci. 37 (2006) 526-536.
   */
  RealType Thermo::getGyrationalVolume() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasGyrationalVolume) {
      Mat3x3d intTensor;
      RealType det;
      Vector3d dummyAngMom;
      RealType sysconstants;
      RealType geomCnst;
      RealType volume;

      geomCnst = 3.0 / 2.0;
      /* Get the inertial tensor and angular momentum for free*/
      getInertiaTensor(intTensor, dummyAngMom);

      det = intTensor.determinant();
      sysconstants =
          geomCnst / (RealType)(info_->getNGlobalIntegrableObjects());
      volume =
          4.0 / 3.0 * Constants::PI * pow(sysconstants, geomCnst) * sqrt(det);

      snap->setGyrationalVolume(volume);
    }
    return snap->getGyrationalVolume();
  }

  void Thermo::getGyrationalVolume(RealType& volume, RealType& detI) {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!(snap->hasInertiaTensor && snap->hasGyrationalVolume)) {
      Mat3x3d intTensor;
      Vector3d dummyAngMom;
      RealType sysconstants;
      RealType geomCnst;

      geomCnst = 3.0 / 2.0;
      /* Get the inertia tensor and angular momentum for free*/
      this->getInertiaTensor(intTensor, dummyAngMom);

      detI = intTensor.determinant();
      sysconstants =
          geomCnst / (RealType)(info_->getNGlobalIntegrableObjects());
      volume =
          4.0 / 3.0 * Constants::PI * pow(sysconstants, geomCnst) * sqrt(detI);
      snap->setGyrationalVolume(volume);
    } else {
      volume = snap->getGyrationalVolume();
      detI   = snap->getInertiaTensor().determinant();
    }
    return;
  }

  RealType Thermo::getTaggedAtomPairDistance() {
    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    Globals* simParams     = info_->getSimParams();

    if (simParams->haveTaggedAtomPair() &&
        simParams->havePrintTaggedPairDistance()) {
      if (simParams->getPrintTaggedPairDistance()) {
        pair<int, int> tap = simParams->getTaggedAtomPair();
        Vector3d pos1, pos2, rab;

#ifdef IS_MPI
        int mol1 = info_->getGlobalMolMembership(tap.first);
        int mol2 = info_->getGlobalMolMembership(tap.second);

        int proc1 = info_->getMolToProc(mol1);
        int proc2 = info_->getMolToProc(mol2);

        RealType data[3];
        if (proc1 == worldRank) {
          StuntDouble* sd1 = info_->getIOIndexToIntegrableObject(tap.first);
          pos1             = sd1->getPos();
          data[0]          = pos1.x();
          data[1]          = pos1.y();
          data[2]          = pos1.z();
          MPI_Bcast(data, 3, MPI_REALTYPE, proc1, MPI_COMM_WORLD);
        } else {
          MPI_Bcast(data, 3, MPI_REALTYPE, proc1, MPI_COMM_WORLD);
          pos1 = Vector3d(data);
        }

        if (proc2 == worldRank) {
          StuntDouble* sd2 = info_->getIOIndexToIntegrableObject(tap.second);
          pos2             = sd2->getPos();
          data[0]          = pos2.x();
          data[1]          = pos2.y();
          data[2]          = pos2.z();
          MPI_Bcast(data, 3, MPI_REALTYPE, proc2, MPI_COMM_WORLD);
        } else {
          MPI_Bcast(data, 3, MPI_REALTYPE, proc2, MPI_COMM_WORLD);
          pos2 = Vector3d(data);
        }
#else
        StuntDouble* at1 = info_->getIOIndexToIntegrableObject(tap.first);
        StuntDouble* at2 = info_->getIOIndexToIntegrableObject(tap.second);
        pos1             = at1->getPos();
        pos2             = at2->getPos();
#endif
        rab = pos2 - pos1;
        currSnapshot->wrapVector(rab);
        return rab.length();
      }
      return 0.0;
    }
    return 0.0;
  }

  RealType Thermo::getHullVolume() {
#ifdef HAVE_QHULL
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    if (!snap->hasHullVolume) {
      Hull* surfaceMesh_;

      Globals* simParams   = info_->getSimParams();
      const std::string ht = simParams->getHULL_Method();

      if (ht == "Convex") {
        surfaceMesh_ = new ConvexHull();
      } else if (ht == "AlphaShape") {
        surfaceMesh_ = new AlphaHull(simParams->getAlpha());
      } else {
        return 0.0;
      }

      // Build a vector of stunt doubles to determine if they are
      // surface atoms
      std::vector<StuntDouble*> localSites_;
      Molecule* mol;
      StuntDouble* sd;
      SimInfo::MoleculeIterator i;
      Molecule::IntegrableObjectIterator j;

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (sd = mol->beginIntegrableObject(j); sd != NULL;
             sd = mol->nextIntegrableObject(j)) {
          localSites_.push_back(sd);
        }
      }

      // Compute surface Mesh
      surfaceMesh_->computeHull(localSites_);
      snap->setHullVolume(surfaceMesh_->getVolume());

      delete surfaceMesh_;
    }

    return snap->getHullVolume();
#else
    return 0.0;
#endif
  }

}  // namespace OpenMD
