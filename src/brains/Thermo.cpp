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

  /**
   * @brief Constructs a Thermo helper bound to a SimInfo instance.
   *
   * Constructs the (optional) background VelocityField and caches whether
   * peculiar velocities should be used.  When a velocity field is active,
   * kinetic energies, temperatures, the kinetic part of the pressure tensor,
   * and the convective heat flux are all evaluated in the frame co-moving
   * with the local background flow.
   */
  Thermo::Thermo(SimInfo* info) : info_(info) {
    velField_              = std::make_unique<VelocityField>(info);
    usePeculiarVelocities_ = velField_->isActive();
  }

  /**
   * @brief Returns the peculiar (flow-subtracted) velocity at a position.
   *
   * When no background velocity field is active, @p vel is returned
   * unchanged.  Otherwise the local background flow velocity, sampled at
   * @p pos, is subtracted.
   *
   * @param vel lab-frame velocity (Angstrom/fs)
   * @param pos position at which to sample the background flow (Angstrom)
   * @return peculiar velocity in the co-moving frame (Angstrom/fs)
   */
  Vector3d Thermo::peculiarVelocity(const Vector3d& vel, const Vector3d& pos) {
    if (usePeculiarVelocities_) { return vel - velField_->getVelocity(pos); }
    return vel;
  }

  /**
   * @brief Returns the (peculiar) body-frame angular velocity of a
   * directional StuntDouble.
   *
   * The body-frame angular velocity is first obtained from the body-frame
   * angular momentum, \f$\omega_b = I^{-1} J\f$.  The local angular velocity
   * of the background flow, @p flowOmega, is then removed in the lab frame
   * and the result is transformed back to the body frame.  Using the
   * orthogonality of the rotation matrix \f$A\f$ this round trip reduces to
   * \f[
   *   \omega_b^{\text{pec}} = A\,(A^{T}\omega_b - \Omega_{\text{flow}})
   *                         = \omega_b - A\,\Omega_{\text{flow}}.
   * \f]
   *
   * @note @p flowOmega is the *angular velocity* of the background
   * flow, i.e. one half of the vorticity \f$\tfrac{1}{2}\nabla\times
   * v\f$ for an incompressible field. Functions calling this routine
   * should not use vorticity without modification.
   *
   * @param sd directional StuntDouble
   * @param flowOmega lab-frame angular velocity of the background flow
   * @return peculiar body-frame angular velocity
   */
  Vector3d Thermo::bodyAngularVelocity(StuntDouble* sd,
                                       const Vector3d& flowOmega) {
    Vector3d J = sd->getJ();
    Mat3x3d I  = sd->getI();
    Vector3d omegaB;

    if (sd->isLinear()) {
      int i     = sd->linearAxis();
      int j     = (i + 1) % 3;
      int k     = (i + 2) % 3;
      omegaB[i] = 0.0;
      omegaB[j] = J[j] / I(j, j);
      omegaB[k] = J[k] / I(k, k);
    } else {
      omegaB[0] = J[0] / I(0, 0);
      omegaB[1] = J[1] / I(1, 1);
      omegaB[2] = J[2] / I(2, 2);
    }

    return omegaB - sd->getA() * flowOmega;
  }

  /**
   * @brief Returns the translational kinetic energy of the system (kcal/mol).
   *
   * Uses peculiar velocities when a background velocity field is active so
   * that the streaming contribution is excluded.
   */
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
          vel  = peculiarVelocity(sd->getVel(), sd->getPos());

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

  /**
   * @brief Returns the rotational kinetic energy of the system (kcal/mol).
   *
   * When a background velocity field is active, the local flow rotation is
   * removed from each directional object's angular velocity (see
   * bodyAngularVelocity) before the kinetic energy is accumulated.
   */
  RealType Thermo::getRotationalKinetic() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasRotationalKineticEnergy) {
      SimInfo::MoleculeIterator miter;
      vector<StuntDouble*>::iterator iiter;
      Molecule* mol;
      StuntDouble* sd;
      Vector3d flowOmega(0.0);
      RealType kinetic(0.0);

      if (usePeculiarVelocities_) { flowOmega = 0.5 * velField_->getVorticity(); }

      for (mol = info_->beginMolecule(miter); mol != NULL;
           mol = info_->nextMolecule(miter)) {
        for (sd = mol->beginIntegrableObject(iiter); sd != NULL;
             sd = mol->nextIntegrableObject(iiter)) {
          if (sd->isDirectional()) {
            Vector3d omegaB = bodyAngularVelocity(sd, flowOmega);
            kinetic += dot(omegaB, sd->getI() * omegaB);
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

  /**
   * @brief Returns the total kinetic energy (translational + rotational +
   * electronic) of the system in kcal/mol.
   */
  RealType Thermo::getKinetic() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasKineticEnergy) {
      RealType ke = getTranslationalKinetic() + getRotationalKinetic() +
                    getElectronicKinetic();

      snap->setKineticEnergy(ke);
    }
    return snap->getKineticEnergy();
  }

  /**
   * @brief Returns the total potential energy (kcal/mol).
   *
   * ForceManager computes the potential and stores it in the Snapshot; this
   * accessor simply reports it.
   */
  RealType Thermo::getPotential() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    return snap->getPotentialEnergy();
  }

  /**
   * @brief Returns the per-selection potential energies.
   *
   * ForceManager computes the selection potentials and stores them in the
   * Snapshot; this accessor simply reports them.
   */
  potVec Thermo::getSelectionPotentials() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    return snap->getSelectionPotentials();
  }

  /**
   * @brief Returns the total energy (kinetic + potential) in kcal/mol.
   */
  RealType Thermo::getTotalEnergy() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasTotalEnergy) {
      snap->setTotalEnergy(this->getKinetic() + this->getPotential());
    }

    return snap->getTotalEnergy();
  }

  /**
   * @brief Returns the nuclear (translational + rotational) temperature in K.
   *
   * See getElectronicTemperature for the electronic portion.  When a
   * background velocity field is active this is the peculiar temperature,
   * since the underlying kinetic energies use peculiar velocities.  With no
   * degrees of freedom the temperature is not well defined and is set to
   * zero.
   */
  RealType Thermo::getTemperature() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasTemperature) {
      RealType nuclearKE =
          this->getTranslationalKinetic() + this->getRotationalKinetic();

      RealType temperature {};

      if (info_->getNdf() > 0)
        temperature = (2.0 * nuclearKE) / (info_->getNdf() * Constants::kb);
      else
        temperature = 0.0;

      snap->setTemperature(temperature);
    }

    return snap->getTemperature();
  }

  /**
   * @brief Returns the electronic (fluctuating-charge) kinetic energy in
   * kcal/mol.
   */
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

  /**
   * @brief Returns the electronic (fluctuating-charge) temperature in K.
   */
  RealType Thermo::getElectronicTemperature() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasElectronicTemperature) {
      RealType eTemp = (2.0 * this->getElectronicKinetic()) /
                       (info_->getNFluctuatingCharges() * Constants::kb);

      snap->setElectronicTemperature(eTemp);
    }

    return snap->getElectronicTemperature();
  }

  /**
   * @brief Returns the total net charge on the system (electrons).
   */
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

  /**
   * @brief Returns the instantaneous charge momentum in kcal fs / e / mol.
   */
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

  /**
   * @brief Returns the current density of the system.
   *
   * The first entry is the total current density; the remaining entries are
   * per-atom-type contributions, in the order given by
   * SimInfo::getSimulatedAtomTypes().  When a background velocity field is
   * active, peculiar velocities are used, so the returned quantity is the
   * conduction current in the frame co-moving with the local flow rather
   * than the lab-frame current.
   */
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
        Vector3d v = peculiarVelocity(atom->getVel(), atom->getPos());

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

  /**
   * @brief Returns the system volume in Ang^3.
   */
  RealType Thermo::getVolume() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    return snap->getVolume();
  }

  /**
   * @brief Returns the instantaneous pressure in atm (current snapshot).
   */
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

  /**
   * @brief Returns the instantaneous pressure in atm for a given snapshot.
   */
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

  /**
   * @brief Returns the pressure tensor in amu*fs^-2*Ang^-1 (current snapshot).
   *
   * Derived via the virial theorem description in Paci & Marchi,
   * J. Phys. Chem. 1996, 100, 4314-4322.  When a background velocity field is
   * active, the kinetic (streaming) part uses peculiar velocities so that the
   * bulk flow does not contribute a spurious ram pressure.
   */
  Mat3x3d Thermo::getPressureTensor() {
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
          vcom = peculiarVelocity(sd->getVel(), sd->getPos());
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

  /**
   * @brief Returns the pressure tensor in amu*fs^-2*Ang^-1 for a given
   * snapshot.
   *
   * Derived via the virial theorem description in Paci & Marchi,
   * J. Phys. Chem. 1996, 100, 4314-4322.  As with the current-snapshot
   * overload, the kinetic part uses peculiar velocities when a background
   * velocity field is active.
   */
  Mat3x3d Thermo::getPressureTensor(Snapshot* snap) {
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
          vcom = peculiarVelocity(sd->getVel(snap), sd->getPos(snap));
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

  /**
   * @brief Accumulates and returns the simulation box dipole moment in C*m.
   */
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

  /**
   * @brief Accumulates and returns the simulation box quadrupole moment.
   */
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

  /**
   * @brief Returns the heat flux vector for the system in amu/fs^3.
   *
   * The convective part (Jc) is accumulated over local integrable objects
   * and reduced across processors; the conductive part (Jv) is already
   * globally reduced in the ForceManager.  When a background velocity field
   * is active, peculiar translational and angular velocities are used so that
   * the flux is measured in the frame co-moving with the local flow.
   */
  Vector3d Thermo::getHeatFlux() {
    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    SimInfo::MoleculeIterator miter;
    vector<StuntDouble*>::iterator iiter;
    Molecule* mol;
    StuntDouble* sd;
    RigidBody::AtomIterator ai;
    Atom* atom;
    Vector3d vel;
    Vector3d flowOmega(0.0);
    RealType mass;
    RealType kinetic;
    RealType potential;
    RealType eatom;
    // Convective portion of the heat flux
    Vector3d heatFluxJc = V3Zero;

    if (usePeculiarVelocities_) { flowOmega = 0.5 * velField_->getVorticity(); }

    /* Calculate convective portion of the heat flux */
    for (mol = info_->beginMolecule(miter); mol != NULL;
         mol = info_->nextMolecule(miter)) {
      for (sd = mol->beginIntegrableObject(iiter); sd != NULL;
           sd = mol->nextIntegrableObject(iiter)) {
        mass = sd->getMass();
        vel  = peculiarVelocity(sd->getVel(), sd->getPos());

        kinetic = mass * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);

        if (sd->isDirectional()) {
          Vector3d omegaB = bodyAngularVelocity(sd, flowOmega);
          kinetic += dot(omegaB, sd->getI() * omegaB);
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

  /**
   * @brief Returns the center-of-mass velocity of the whole system.
   *
   * @note When a background velocity field is active, the returned (and
   * cached) value is the mass-weighted *peculiar* COM velocity.  Callers that
   * read Snapshot::getCOMvel() and expect the lab-frame drift velocity should
   * be aware of this.
   */
  Vector3d Thermo::getComVel() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasCOMvel) {
      SimInfo::MoleculeIterator i;
      Molecule* mol;

      Vector3d vel;
      Vector3d comVel(0.0);
      RealType totalMass(0.0);

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        RealType mass = mol->getMass();
        totalMass += mass;
        vel = peculiarVelocity(mol->getComVel(), mol->getCom());

        comVel += mass * vel;
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

  /**
   * @brief Returns the center of mass of the whole system.
   */
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
   * @brief Returns center of mass and center-of-mass velocity in one call.
   *
   * @param[out] com the center of mass of the whole system
   * @param[out] comVel the center-of-mass velocity of the whole system
   *
   * @note As with getComVel, the returned/cached @p comVel is the peculiar
   * COM velocity when a background velocity field is active.
   */
  void Thermo::getComAll(Vector3d& com, Vector3d& comVel) {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!(snap->hasCOM && snap->hasCOMvel)) {
      SimInfo::MoleculeIterator i;
      Molecule* mol;

      RealType totalMass(0.0);
      Vector3d vel;

      com    = 0.0;
      comVel = 0.0;

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        RealType mass = mol->getMass();
        totalMass += mass;
        com += mass * mol->getCom();

        vel = peculiarVelocity(mol->getComVel(), mol->getCom());

        comVel += mass * vel;
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
   * @brief Returns the inertia tensor and angular momentum for the entire
   * system.
   *
   * @param[out] inertiaTensor the inertia tensor
   *   \f[
   *     I = \begin{bmatrix} I_{xx} & -I_{xy} & -I_{xz} \\
   *                        -I_{yx} &  I_{yy} & -I_{yz} \\
   *                        -I_{zx} & -I_{yz} &  I_{zz} \end{bmatrix}
   *   \f]
   * @param[out] angularMomentum the angular momentum vector
   *
   * @note The COM velocity returned by getComAll is peculiar when a
   * background velocity field is active, so per-object velocities are also
   * taken peculiar here to keep the angular momentum in a single frame.
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
          v = peculiarVelocity(sd->getVel(), sd->getPos()) - comVel;

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

  /**
   * @brief Returns the axis-aligned bounding box for the current system.
   */
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

  /**
   * @brief Returns the angular momentum of the system.
   *
   * @note Uses peculiar molecular COM velocities when a background velocity
   * field is active, consistent with the peculiar COM velocity returned by
   * getComAll.
   */
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
        thisp    = (peculiarVelocity(mol->getComVel(), mol->getCom()) - comVel) *
                thisMass;

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
   * @brief Returns the volume of the system estimated from an ellipsoid with
   * semi-axes based on the radius of gyration.
   *
   * \f$V = \tfrac{4}{3}\pi R_1 R_2 R_3\f$ with \f$R_i = \sqrt{C\,I_i/N}\f$,
   * which reduces to \f$V = \tfrac{4}{3}\pi (C/N)^{3/2}\sqrt{\det I}\f$.
   * See S.E. Baltazar et al., Comp. Mat. Sci. 37 (2006) 526-536.
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

  /**
   * @brief Overloaded gyrational volume that also returns det(I) so that
   * dV/dr can be calculated.
   *
   * @param[out] volume the gyrational volume
   * @param[out] detI the determinant of the inertia tensor
   */
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

  /**
   * @brief Returns the distance between a tagged pair of atoms in Angstroms.
   *
   * Returns 0 when no tagged pair is configured or when printing of the
   * tagged-pair distance is disabled.
   */
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

  /**
   * @brief Returns the volume of the convex or alpha-shape hull enclosing the
   * system (requires QHULL).  Returns 0 when QHULL is unavailable or no valid
   * hull method is configured.
   */
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
