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
 */

#include "integrators/LHDForceModifier.hpp"

#include <cmath>

#include "primitives/Molecule.hpp"
#include "types/GayBerneAdapter.hpp"
#include "types/LennardJonesAdapter.hpp"
#include "utils/Constants.hpp"
#include "utils/ElementsTable.hpp"

namespace OpenMD {

  LHDForceModifier::LHDForceModifier(SimInfo* info) :
      ForceModifier {info}, maxIterNum_ {6}, forceTolerance_ {1e-6},
      simParams_ {info->getSimParams()},
      randNumGen_ {info->getRandomNumberGenerator()} {
    dt_  = simParams_->getDt();
    dt2_ = 0.5 * dt_;

    if (!simParams_->haveTargetTemp()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "LHDForceModifier: a targetTemp is required.\n");
      painCave.isFatal = 1;
      simError();
    }
    if (!simParams_->haveViscosity()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "LHDForceModifier: a viscosity is required.\n");
      painCave.isFatal = 1;
      simError();
    }
    // Convert the input viscosity (Poise) to OpenMD internal units, exactly as
    // Sphere::getHydroProp does, so the resistance R = M^{-1} is consistent
    // with LDForceModifier's Xitt. Folding the factor in here means the drag
    // (linear in R) and the random force (chol(R), so sqrt of the factor) both
    // scale correctly and the fluctuation-dissipation balance is preserved.
    viscosity_ = Constants::viscoConvert * simParams_->getViscosity();
    kT_        = Constants::kb * simParams_->getTargetTemp();

    velField_ = std::make_unique<VelocityField>(info);
    veloMunge_ = std::make_unique<Velocitizer>(info_);

    // Build one coupled group per molecule. Bead radii are fixed, so the
    // RPYMobility is constructed once and only its positions are refreshed.
    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ii;
    for (Molecule* mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      HydroMolecule hm;
      std::vector<RealType> radii;
      for (StuntDouble* sd = mol->beginIntegrableObject(ii); sd != NULL;
           sd = mol->nextIntegrableObject(ii)) {
        if (!sd->isAtom()) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "LHDForceModifier: only spherical atoms are supported as\n"
                   "\thydrodynamic beads (found a non-atom integrable object).\n");
          painCave.isFatal = 1;
          simError();
        }
        hm.beads.push_back(sd);
        hm.masses.push_back(sd->getMass());
        radii.push_back(beadRadius(static_cast<Atom*>(sd)));
      }
      if (!hm.beads.empty()) {
        hm.mobility = std::make_unique<RPYMobility>(radii, viscosity_);
        molecules_.push_back(std::move(hm));
      }
    }
  }

  RealType LHDForceModifier::beadRadius(Atom* atom) const {
    AtomType* atomType = atom->getAtomType();

    GayBerneAdapter gba = GayBerneAdapter(atomType);
    if (gba.isGayBerne()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "LHDForceModifier: Gay-Berne (non-spherical) atoms are not\n"
               "\tsupported by the translation-only bead model.\n");
      painCave.isFatal = 1;
      simError();
    }

    LennardJonesAdapter lja = LennardJonesAdapter(atomType);
    if (lja.isLennardJones()) return lja.getSigma() / 2.0;

    std::vector<AtomType*> atChain = atomType->allYourBase();
    for (std::vector<AtomType*>::iterator i = atChain.begin();
         i != atChain.end(); ++i) {
      int aNum = etab.GetAtomicNum((*i)->getName().c_str());
      if (aNum != 0) return etab.GetVdwRad(aNum);
    }

    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "LHDForceModifier: could not determine a hydrodynamic radius for\n"
             "\tatom type %s.\n",
             atomType->getName().c_str());
    painCave.isFatal = 1;
    simError();
    return 0.0;
  }

  void LHDForceModifier::modifyForces() {
    const RealType eConv = Constants::energyConvert;
    bool useFlow         = velField_->isActive();
    Mat3x3d E            = useFlow ? velField_->getRateOfStrain() : Mat3x3d(0.0);

    std::size_t molIndex = 0;
    for (HydroMolecule& hm : molecules_) {
      std::size_t N = hm.beads.size();
      RPYMobility& mob = *hm.mobility;

      // --- gather state ---
      std::vector<Vector3d> pos(N), vel(N), ambient(N, V3Zero);
      for (std::size_t i = 0; i < N; ++i) {
        pos[i] = hm.beads[i]->getPos();
        vel[i] = hm.beads[i]->getVel();
        if (useFlow) ambient[i] = velField_->getVelocity(pos[i]);
      }

      // --- rebuild mobility / resistance / Cholesky at this configuration ---
      // RPYC keeps M (and R = M^{-1}) positive definite by construction, so a
      // failure here is not an overlap artifact to tolerate: it means the
      // configuration handed in is bad (NaN / runaway coordinates) or the bead
      // radii / viscosity are unphysical. Fail loudly rather than let clamped,
      // FDT-violating noise enter the trajectory.
      if (!mob.update(pos)) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "LHDForceModifier: the resistance tensor for molecule %lu\n"
                 "\t(%lu beads) is not positive definite. The RPYC mobility is\n"
                 "\tSPD by construction, so this indicates invalid coordinates\n"
                 "\t(NaN or overflow) or unphysical bead radii / viscosity.\n",
                 static_cast<unsigned long>(molIndex),
                 static_cast<unsigned long>(N));
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }
      ++molIndex;

      // effective ambient velocity (affine flow + dipolar disturbance)
      std::vector<Vector3d> vEff =
          useFlow ? mob.effectiveAmbient(pos, ambient, E) : ambient;

      // --- correlated random force, applied once ---
      std::vector<RealType> Z(3 * N);
      for (std::size_t k = 0; k < 3 * N; ++k) Z[k] = normal_(*randNumGen_);
      std::vector<Vector3d> Frand = mob.randomForce(Z, kT_, dt_);
      for (std::size_t i = 0; i < N; ++i) hm.beads[i]->addFrc(Frand[i]);

      // --- self-consistent coupled friction solve ---
      // velocity is known at the half step; the friction needs the full-step
      // velocity, which depends on the friction. Iterate to convergence.
      std::vector<Vector3d> frc(N), velStep(N), Ffric(N, V3Zero), oldF(N);
      for (std::size_t i = 0; i < N; ++i) {
        frc[i]     = hm.beads[i]->getFrc();  // conservative + random
        velStep[i] = vel[i] + (dt2_ / hm.masses[i] * eConv) * frc[i];
      }

      for (int k = 0; k < maxIterNum_; ++k) {
        oldF = Ffric;
        // f_i = sum_j R_ij ( vEff_j - velStep_j )
        Ffric = mob.dragForce(vEff, velStep);
        for (std::size_t i = 0; i < N; ++i)
          velStep[i] =
              vel[i] + (dt2_ / hm.masses[i] * eConv) * (frc[i] + Ffric[i]);

        // converged when every bead's friction force has stopped changing
        // direction/magnitude (fdot -> 1)
        RealType worst = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
          RealType f2 = Ffric[i].lengthSquare();
          if (f2 < 1.0e-12) continue;  // negligible drag on this bead
          RealType fdot = dot(Ffric[i], oldF[i]) / f2;
          worst         = std::max(worst, std::fabs(1.0 - fdot));
        }
        if (worst <= forceTolerance_) break;
      }

      for (std::size_t i = 0; i < N; ++i) hm.beads[i]->addFrc(Ffric[i]);
    }

    // Drift removal only makes sense without an imposed flow; with a background
    // flow it would cancel the net advection the flow imparts.
    if (!useFlow) {
      if (simParams_->getConserveLinearMomentum()) veloMunge_->removeComDrift();
      if (!simParams_->getUsePeriodicBoundaryConditions() &&
          simParams_->getConserveAngularMomentum())
        veloMunge_->removeAngularDrift();
    }
  }
}  // namespace OpenMD
