/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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

#include "applications/dynamicProps/ChargeKineticCorrFunc.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "primitives/Molecule.hpp"


namespace OpenMD {
  ChargeKineticCorrFunc::ChargeKineticCorrFunc(SimInfo* info,
                                               const std::string& filename,
                                               const std::string& sele1,
                                               const std::string& sele2,
                                               const RealType cutoff)
    : ObjectCCF<RealType>(info, filename, sele1, sele2,
                          DataStorage::dslFlucQPosition |
                          DataStorage::dslVelocity), rcut_(cutoff) {

    setCorrFuncType("Charge - Kinetic Cross Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".QKcorr");

    charges_.resize(nFrames_);
    kinetic_.resize(nFrames_);

    sumCharge_ = 0;
    sumKinetic_ = 0;
    chargeCount_ = 0;
    kineticCount_ = 0;
  }

  void ChargeKineticCorrFunc::validateSelection(SelectionManager& seleMan) {
    StuntDouble* sd;
    int i;

    for (sd = seleMan.beginSelected(i); sd != NULL;
         sd = seleMan.nextSelected(i)) {

      Atom* atom = static_cast<Atom*>(sd);
      AtomType* atomType = atom->getAtomType();
      FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);

      if (!sd->isDirectional() && !fqa.isFluctuatingCharge()) {
        sprintf(painCave.errMsg,
                "ChargeKineticCorrFunc::validateSelection Error: selection "
                "%d (%s)\n"
                "\t is not a Directional object\n", sd->getGlobalIndex(),
                sd->getType().c_str() );
        painCave.isFatal = 1;
        simError();
      }
    }
  }


  int ChargeKineticCorrFunc::computeProperty1(int frame, StuntDouble* sd) {
    RealType q = 0.0;
    pos1_ = sd->getPos();
    Atom* atom = static_cast<Atom*>(sd);

    AtomType* atomType = atom->getAtomType();

    FixedChargeAdapter fca = FixedChargeAdapter(atomType);
    if ( fca.isFixedCharge() ) {
      q += fca.getCharge();
    }

    FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
    if ( fqa.isFluctuatingCharge() ) {
      q += atom->getFlucQPos();
    }

    propertyTemp = q;
    charges_[frame].push_back( propertyTemp );
    sumCharge_ += propertyTemp;
    chargeCount_++;
    return charges_[frame].size() - 1;
  }

  int ChargeKineticCorrFunc::computeProperty2(int frame, StuntDouble* sd) {
    RealType kinetic(0.0);
    pos2_ = sd->getPos();

    RealType mass = sd->getMass();
    Vector3d vel = sd->getVel();

    kinetic += mass * (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);

    Vector3d angMom;
    Mat3x3d I;
    int i, j, k;

    if (sd->isDirectional()) {
      angMom = sd->getJ();
      I = sd->getI();

      if (sd->isLinear()) {
        i = sd->linearAxis();
        j = (i + 1) % 3;
        k = (i + 2) % 3;
        kinetic += angMom[j] * angMom[j] / I(j, j)
          + angMom[k] * angMom[k] / I(k, k);
      } else {
        kinetic += angMom[0]*angMom[0]/I(0, 0)
          + angMom[1]*angMom[1]/I(1, 1)
          + angMom[2]*angMom[2]/I(2, 2);
      }
    }
    propertyTemp = 0.5 * kinetic;
    kinetic_[frame].push_back( propertyTemp );
    sumKinetic_ += propertyTemp;
    kineticCount_++;
    return kinetic_[frame].size() - 1;
  }

  RealType ChargeKineticCorrFunc::calcCorrVal(int frame1, int frame2,
                                              int id1, int id2) {
    Vector3d sep = pos1_ - pos2_;
    RealType distance = sep.length();
    if (distance <= rcut_) return 0;
    else return charges_[frame1][id1] * kinetic_[frame2][id2] ;
  }

  void ChargeKineticCorrFunc::postCorrelate() {
    //gets the average of the charges
    sumCharge_ /= RealType(chargeCount_);

    //gets the average of the kinetic energy
    sumKinetic_ /= RealType(kineticCount_);

    RealType correlationOfAverages_ = sumCharge_ * sumKinetic_;
    for (unsigned int i =0 ; i < nTimeBins_; ++i) {
      if (count_[i] > 0) {

        histogram_[i] /= RealType(count_[i]);

        // The  correlation of the averages is subtracted
        // from the correlation value:
        histogram_[i] -= correlationOfAverages_;
      } else {
        histogram_[i] = 0;
      }
    }
  }
}
