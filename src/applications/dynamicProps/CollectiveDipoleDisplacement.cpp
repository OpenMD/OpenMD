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

#include "applications/dynamicProps/CollectiveDipoleDisplacement.hpp"

#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/Revision.hpp"

using namespace std;
namespace OpenMD {
  CollectiveDipoleDisplacement::CollectiveDipoleDisplacement(
      SimInfo* info, const string& filename, const string& sele1,
      const std::string& sele2) :
      SystemACF<Vector3d>(
          info, filename, sele1, sele2,
          DataStorage::dslPosition | DataStorage::dslFlucQPosition) {
    setCorrFuncType("Collective Dipole Displacement Function");
    setOutputName(getPrefix(dumpFilename_) + ".ddisp");

    std::stringstream label;
    label << "<|Mtrans(t)-Mtrans(0)|^2>\t"
          << "<|Mtot(t)-Mtot(0)|^2>\t"
          << "<|Mrot(t)-Mrot(0)|^2>";
    const std::string labelString = label.str();
    setLabelString(labelString);

    CRcm_.resize(nFrames_, V3Zero);
    CRtot_.resize(nFrames_, V3Zero);
    CRrot_.resize(nFrames_, V3Zero);

    // We'll need thermo to compute the volume:
    thermo_ = new Thermo(info_);
  }

  void CollectiveDipoleDisplacement::computeProperty1(int frame) {
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::AtomIterator ai;
    Atom* atom;
    AtomType* atype;

    RealType q, qtot;
    RealType m, mtot;
    Vector3d r(0.0), rcm(0.0), rcq(0.0);

    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      qtot = 0.0;
      mtot = 0.0;
      rcm *= 0.0;
      rcq *= 0.0;

      for (atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
        q                      = 0.0;
        atype                  = atom->getAtomType();
        FixedChargeAdapter fca = FixedChargeAdapter(atype);
        if (fca.isFixedCharge()) q = fca.getCharge();
        FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
        if (fqa.isFluctuatingCharge()) q += atom->getFlucQPos();

        r = atom->getPos();
        m = atom->getMass();

        qtot += q;
        mtot += m;

        rcm += r * m;
        rcq += r * q;
      }

      rcm /= mtot;

      if (qtot <= std::numeric_limits<RealType>::min()) {
        rcq = rcm;
      } else {
        rcq /= qtot;
      }

      CRcm_[frame] += qtot * rcm;
      CRtot_[frame] += qtot * rcq;
      CRrot_[frame] += qtot * (rcq - rcm);
      count_[frame]++;
    }

    RealType vol = thermo_->getVolume();

    CRcm_[frame] /= (vol * Constants::chargeDensityConvert);
    CRtot_[frame] /= (vol * Constants::chargeDensityConvert);
    CRrot_[frame] /= (vol * Constants::chargeDensityConvert);
  }

  Vector3d CollectiveDipoleDisplacement::calcCorrVal(int frame1, int frame2) {
    Vector3d diff;
    RealType dcm, dtot, drot;
    diff = CRcm_[frame2] - CRcm_[frame1];
    dcm  = diff.lengthSquare();
    diff = CRtot_[frame2] - CRtot_[frame1];
    dtot = diff.lengthSquare();
    diff = CRrot_[frame2] - CRrot_[frame1];
    drot = diff.lengthSquare();

    return Vector3d(dcm, dtot, drot);
  }
}  // namespace OpenMD
