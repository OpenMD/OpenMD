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

#include "perturbations/MagneticField.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/MultipoleAdapter.hpp"
#include "primitives/Molecule.hpp"
#include "nonbonded/NonBondedInteraction.hpp"
#include "utils/Constants.hpp"

namespace OpenMD {

  MagneticField::MagneticField(SimInfo* info) : initialized(false),
                                                doMagneticField(false),
                                                info_(info) {
    simParams = info_->getSimParams();
  }

  void MagneticField::initialize() {

    std::vector<RealType> mf;

    if (simParams->haveMagneticField()) {
      doMagneticField = true;
      mf = simParams->getMagneticField();
    }
    if (mf.size() != 3) {
      sprintf(painCave.errMsg,
              "MagneticField: Incorrect number of parameters specified.\n"
              "\tthere should be 3 parameters, but %lu were specified.\n", 
	      mf.size());
      painCave.isFatal = 1;
      simError();
    }
    MF.x() = mf[0];
    MF.y() = mf[1];
    MF.z() = mf[2];

    initialized = true;
  }

  void MagneticField::applyPerturbation() {

    if (!initialized) initialize();

    SimInfo::MoleculeIterator i;
    Molecule::AtomIterator  j;
    Molecule* mol;
    Atom* atom;
    AtomType* atype;

    int l,m,n;
    RealType C;
    Vector3d v;
    Vector3d f;
    Vector3d r;
    Vector3d t;
    Vector3d D;
    Vector3d AngMomentum;
    Vector3d omega;
    Mat3x3d I;
    bool isCharge;

    if (doMagneticField) {

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {

	for (atom = mol->beginAtom(j); atom != NULL;
	     atom = mol->nextAtom(j)) {

          isCharge = false;
          C = 0.0;

          atype = atom->getAtomType();
          r = atom->getPos();
          v = atom->getVel();

	  FixedChargeAdapter fca = FixedChargeAdapter(atype);
	  if ( fca.isFixedCharge() ) {
	    isCharge = true;
	    C = fca.getCharge();
	  }

          C *= Constants::chargeFieldConvert;

	  if (isCharge) {
	    f = cross(v, MF) * C * Constants::magneticFieldConvert;
	    atom->addFrc(f);
	  }

          MultipoleAdapter ma = MultipoleAdapter(atype);
          if (ma.isDipole() ) {
            D = atom->getDipole() * Constants::dipoleFieldConvert;

            t = cross(D, cross(v, MF));
            atom->addTrq(t);

            AngMomentum = atom->getJ();
            I = atom->getI();
            if(atom->isLinear()) {
              l = atom->linearAxis();
              m = (l + 1) % 3;
              n = (l + 2) % 3;
              omega[l] = 0;
              omega[m] = AngMomentum[m] / I(m,m);
              omega[n] = AngMomentum[n] / I(n,n);
            } else {
              omega = I.inverse() * AngMomentum;
            }

            f=cross(cross(omega, D), MF);
            atom->addFrc(f);
          }
        }
      }
    }
  }
}
