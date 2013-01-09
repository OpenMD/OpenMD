/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
#include "perturbations/ElectricField.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "types/MultipoleAdapter.hpp"
#include "primitives/Molecule.hpp"
#include "nonbonded/NonBondedInteraction.hpp"

namespace OpenMD {

  ElectricField::ElectricField(SimInfo* info) : info_(info), 
						doElectricField(false), 
						doParticlePot(false),
						initialized(false) {
    simParams = info_->getSimParams();
  }

  void ElectricField::initialize() {
    if (simParams->haveElectricField()) {
      doElectricField = true;
      EF = simParams->getElectricField();
    }   
    int storageLayout_ = info_->getSnapshotManager()->getStorageLayout();
    if (storageLayout_ & DataStorage::dslParticlePot) doParticlePot = true;
    initialized = true;
  }

  void ElectricField::applyPerturbation() {
    if (!initialized) initialize();

    SimInfo::MoleculeIterator i;
    Molecule::AtomIterator  j;
    Molecule* mol;
    Atom* atom;
    potVec longRangePotential(0.0);
    Vector3d dip;
    Vector3d trq;
    Vector3d EFfrc;				
    Vector3d pos;
    RealType chrg;
    RealType pot, fieldPot;
    RealType chrgToKcal = 23.0609;
    RealType debyeToKcal = 4.8018969509;
    bool isCharge;

    if (doElectricField) {
      fieldPot = 0.0;

      for (mol = info_->beginMolecule(i); mol != NULL; mol = info_->nextMolecule(i)) {      
	for (atom = mol->beginAtom(j); atom != NULL;
	     atom = mol->nextAtom(j)) {
	  isCharge = false;
	  chrg = 0.0;

	  FixedChargeAdapter fca = FixedChargeAdapter(atom->getAtomType());
	  if ( fca.isFixedCharge() ) {
	    isCharge = true;
	    chrg = fca.getCharge();
	  }
	  
          FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atom->getAtomType());
          if ( fqa.isFluctuatingCharge() ) {
	    isCharge = true;
            chrg += atom->getFlucQPos();
          }
	  
	  if (isCharge) {
	    EFfrc = EF*chrg;
	    EFfrc *= chrgToKcal;
	    atom->addFrc(EFfrc);
	    // ad-hoc choice of the origin for potential calculation
	    pos = atom->getPos();
	    pot = -dot(pos, EFfrc);
	    if (doParticlePot) {      
	      atom->addParticlePot(pot);
	    }
	    fieldPot += pot;
	  }
	    
	  MultipoleAdapter ma = MultipoleAdapter(atom->getAtomType());
	  if (ma.isDipole() ) {
            Vector3d dipole = atom->getDipole();
	    dipole *= debyeToKcal;

	    trq = cross(dipole, EF);
	    atom->addTrq(trq);

	    pot = -dot(dipole, EF);
	    if (doParticlePot) {      
	      atom->addParticlePot(pot);
	    }
	    fieldPot += pot;
	  }
	}
      }
#ifdef IS_MPI
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &fieldPot, 1, MPI::REALTYPE, 
                                MPI::SUM);
#endif
      Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
      longRangePotential = snap->getLongRangePotentials();
      longRangePotential[ELECTROSTATIC_FAMILY] += fieldPot;
      snap->setLongRangePotential(longRangePotential);
    }
  }

}
