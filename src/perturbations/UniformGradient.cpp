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

#include "perturbations/UniformGradient.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "types/MultipoleAdapter.hpp"
#include "primitives/Molecule.hpp"
#include "nonbonded/NonBondedInteraction.hpp"
#include "utils/Constants.hpp"
#ifdef IS_MPI
#include "mpi.h"
#endif

namespace OpenMD {
  
  UniformGradient::UniformGradient(SimInfo* info) : initialized(false),
						    doUniformGradient(false), 
						    doParticlePot(false),
						    info_(info) {
    simParams = info_->getSimParams();
  }
  
  void UniformGradient::initialize() {

    bool haveA = false;
    bool haveB = false;
    bool haveG = false;
    
    if (simParams->haveUniformGradientDirection1()) {
      std::vector<RealType> d1 = simParams->getUniformGradientDirection1();
      if (d1.size() != 3) {
        sprintf(painCave.errMsg,
                "uniformGradientDirection1: Incorrect number of parameters\n"
                "\tspecified. There should be 3 parameters, but %lu were\n"
                "\tspecified.\n", d1.size());
        painCave.isFatal = 1;
        simError();      
      }
      a_.x() = d1[0];
      a_.y() = d1[1];
      a_.z() = d1[2];

      a_.normalize();
      haveA = true;      
    }
    
    if (simParams->haveUniformGradientDirection2()) {
      std::vector<RealType> d2 = simParams->getUniformGradientDirection2();
      if (d2.size() != 3) {
        sprintf(painCave.errMsg,
                "uniformGradientDirection2: Incorrect number of parameters\n"
                "\tspecified. There should be 3 parameters, but %lu were\n"
                "\tspecified.\n", d2.size());
        painCave.isFatal = 1;
        simError();      
      }
      b_.x() = d2[0];
      b_.y() = d2[1];
      b_.z() = d2[2];

      b_.normalize();
      haveB = true;      
    }

    if (simParams->haveUniformGradientStrength()) {
      g_ = simParams->getUniformGradientStrength();
      haveG = true;
    }

    if (haveA && haveB && haveG) {
      doUniformGradient = true;
      cpsi_ = dot(a_, b_);
      
      Grad_(0,0) = 2.0 * (a_.x()*b_.x() - cpsi_ / 3.0);
      Grad_(0,1) = a_.x()*b_.y() + a_.y()*b_.x();
      Grad_(0,2) = a_.x()*b_.z() + a_.z()*b_.x();
      Grad_(1,0) = Grad_(0,1);
      Grad_(1,1) = 2.0 * (a_.y()*b_.y() - cpsi_ / 3.0);
      Grad_(1,2) = a_.y()*b_.z() + a_.z()*b_.y();
      Grad_(2,0) = Grad_(0,2);
      Grad_(2,1) = Grad_(1,2);
      Grad_(2,2) = 2.0 * (a_.z()*b_.z() - cpsi_ / 3.0);

      Grad_ *= g_ / 2.0;

    } else {
      if (!haveA) {
        sprintf(painCave.errMsg,
                "UniformGradient: uniformGradientDirection1 not specified.\n");
        painCave.isFatal = 1;
        simError();
      }
      if (!haveB) {
        sprintf(painCave.errMsg,
                "UniformGradient: uniformGradientDirection2 not specified.\n");
        painCave.isFatal = 1;
        simError();
      }
      if (!haveG) {
        sprintf(painCave.errMsg,
                "UniformGradient: uniformGradientStrength not specified.\n");
        painCave.isFatal = 1;
        simError();
      }
    }    
    
    int storageLayout_ = info_->getSnapshotManager()->getStorageLayout();
    if (storageLayout_ & DataStorage::dslParticlePot) doParticlePot = true;
    initialized = true;
  }
  
  void UniformGradient::applyPerturbation() {

    if (!initialized) initialize();

    SimInfo::MoleculeIterator i;
    Molecule::AtomIterator  j;
    Molecule* mol;
    Atom* atom;
    AtomType* atype;
    potVec longRangePotential(0.0);

    RealType C;
    Vector3d D;
    Mat3x3d Q;

    RealType U;
    RealType fPot;
    Vector3d t;
    Vector3d f;

    Vector3d r;
    Vector3d EF;

    bool isCharge;

    if (doUniformGradient) {

      U = 0.0;
      fPot = 0.0;

      for (mol = info_->beginMolecule(i); mol != NULL; 
           mol = info_->nextMolecule(i)) {      

	for (atom = mol->beginAtom(j); atom != NULL;
	     atom = mol->nextAtom(j)) {

          isCharge = false;
	  C = 0.0;
          
          atype = atom->getAtomType();

          r = atom->getPos();
	  info_->getSnapshotManager()->getCurrentSnapshot()->wrapVector(r);

          EF = Grad_ * r;
          
          atom->addElectricField(EF * Constants::chargeFieldConvert);
          
	  FixedChargeAdapter fca = FixedChargeAdapter(atype);
	  if ( fca.isFixedCharge() ) {
	    isCharge = true;
	    C = fca.getCharge();
	  }
	  
          FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
          if ( fqa.isFluctuatingCharge() ) {
	    isCharge = true;
            C += atom->getFlucQPos();
            atom->addFlucQFrc( dot(r, EF) 
                               * Constants::chargeFieldConvert );
          }
	  
	  if (isCharge) {
	    f = EF * C * Constants::chargeFieldConvert;
	    atom->addFrc(f);

	    U = -dot(r, f);
	    if (doParticlePot) {      
	      atom->addParticlePot(U);
	    }
	    fPot += U;
	  }
	    
          MultipoleAdapter ma = MultipoleAdapter(atype);
	  if (ma.isDipole() ) {
            D = atom->getDipole() * Constants::dipoleFieldConvert;
            
            f = D * Grad_;
            atom->addFrc(f);

	    t = cross(D, EF);
	    atom->addTrq(t);

	    U = -dot(D, EF);
	    if (doParticlePot) {      
	      atom->addParticlePot(U);
	    }
	    fPot += U;
	  }

          if (ma.isQuadrupole() ) {
            Q = atom->getQuadrupole() * Constants::dipoleFieldConvert;
            
            t = 2.0 * mCross(Q, Grad_);
            atom->addTrq(t);

            U = -doubleDot(Q, Grad_);
	    if (doParticlePot) {      
	      atom->addParticlePot(U);
	    }
	    fPot += U;
          }
	}
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, &fPot, 1, MPI_REALTYPE, 
                    MPI_SUM, MPI_COMM_WORLD);
#endif

      Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
      longRangePotential = snap->getLongRangePotentials();
      longRangePotential[ELECTROSTATIC_FAMILY] += fPot;
      snap->setLongRangePotential(longRangePotential);
    }
  }
}
