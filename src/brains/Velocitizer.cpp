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

#include "brains/Velocitizer.hpp"
#include "brains/Thermo.hpp"
#include "math/SquareMatrix3.hpp"
#include "utils/Constants.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/StuntDouble.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "flucq/FluctuatingChargeConstraints.hpp"


#ifndef IS_MPI
#include "math/SeqRandNumGen.hpp"
#else
#include "math/ParallelRandNumGen.hpp"
#endif

namespace OpenMD {

  Velocitizer::Velocitizer(SimInfo* info) : info_(info), thermo_(info) {

    globals_ = info->getSimParams();

#ifndef IS_MPI
    if (globals_->haveSeed()) {
      int seedValue = globals_->getSeed();
      randNumGen_ = new SeqRandNumGen(seedValue);
    } else {
      randNumGen_ = new SeqRandNumGen();
    }
#else
    if (globals_->haveSeed()) {
      int seedValue = globals_->getSeed();
      randNumGen_ = new ParallelRandNumGen(seedValue);
    } else {
      randNumGen_ = new ParallelRandNumGen();
    }
#endif
  }

  Velocitizer::~Velocitizer() {
    delete randNumGen_;
  }

  void Velocitizer::scale(RealType lambda) {
    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ioi;
    Molecule * mol;
    StuntDouble * sd;
    Vector3d v, j;

    for( mol = info_->beginMolecule(mi); mol != NULL;
	 mol = info_->nextMolecule(mi) ) {

      for( sd = mol->beginIntegrableObject(ioi); sd != NULL;
	   sd = mol->nextIntegrableObject(ioi) ) {

	v = sd->getVel();
        v *= lambda;
	sd->setVel(v);

	if (sd->isDirectional()) {
          j = sd->getJ();
          j *= lambda;
	  sd->setJ(j);
	}
      }
    }

    removeComDrift();

    // Remove angular drift if we are not using periodic boundary
    // conditions:

    if(!globals_->getUsePeriodicBoundaryConditions()) removeAngularDrift();
  }

  void Velocitizer::randomize(RealType temperature) {
    Vector3d v;
    Vector3d j;
    Mat3x3d I;
    int l, m, n;
    Vector3d vdrift;
    RealType vbar;
    RealType jbar;
    RealType av2;
    RealType kebar;

    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ioi;
    Molecule * mol;
    StuntDouble * sd;

    kebar = Constants::kB * temperature * info_->getNdfRaw() /
      (2.0 * info_->getNdf());
    for( mol = info_->beginMolecule(mi); mol != NULL;
	 mol = info_->nextMolecule(mi) ) {

      for( sd = mol->beginIntegrableObject(ioi); sd != NULL;
	   sd = mol->nextIntegrableObject(ioi) ) {

	// uses equipartition theory to solve for vbar in angstrom/fs

	av2 = 2.0 * kebar / sd->getMass();
	vbar = sqrt(av2);

	// picks random velocities from a gaussian distribution
	// centered on vbar

	for( int k = 0; k < 3; k++ ) {
	  v[k] = vbar * randNumGen_->randNorm(0.0, 1.0);
	}
	sd->setVel(v);

	if (sd->isDirectional()) {
	  I = sd->getI();

	  if (sd->isLinear()) {
	    l = sd->linearAxis();
	    m = (l + 1) % 3;
	    n = (l + 2) % 3;

	    j[l] = 0.0;
	    jbar = sqrt(2.0 * kebar * I(m, m));
	    j[m] = jbar * randNumGen_->randNorm(0.0, 1.0);
	    jbar = sqrt(2.0 * kebar * I(n, n));
	    j[n] = jbar * randNumGen_->randNorm(0.0, 1.0);
	  } else {
	    for( int k = 0; k < 3; k++ ) {
	      jbar = sqrt(2.0 * kebar * I(k, k));
	      j[k] = jbar *randNumGen_->randNorm(0.0, 1.0);
	    }
	  }

	  sd->setJ(j);
	}
      }
    }

    removeComDrift();

    // Remove angular drift if we are not using periodic boundary
    // conditions:

    if(!globals_->getUsePeriodicBoundaryConditions()) removeAngularDrift();
  }


  void Velocitizer::randomizeChargeVelocity(RealType temperature) {
    RealType aw2;
    RealType kebar;
    RealType wbar;

    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ioi;
    Molecule * mol;
    StuntDouble * sd;
    FluctuatingChargeParameters* fqParams;
    FluctuatingChargeConstraints* fqConstraints;

    Globals* simParams = info_->getSimParams();
    fqParams = simParams->getFluctuatingChargeParameters();

    fqConstraints = new FluctuatingChargeConstraints(info_);
    fqConstraints->setConstrainRegions(fqParams->getConstrainRegions());

    int nConstrain =  fqConstraints->getNumberOfFlucQConstraints(); // no of constraints in charge
    int dfRaw =  fqConstraints->getNumberOfFlucQAtoms(); // no of FlucQ freedom
    int dfActual = dfRaw - nConstrain;
    kebar = dfRaw * Constants::kb * temperature / (2 * dfActual);

    for( mol = info_->beginMolecule(mi); mol != NULL;
   mol = info_->nextMolecule(mi) ) {

      for( sd = mol->beginIntegrableObject(ioi); sd != NULL;
     sd = mol->nextIntegrableObject(ioi) ) {


     if(sd->isAtom()){
       Atom * atom = static_cast<Atom*>(sd);
       AtomType* atomType = atom->getAtomType();
       FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
       if ( fqa.isFluctuatingCharge() ) {
         // uses equipartition theory to solve for vbar in angstrom/fs

         aw2 = 2.0 * kebar / atom->getChargeMass();
         wbar = sqrt(aw2);

         // picks random velocities from a gaussian distribution
         // centered on vbar
         atom-> setFlucQVel(wbar * randNumGen_->randNorm(0.0, 1.0));
       }
     }

// randomization of the charge velocities for atoms in the rigidbody

     if(sd->isRigidBody()){
       RigidBody* rigidbody = static_cast<RigidBody*>(sd);
       vector <Atom*> atomList;
       atomList = rigidbody->getAtoms();
       vector <Atom*>::iterator atom_iterator;
       for(size_t i = 0; i < atomList.size(); ++i){
         Atom *atom = atomList[i];
         AtomType* atomType = atom->getAtomType();
         FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
         if ( fqa.isFluctuatingCharge() ) {
           // uses equipartition theory to solve for vbar in angstrom/fs
           aw2 = 2.0 * kebar / atom->getChargeMass();
           wbar = sqrt(aw2);
           // picks random velocities from a gaussian distribution
           // centered on vbar
           atom-> setFlucQVel(wbar * randNumGen_->randNorm(0.0, 1.0));
         }
       }
     }
   }
 }
fqConstraints->applyConstraintsOnChargeVelocities();
}


  void Velocitizer::removeComDrift() {
    // Get the Center of Mass drift velocity.
    Vector3d vdrift = thermo_.getComVel();

    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ioi;
    Molecule * mol;
    StuntDouble * sd;

    //  Corrects for the center of mass drift.
    // sums all the momentum and divides by total mass.
    for( mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi) ) {

      for( sd = mol->beginIntegrableObject(ioi); sd != NULL;
	   sd = mol->nextIntegrableObject(ioi) ) {

	sd->setVel(sd->getVel() - vdrift);

      }
    }
  }

  void Velocitizer::removeAngularDrift() {
    // Get the Center of Mass drift velocity.

    Vector3d vdrift;
    Vector3d com;

    thermo_.getComAll(com, vdrift);

    Mat3x3d inertiaTensor;
    Vector3d angularMomentum;
    Vector3d omega;

    thermo_.getInertiaTensor(inertiaTensor, angularMomentum);

    // We now need the inverse of the inertia tensor.
    inertiaTensor = inertiaTensor.inverse();
    omega = inertiaTensor * angularMomentum;
    
    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ioi;
    Molecule* mol;
    StuntDouble* sd;
    Vector3d tempComPos;

    // Corrects for the center of mass angular drift by summing all
    // the angular momentum and dividing by the total mass.

    for( mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi) ) {

      for( sd = mol->beginIntegrableObject(ioi); sd != NULL;
           sd = mol->nextIntegrableObject(ioi) ) {

        tempComPos = sd->getPos() - com;
        sd->setVel((sd->getVel() - vdrift) - cross(omega, tempComPos));
      }
    }
  }
}
