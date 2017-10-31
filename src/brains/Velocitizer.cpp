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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include "brains/Velocitizer.hpp"
#include "brains/Thermo.hpp"
#include "math/SquareMatrix3.hpp"
#include "utils/Constants.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/StuntDouble.hpp"


#ifndef IS_MPI
#include "math/SeqRandNumGen.hpp"
#else
#include "math/ParallelRandNumGen.hpp"
#endif

namespace OpenMD {
  
  Velocitizer::Velocitizer(SimInfo* info) : info_(info), thermo(info) {
    

    Globals * simParams = info->getSimParams();
    
#ifndef IS_MPI
    if (simParams->haveSeed()) {
      int seedValue = simParams->getSeed();
      randNumGen_ = new SeqRandNumGen(seedValue);
    }else {
      randNumGen_ = new SeqRandNumGen();
    }    
#else
    if (simParams->haveSeed()) {
      int seedValue = simParams->getSeed();
      randNumGen_ = new ParallelRandNumGen(seedValue);
    }else {
      randNumGen_ = new ParallelRandNumGen();
    }    
#endif 
  }
  
  Velocitizer::~Velocitizer() {
    delete randNumGen_;
  }

  void Velocitizer::scale(RealType lambda) {
   
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule * mol;
    StuntDouble * sd;
    Vector3d vcom, aJ;
    Mat3x3d I;
    int l, m, n;

    std::cerr << "V_lambda = " << lambda << endl;
    
    for( mol = info_->beginMolecule(i); mol != NULL;
	 mol = info_->nextMolecule(i) ) {

      for( sd = mol->beginIntegrableObject(j); sd != NULL;
	   sd = mol->nextIntegrableObject(j) ) {

	vcom = sd->getVel();
	std::cerr << "vcom = " << vcom << endl;
	
	// scale each component of the velocity by lambda
	for (int i = 0; i < vcom.size(); i++) {
	  vcom[i] = lambda * vcom[i];
	}
	std::cerr << "vcom = " << vcom << endl;

	sd->setVel(vcom);

	std::cerr << "setVel = " << sd->getVel() << endl << endl;
	std::cerr << "V_KE = " << thermo.getKinetic() << endl;
	
	if (sd->isDirectional()) {
	  I = sd->getI();
	  
	  if (sd->isLinear()) {
	    l = sd->linearAxis();
	    m = (l + 1) % 3;
	    n = (l + 2) % 3;
	    
	    aJ[l] = aJ[l];
	    aJ[m] = lambda * aJ[m];
	    aJ[n] = lambda * aJ[n];
	  }
	  else {
	    for( int k = 0; k < 3; k++ ) {
	      aJ[k] = lambda * aJ[k];
	    }
	  }
	  
	  sd->setJ(aJ);
	}
	
      }
    }
    
  }
      
  void Velocitizer::randomize(RealType temperature) {
    Vector3d aVel;
    Vector3d aJ;
    Mat3x3d I;
    int l, m, n;
    Vector3d vdrift;
    RealType vbar;
    RealType jbar;
    RealType av2;
    RealType kebar;
    
    Globals * simParams = info_->getSimParams();
    
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule * mol;
    StuntDouble * sd;

    kebar = Constants::kB * temperature * info_->getNdfRaw() /
      (2.0 * info_->getNdf());

    for( mol = info_->beginMolecule(i); mol != NULL;
	 mol = info_->nextMolecule(i) ) {

      for( sd = mol->beginIntegrableObject(j); sd != NULL;
	   sd = mol->nextIntegrableObject(j) ) {
	
	// uses equipartition theory to solve for vbar in angstrom/fs
	
	av2 = 2.0 * kebar / sd->getMass();
	vbar = sqrt(av2);
	
	// picks random velocities from a gaussian distribution
	// centered on vbar
	
	for( int k = 0; k < 3; k++ ) {
	  aVel[k] = vbar * randNumGen_->randNorm(0.0, 1.0);
	}
	sd->setVel(aVel);
	
	if (sd->isDirectional()) {
	  I = sd->getI();
	  
	  if (sd->isLinear()) {
	    l = sd->linearAxis();
	    m = (l + 1) % 3;
	    n = (l + 2) % 3;
	    
	    aJ[l] = 0.0;
	    jbar = sqrt(2.0 * kebar * I(m, m));
	    aJ[m] = jbar * randNumGen_->randNorm(0.0, 1.0);
	    jbar = sqrt(2.0 * kebar * I(n, n));
	    aJ[n] = jbar * randNumGen_->randNorm(0.0, 1.0);
	  } else {
	    for( int k = 0; k < 3; k++ ) {
	      jbar = sqrt(2.0 * kebar * I(k, k));
	      aJ[k] = jbar *randNumGen_->randNorm(0.0, 1.0);
	    }
	  }
	  
	  sd->setJ(aJ);
	}
      }
    }
        
    removeComDrift();

    // Remove angular drift if we are not using periodic boundary
    // conditions:

    if(!simParams->getUsePeriodicBoundaryConditions()) removeAngularDrift();    
  }
    
  void Velocitizer::removeComDrift() {
    // Get the Center of Mass drift velocity.
    Vector3d vdrift = thermo.getComVel();
    
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule * mol;
    StuntDouble * sd;
    
    //  Corrects for the center of mass drift.
    // sums all the momentum and divides by total mass.
    for( mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i) ) {

      for( sd = mol->beginIntegrableObject(j); sd != NULL;
	   sd = mol->nextIntegrableObject(j) ) {

	sd->setVel(sd->getVel() - vdrift);

      }
    }    
  }
    
  void Velocitizer::removeAngularDrift() {
    // Get the Center of Mass drift velocity.
      
    Vector3d vdrift;
    Vector3d com; 
      
    thermo.getComAll(com, vdrift);
         
    Mat3x3d inertiaTensor;
    Vector3d angularMomentum;
    Vector3d omega;
               
    thermo.getInertiaTensor(inertiaTensor, angularMomentum);

    // We now need the inverse of the inertia tensor.
    inertiaTensor = inertiaTensor.inverse();
    omega = inertiaTensor * angularMomentum;
    
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule* mol;
    StuntDouble* sd;
    Vector3d tempComPos;
    
    // Corrects for the center of mass angular drift by summing all
    // the angular momentum and dividing by the total mass.

    for( mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i) ) {

      for( sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j) ) {

        tempComPos = sd->getPos() - com;
        sd->setVel((sd->getVel() - vdrift) - cross(omega, tempComPos));
        
      }
    }   
  }   
}
