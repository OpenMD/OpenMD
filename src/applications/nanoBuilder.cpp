#include <iostream>

#include <cstdlib>
#include <cstring>
#include <cmath>


#include "utils/simError.h"
#include "brains/SimInfo.hpp"
#include "io/ReadWrite.hpp"

#include "latticeBuilder.hpp"
#include "applications/MoLocator.hpp"
#include "sysBuild.hpp"
#include "applications/nanoBuilder.hpp"

nanoBuilder::nanoBuilder(int &hasError){
  int Errors;
  int foundCore,foundShell;
  int i;


  
  //Zero variables
  particleRadius = 0.0;
  coreRadius = 0.0;
  vacancyFraction = 0.0;
  vacancyRadius = 0.0;
  shellRadius = 0.0;
  latticeSpacing = 0.0;

  buildNmol = 0;

  nCoreMolecules = 0;
  nShellMolecules = 0;

  atomCount = 0;
  coreAtomCount = 0;
  shellAtomCount = 0;

  

  moleculeCount = 0; 
  foundCore  = 0;
  foundShell = 0;
  totalMolecules = 0;
  coreHasOrientation = 0;
  shellHasOrientation = 0;
  nInterface = 0;
  nMol = 0;

  hasError = 0;
  Errors = 0;
 
  //Initialize class members from bsInfo struct that sysbuilder provides.
  isRandom        = bsInfo.isRandomParticle;
  hasVacancies    = bsInfo.hasVacancies;
  latticeType     = bsInfo.latticeType;
  particleRadius  = bsInfo.particleRadius;
  coreRadius      = bsInfo.coreRadius;
  vacancyFraction = bsInfo.vacancyFraction;
  latticeSpacing  = bsInfo.latticeSpacing;
  soluteX         = bsInfo.soluteX; //Mole fraction for random particle.





  for (i=0;bsInfo.nComponents;i++){
    if( !strcmp( bsInfo.compStamps[i]->getID(),bsInfo.coreName )){
      foundCore = 1;
      coreStamp = bsInfo.compStamps[i];
      nCoreMolecules = bsInfo.componentsNmol[i];
    }
    if( !strcmp( bsInfo.compStamps[i]->getID(),bsInfo.shellName)){
      foundShell = 1;
      shellStamp = bsInfo.compStamps[i];
      nShellMolecules = bsInfo.componentsNmol[i];

    }    

  }



  if( !foundCore ){
    hasError = 1;
    return;
  }
  if( !foundShell ){
    hasError = 1;
    return;
  }



  Errors = sanityCheck();

  if (Errors){
    hasError = 1;
    return;
  }







  nCoreModelAtoms  = coreStamp->getNAtoms();
  nShellModelAtoms = shellStamp->getNAtoms();
  

  // We assume that if the core or shell model has more then one atom
  // the model has an orientational component...
  if (nCoreModelAtoms > 1)    coreHasOrientation = 1;
  if (nShellModelAtoms > 1)   shellHasOrientation = 1;
  
  maxModelNatoms = std::max(nCoreModelAtoms,nShellModelAtoms);
  
  /* If we specify a number of atoms in bass, we will try to build a nanopartice
     with that number.
  */
 

  if ((nShellMolecules != 0) && (nCoreMolecules != 0)){
    totalMolecules = nShellMolecules + nCoreMolecules;
    nCells = ceil(pow((double)totalMolecules/4.0, 1/3));
    buildNmol = 1;
  }
  else {
    nCells = 2.0 * particleRadius/latticeSpacing;
    shellRadius = particleRadius - coreRadius;
  }



 
  // Initialize random seed
  srand48( RAND_SEED );

  
}


nanoBuilder::~nanoBuilder(){
}


// Checks to make sure we aren't doing something the builder can't do.
int nanoBuilder::sanityCheck(void){

  // Right now we only do bimetallic nanoparticles  
  if (bsInfo.nComponents > 2) return 1;

  //Check for vacancies and random
  if (hasVacancies && isRandom) return 1;

  // make sure we aren't trying to build a core larger then the total particle size
  if ((coreRadius >= particleRadius) && (particleRadius != 0)) return 1;

  // we initialize the lattice spacing to be 0.0, if the lattice spacing is still 0.0
  // we have a problem
  if (latticeSpacing == 0.0) return 1;

  // Check to see if we are specifing the number of atoms in the particle correctly.
  if ((nShellMolecules == 0) && (nCoreMolecules != 0)){
    cerr << "nShellParticles is zero and nCoreParticles != 0" << "\n";
    return 1;
  } 
  // Make sure there are more then two components if we are building a randomly mixed particle.
  if ((bsInfo.nComponents < 2) && (isRandom)){
    cerr << "Two Components are needed to build a random particle." << "\n";
  }
  // Make sure both the core and shell models specify a target nmol.
  if ((nShellMolecules != 0) && (nCoreMolecules == 0)){
    cerr << "nCoreParticles is zero and nShellParticles != 0" << "\n";
    return 1;
  } 
  
  return 0;

}

    

int nanoBuilder::buildNanoParticle( void ){

  int ix;
  int iy;
  int iz;
  double *rx;
  double *ry;
  double *rz;
  double pos[3];
  double A[3][3];
  double HmatI[3][3];
  
  int nCellSites;
  int iref;
  int appNMols;
  int latticeCount = 0;
 
  int nAtoms;
  int nCoreAtomCounter = 0;
  int nShellAtomCounter = 0;
  int hasError;

  int i, j;

  int interfaceIndex = 0;
  double dist;
  double distsq;
  int latticeNpoints;
  int shesActualSizetoMe = 0;

  DumpWriter* writer;
  SimInfo* simnfo;
  SimState* theConfig;

  Lattice *myLattice;
  MoLocator *coreLocate;
  MoLocator *shellLocate;


  Atom** atoms;

  hasError = 0;

  myLattice = new Lattice(FCC_LATTICE_TYPE,latticeSpacing);
  /*
  latticeNpoints = myLattice.getNpoints();

  // Initializd atom vector to approximate size. 
  switch (buildType){

  case BUILD_NMOL_PARTICLE:

    break;
  case BUILD_CORE_SHELL_VACANCY:
   // Make space in the vector for all atoms except the last full cells
    // We will have to add at most (latticeNpoints-1)^3 to vector 
    appNMols = latticeNPoints * pow((double)(nCells - 1),3);
    moleculeVector.pushBack();

  default:
    // Make space in the vector for all atoms except the last full cells
    // We will have to add at most (latticeNpoints-1)^3 to vector 
    appNMols = latticeNPoints * pow((double)(nCells - 1),3);

  }
  */
  



  // Create molocator and atom arrays.
  coreLocate  = new MoLocator(coreStamp);
  shellLocate = new MoLocator(shellStamp);
 
  




  for(iz=-nCells;iz < nCells;iz++){ 
    for(iy=-nCells;iy<nCells;iy++){      
      for(ix=-nCells;ix<nCells;ix++){
	nCellSites = myLattice->getLatticePoints(&rx,&ry,&rz,
						 ix,iy,iz);
	for (iref=1;iref<nCellSites;iref++){
	  latticeCount++;
	  
	  pos[0] = rx[iref];
	  pos[1] = ry[iref];
	  pos[2] = rz[iref];
	  
	  distsq = rx[iref]*rx[iref] + ry[iref]*ry[iref] +rz[iref]*rz[iref];
	  dist = sqrt(distsq);
	  
	  switch(buildType){

	  case BUILD_CORE_SHELL:
	    nanoBuilder::buildWithCoreShell(dist,pos);
	    break;
	  case BUILD_CORE_SHELL_VACANCY:
	    nanoBuilder::buildWithVacancies(dist,pos);
	    break;

	  case BUILD_RANDOM_PARTICLE:
	    nanoBuilder::buildRandomlyMixed(dist,pos);
	    break;
	  case BUILD_NMOL_PARTICLE:
	    nanoBuilder::buildNmolParticle(dist,pos);
	  }
	}
      }
    }
  }



  // Create vacancies
  if (hasVacancies) buildVacancies();

  // Find the size of the atom vector not including Null atoms
  for (i=0;i<moleculeVector.size();i++){
    if (! moleculeVector[i].isVacancy){
      shesActualSizetoMe++;
      nAtoms = moleculeVector[i].myStamp->getNAtoms();
    }
  }

// Make a random particle.
  if (isRandom){
    placeRandom(shesActualSizetoMe); 

  // Loop back thru and count natoms since they may have changed
    for (i=0;i<moleculeVector.size();i++){
      if (! moleculeVector[i].isVacancy){
	shesActualSizetoMe++;
	nAtoms = moleculeVector[i].myStamp->getNAtoms();
      }
    }
  }


  // set up the SimInfo object

  simnfo = new SimInfo();
  simnfo->n_atoms = nAtoms;

  theConfig = simnfo->getConfiguration();
  theConfig->createArrays( nAtoms );
  simnfo->atoms = new Atom*[nAtoms];
  atoms = simnfo->atoms;
 

  shesActualSizetoMe = 0;
  /*  Use the information from the molecule vector to place the atoms.
   */
  for (i= 0;i<moleculeVector.size();i++){
    if (! moleculeVector[i].isVacancy) {
      orientationMunger( A );
      if( moleculeVector[i].isCore){
	nCoreAtomCounter += nCoreModelAtoms;
	coreLocate->placeMol(moleculeVector[i].pos,A,atoms,nShellAtomCounter, theConfig);
      }
      else {
	nShellAtomCounter += nShellModelAtoms;
	shellLocate->placeMol(moleculeVector[i].pos,A,atoms,nCoreAtomCounter, theConfig);
      }
      shesActualSizetoMe++;
    }
  }


  //      shellLocate.placeMol(pos, A, moleculeVector,shellAtomCount);

  for (i=0;i<3;i++) 
    for (j=0; j<3; j++) 
      simnfo->Hmat[i][j] = 0.0;

  simnfo->Hmat[0][0] = 1.0;
  simnfo->Hmat[1][1] = 1.0;
  simnfo->Hmat[2][2] = 1.0;
  

  
  sprintf( simnfo->sampleName, "%s.dump", bsInfo.outPrefix );
  sprintf( simnfo->finalName, "%s.init", bsInfo.outPrefix );

  // set up the writer and write out
  
  writer = new DumpWriter( simnfo );
  writer->writeFinal(0.0);

    // clean up  

  delete[] myLattice;

  return hasError;
} 

// Begin Builder routines------------------------------->

/* Builds a standard core-shell nanoparticle.
*/
void nanoBuilder::buildWithCoreShell(double dist, double pos[3]){

  
  if ( dist <= particleRadius ){
    moleculeVector.push_back(myMol);
    
    if (dist <= coreRadius){
      coreAtomCount += nCoreModelAtoms;
      moleculeVector[moleculeCount].pos[0] = pos[0]; 
      moleculeVector[moleculeCount].pos[1] = pos[1]; 
      moleculeVector[moleculeCount].pos[2] = pos[2]; 
      moleculeVector[moleculeCount].myStamp = coreStamp;
      moleculeVector[moleculeCount].isCore = 1;
      moleculeVector[moleculeCount].isShell = 0;
      
    }
    // Place shell
    else{
      shellAtomCount += nShellModelAtoms;
      moleculeVector[moleculeCount].pos[0] = pos[0]; 
      moleculeVector[moleculeCount].pos[1] = pos[1]; 
      moleculeVector[moleculeCount].pos[2] = pos[2]; 
      moleculeVector[moleculeCount].myStamp = shellStamp;
      moleculeVector[moleculeCount].isCore = 0;
      moleculeVector[moleculeCount].isShell = 1;
      
    }
    moleculeCount++;
  }	 
  
}
/*
Builds a core-shell nanoparticle and tracks the number of molecules at the
interface between the core-shell. These are recorded in vacancyInterface which is just
an integer vector. 
*/
void nanoBuilder::buildWithVacancies(double dist, double pos[3]){
  if ( dist <= particleRadius ){

    moleculeVector.push_back(myMol);
    if (dist <= coreRadius){
	
      coreAtomCount += nCoreModelAtoms;
      moleculeVector[moleculeCount].pos[0] = pos[0]; 
      moleculeVector[moleculeCount].pos[1] = pos[1]; 
      moleculeVector[moleculeCount].pos[2] = pos[2]; 
      moleculeVector[moleculeCount].myStamp = coreStamp;
      moleculeVector[moleculeCount].isCore = 1;
      moleculeVector[moleculeCount].isShell = 0;

      if ((dist >= coreRadius - vacancyRadius/2.0) && 
	  (dist <= coreRadius + vacancyRadius/2.0)){
	
	vacancyInterface.push_back(moleculeCount);
	nInterface++;
      }
    } else {
      // Place shell
      shellAtomCount += nShellModelAtoms;
      moleculeVector[moleculeCount].pos[0] = pos[0]; 
      moleculeVector[moleculeCount].pos[1] = pos[1]; 
      moleculeVector[moleculeCount].pos[2] = pos[2]; 
      moleculeVector[moleculeCount].myStamp = shellStamp;
      moleculeVector[moleculeCount].isCore = 0;
      moleculeVector[moleculeCount].isShell = 1;

    }
    moleculeCount++;
  }



}

/* Builds a core-shell nanoparticle where the number of core and shell
 molecules is known.
*/
void nanoBuilder::buildNmolParticle(double dist, double pos[3]){
  static int nMolCounter = 0;
  static int nCoreMolCounter = 0;
  
  
  if (nMolCounter < totalMolecules){
    moleculeVector.push_back(myMol);
    if (nCoreMolCounter < nCoreMolecules){
      
      coreAtomCount += nCoreModelAtoms;
      moleculeVector[moleculeCount].pos[0] = pos[0]; 
      moleculeVector[moleculeCount].pos[1] = pos[1]; 
      moleculeVector[moleculeCount].pos[2] = pos[2]; 
      moleculeVector[moleculeCount].myStamp = coreStamp;
      moleculeVector[moleculeCount].isCore = 1;
      moleculeVector[moleculeCount].isShell = 0;
      
      
    } else {
      shellAtomCount += nShellModelAtoms;
      moleculeVector[moleculeCount].pos[0] = pos[0]; 
      moleculeVector[moleculeCount].pos[1] = pos[1]; 
      moleculeVector[moleculeCount].pos[2] = pos[2]; 
      moleculeVector[moleculeCount].myStamp = shellStamp;
      moleculeVector[moleculeCount].isCore = 0;
      moleculeVector[moleculeCount].isShell = 1;
      
	
    }

  }
}


/* Builds a randomly mixed nanoparticle. We build the particle to be 
 entirely the core model, then randomly switch identities after the particle is built.
*/ 
void nanoBuilder::buildRandomlyMixed(double dist, double pos[3]){


  if ( dist <= particleRadius ){
    moleculeCount++;


    moleculeVector[moleculeCount].pos[0] = pos[0]; 
    moleculeVector[moleculeCount].pos[1] = pos[1]; 
    moleculeVector[moleculeCount].pos[2] = pos[2]; 
    moleculeVector[moleculeCount].myStamp = coreStamp;
    moleculeVector[moleculeCount].isCore = 1;
    moleculeVector[moleculeCount].isShell = 0;
    
  }	



}


// -----------------------END Builder routines.



//------------------------Begin Helper routines.
void nanoBuilder::placeRandom(int totalMol){
  int nSolute;
  int nSolvent;
  int i;
  int notfound;
  double solute_x;
  double solvent_x;
  
  int tester;

  nSolute = floor(soluteX * (double)totalMolecules); //CHECK ME
  nSolvent = totalMolecules - nSolute;
  
  solute_x = (double)nSolute/(double)totalMolecules;
  solvent_x = 1.0 - solute_x;
  
  
  

  for(i=0;nSolute-1;i++){
    notfound = 1;
    
    while(notfound){
      
      tester = floor((double)totalMolecules * drand48()); //Pick a molecule
      
      if (moleculeVector[tester].isCore){ // Make sure we select a core atom to change
	  
	moleculeVector[tester].isCore  = 0;
	moleculeVector[tester].isShell = 1;
	moleculeVector[tester].myStamp = shellStamp;
	notfound = 0; //set notfound = false.
      }
	
    }
      
  }
}


void nanoBuilder::buildVacancies(void){
  int i;  
  int* VacancyList; //logical nInterface long.
  int notfound;
  int index = 0;
  int nVacancies;
  int tester;

  if (nInterface != 0){
    nVacancies = floor((double)nInterface * vacancyFraction);

    VacancyList = new int[nInterface];
    
    // make vacancy list all false
    for(i=0;i<nInterface-1;i++){
      VacancyList[i] = 0;
    }
    
    // Build a vacancy list....    
    for(i=0;nVacancies-1;i++){
      notfound = 1;
      while(notfound){
	
	tester = floor((double)nInterface * drand48());
      
	if(! VacancyList[tester]){
	  VacancyList[tester] = 1;
	  notfound = 0;
	}
	
      }
    }
  }
  // Loop through and kill the vacancies from atom vector.

  for (i=0;i<nInterface;i++){
    if (VacancyList[i]){
      moleculeVector[vacancyInterface[i]].isVacancy = 1;   
    } // End Vacancy List
  } // for nInterface
  

    delete[] VacancyList;
}




void nanoBuilder::orientationMunger(double rot[3][3]){

  double theta, phi, psi;
  double cosTheta;

  // select random phi, psi, and cosTheta

  phi = 2.0 * M_PI * drand48();
  psi = 2.0 * M_PI * drand48();
  cosTheta = (2.0 * drand48()) - 1.0; // sample cos -1 to 1

  theta = acos( cosTheta );

  rot[0][0] = (cos(phi) * cos(psi)) - (sin(phi) * cos(theta) * sin(psi));
  rot[0][1] = (sin(phi) * cos(psi)) + (cos(phi) * cos(theta) * sin(psi));
  rot[0][2] = sin(theta) * sin(psi);

  rot[1][0] = -(cos(phi) * sin(psi)) - (sin(phi) * cos(theta) * cos(psi));
  rot[1][1] = -(sin(phi) * sin(psi)) + (cos(phi) * cos(theta) * cos(psi));
  rot[1][2] = sin(theta) * cos(psi);

  rot[2][0] = sin(phi) * sin(theta);
  rot[2][1] = -cos(phi) * sin(theta);
  rot[2][2] = cos(theta);

}


  




