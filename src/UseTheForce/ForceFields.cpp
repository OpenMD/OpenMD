#include <iostream>

using namespace std;


#include <stdlib.h>

#ifdef IS_MPI
#include <mpi.h>
#endif // is_mpi


#ifdef PROFILE
#include "profiling/mdProfile.hpp"
#endif

#include "utils/simError.h"
#include "UseTheForce/ForceFields.hpp"
#include "primitives/Atom.hpp"
#include "UseTheForce/doForces_interface.h"


void ForceFields::calcRcut( void ){

#ifdef IS_MPI
  double tempBig = bigSigma;
  MPI_Allreduce( &tempBig, &bigSigma, 1, MPI_DOUBLE, MPI_MAX,
		 MPI_COMM_WORLD);
#endif  //is_mpi

  //calc rCut and rList

  entry_plug->setDefaultRcut( 2.5 * bigSigma );  
    
}

void ForceFields::setRcut( double LJrcut ) {
  
#ifdef IS_MPI
  double tempBig = bigSigma;
  MPI_Allreduce( &tempBig, &bigSigma, 1, MPI_DOUBLE, MPI_MAX,
                 MPI_COMM_WORLD);
#endif  //is_mpi
  
  if (LJrcut < 2.5 * bigSigma) {
    sprintf( painCave.errMsg,
             "Setting Lennard-Jones cutoff radius to %lf.\n"
             "\tThis value is smaller than %lf, which is\n"
             "\t2.5 * bigSigma, where bigSigma is the largest\n"
             "\tvalue of sigma present in the simulation.\n"
             "\tThis is potentially a problem since the LJ potential may\n"
             "\tbe appreciable at this distance.  If you don't want the\n"
             "\tsmaller cutoff, change the LJrcut variable.\n",
             LJrcut, 2.5*bigSigma);
    painCave.isFatal = 0;
    simError();
  } else {
    sprintf( painCave.errMsg,
             "Setting Lennard-Jones cutoff radius to %lf.\n"
             "\tThis value is larger than %lf, which is\n"
             "\t2.5 * bigSigma, where bigSigma is the largest\n"
             "\tvalue of sigma present in the simulation. This should\n"
             "\tnot be a problem, but could adversely effect performance.\n",
             LJrcut, 2.5*bigSigma);
    painCave.isFatal = 0;
    simError();
  }
  
  //calc rCut and rList
  
  entry_plug->setDefaultRcut( LJrcut );
}

void ForceFields::doForces( int calcPot, int calcStress ){

  int i, j, isError;
  double* frc;
  double* pos;
  double* trq;
  double* A;
  double* u_l;
  double* rc;
  double* massRatio;
  double factor;
  SimState* config;

  Molecule* myMols;
  Atom** myAtoms;
  int numAtom;
  int curIndex;
  double mtot;
  int numMol;
  int numCutoffGroups;
  CutoffGroup* myCutoffGroup;
  vector<CutoffGroup*>::iterator iterCutoff;
  double com[3];
  vector<double> rcGroup;
  
  short int passedCalcPot = (short int)calcPot;
  short int passedCalcStress = (short int)calcStress;

  // forces are zeroed here, before any are accumulated.
  // NOTE: do not rezero the forces in Fortran.

  for(i=0; i<entry_plug->n_atoms; i++){
    entry_plug->atoms[i]->zeroForces();    
  }

#ifdef PROFILE
  startProfile(pro7);
#endif
  
  for(i=0; i<entry_plug->n_mol; i++ ){
    // CalcForces in molecules takes care of mapping rigid body coordinates
    // into atomic coordinates
    entry_plug->molecules[i].calcForces();    
  }

#ifdef PROFILE
  endProfile( pro7 );
#endif

  config = entry_plug->getConfiguration();
  
  frc = config->getFrcArray();
  pos = config->getPosArray();
  trq = config->getTrqArray();
  A   = config->getAmatArray();
  u_l = config->getUlArray();

  if(entry_plug->haveCutoffGroups){
    myMols = entry_plug->molecules;
    numMol = entry_plug->n_mol;
    for(int i  = 0; i < numMol; i++){
        
      numCutoffGroups = myMols[i].getNCutoffGroups();
      for(myCutoffGroup =myMols[i].beginCutoffGroup(iterCutoff); myCutoffGroup != NULL; 
                                                    myCutoffGroup =myMols[i].nextCutoffGroup(iterCutoff)){
        //get center of mass of the cutoff group
        myCutoffGroup->getCOM(com); 

        rcGroup.push_back(com[0]);
        rcGroup.push_back(com[1]);
        rcGroup.push_back(com[2]);
        
      }// end for(myCutoffGroup)
      
    }//end for(int i = 0)

    rc = &rcGroup[0];
  }
  else{
    // center of mass of the group is the same as position of the atom  if cutoff group does not exist
    rc = pos;
  }
  


  isError = 0;
  entry_plug->lrPot = 0.0;

  for (i=0; i<9; i++) {
    entry_plug->tau[i] = 0.0;
  }


#ifdef PROFILE
  startProfile(pro8);
#endif

  doForceLoop( pos,
  		    rc,
		    A,
		    u_l,
		    frc,
		    trq,
		    entry_plug->tau,
		    &(entry_plug->lrPot), 
		    &passedCalcPot,
		    &passedCalcStress,
		    &isError );


#ifdef PROFILE
  endProfile(pro8);
#endif


  if( isError ){
    sprintf( painCave.errMsg,
	     "Error returned from the fortran force calculation.\n" );
    painCave.isFatal = 1;
    simError();
  }

  // scale forces if thermodynamic integration is used
  if (entry_plug->useSolidThermInt || entry_plug->useLiquidThermInt) {
    factor = pow(entry_plug->thermIntLambda, entry_plug->thermIntK);
    for (i=0; i < entry_plug->n_atoms; i++) {
      for (j=0; j< 3; j++) 
	frc[3*i + j] *= factor; 
      if (entry_plug->atoms[i]->isDirectional()) {
	for (j=0; j< 3; j++) 
	  trq[3*i + j] *= factor;
      }
    }
    entry_plug->vRaw = entry_plug->lrPot;
    entry_plug->lrPot *= factor;
  }

  // collect the atomic forces onto rigid bodies
  for(i=0; i<entry_plug->n_mol; i++ ){
    entry_plug->molecules[i].atoms2rigidBodies();
  }

  // do crystal restraint forces for thermodynamic integration
  if (entry_plug->useSolidThermInt){
    entry_plug->lrPot += entry_plug->restraint->Calc_Restraint_Forces(entry_plug->integrableObjects);
    entry_plug->vHarm = entry_plug->restraint->getVharm();
  }
  

#ifdef IS_MPI
  sprintf( checkPointMsg,
	   "returned from the force calculation.\n" );
  MPIcheckPoint();
#endif // is_mpi
  

}


void ForceFields::initFortran(int ljMixPolicy, int useReactionField ){
  
  int isError;
  
  isError = 0;
  initFortranFF( &ljMixPolicy, &useReactionField, &isError );
  
  if(isError){
    sprintf( painCave.errMsg,
	     "ForceField error: There was an error initializing the forceField in fortran.\n" );
    painCave.isFatal = 1;
    simError();
  }

  
#ifdef IS_MPI
  sprintf( checkPointMsg, "ForceField successfully initialized the fortran component list.\n" );
  MPIcheckPoint();
#endif // is_mpi
  
}


void ForceFields::initRestraints(){
  int i;
  // store the initial info.
  // set the omega values to zero
  for (i=0; i<entry_plug->integrableObjects.size(); i++)
    entry_plug->integrableObjects[i]->setZangle( 0.0 );

  entry_plug->restraint->Store_Init_Info(entry_plug->integrableObjects);

}

void ForceFields::dumpzAngle(){

  // store the initial info.
  entry_plug->restraint->Write_zAngle_File(entry_plug->integrableObjects);

}
