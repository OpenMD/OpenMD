#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream>
using namespace std;

#include "brains/SimInfo.hpp"
#define __C
#include "brains/fSimulation.h"
#include "utils/simError.h"
#include "UseTheForce/DarkSide/simulation_interface.h"
#include "UseTheForce/notifyCutoffs_interface.h"

//#include "UseTheForce/fortranWrappers.hpp"

#include "math/MatVec3.h"

#ifdef IS_MPI
#include "brains/mpiSimulation.hpp"
#endif

inline double roundMe( double x ){
  return ( x >= 0 ) ? floor( x + 0.5 ) : ceil( x - 0.5 );
}
	  
inline double min( double a, double b ){
  return (a < b ) ? a : b;
}

SimInfo* currentInfo;

SimInfo::SimInfo(){

  n_constraints = 0;
  nZconstraints = 0;
  n_oriented = 0;
  n_dipoles = 0;
  ndf = 0;
  ndfRaw = 0;
  nZconstraints = 0;
  the_integrator = NULL;
  setTemp = 0;
  thermalTime = 0.0;
  currentTime = 0.0;
  rCut = 0.0;
  rSw = 0.0;

  haveRcut = 0;
  haveRsw = 0;
  boxIsInit = 0;
  
  resetTime = 1e99;

  orthoRhombic = 0;
  orthoTolerance = 1E-6;
  useInitXSstate = true;

  usePBC = 0;
  useDirectionalAtoms = 0;
  useLennardJones = 0; 
  useElectrostatics = 0;
  useCharges = 0;
  useDipoles = 0;
  useSticky = 0;
  useGayBerne = 0;
  useEAM = 0;
  useShapes = 0;
  useFLARB = 0;

  useSolidThermInt = 0;
  useLiquidThermInt = 0;

  haveCutoffGroups = false;

  excludes = Exclude::Instance();

  myConfiguration = new SimState();

  has_minimizer = false;
  the_minimizer =NULL;

  ngroup = 0;

}


SimInfo::~SimInfo(){

  delete myConfiguration;

  map<string, GenericData*>::iterator i;
  
  for(i = properties.begin(); i != properties.end(); i++)
    delete (*i).second;

}

void SimInfo::setBox(double newBox[3]) {
  
  int i, j;
  double tempMat[3][3];

  for(i=0; i<3; i++) 
    for (j=0; j<3; j++) tempMat[i][j] = 0.0;;

  tempMat[0][0] = newBox[0];
  tempMat[1][1] = newBox[1];
  tempMat[2][2] = newBox[2];

  setBoxM( tempMat );

}

void SimInfo::setBoxM( double theBox[3][3] ){
  
  int i, j;
  double FortranHmat[9]; // to preserve compatibility with Fortran the
                         // ordering in the array is as follows:
                         // [ 0 3 6 ]
                         // [ 1 4 7 ]
                         // [ 2 5 8 ]
  double FortranHmatInv[9]; // the inverted Hmat (for Fortran);

  if( !boxIsInit ) boxIsInit = 1;

  for(i=0; i < 3; i++) 
    for (j=0; j < 3; j++) Hmat[i][j] = theBox[i][j];
  
  calcBoxL();
  calcHmatInv();

  for(i=0; i < 3; i++) {
    for (j=0; j < 3; j++) {
      FortranHmat[3*j + i] = Hmat[i][j];
      FortranHmatInv[3*j + i] = HmatInv[i][j];
    }
  }

  setFortranBox(FortranHmat, FortranHmatInv, &orthoRhombic);
 
}
 

void SimInfo::getBoxM (double theBox[3][3]) {

  int i, j;
  for(i=0; i<3; i++) 
    for (j=0; j<3; j++) theBox[i][j] = Hmat[i][j];
}


void SimInfo::scaleBox(double scale) {
  double theBox[3][3];
  int i, j;

  // cerr << "Scaling box by " << scale << "\n";

  for(i=0; i<3; i++) 
    for (j=0; j<3; j++) theBox[i][j] = Hmat[i][j]*scale;

  setBoxM(theBox);

}

void SimInfo::calcHmatInv( void ) {
  
  int oldOrtho;
  int i,j;
  double smallDiag;
  double tol;
  double sanity[3][3];

  invertMat3( Hmat, HmatInv );

  // check to see if Hmat is orthorhombic
  
  oldOrtho = orthoRhombic;

  smallDiag = fabs(Hmat[0][0]);
  if(smallDiag > fabs(Hmat[1][1])) smallDiag = fabs(Hmat[1][1]);
  if(smallDiag > fabs(Hmat[2][2])) smallDiag = fabs(Hmat[2][2]);
  tol = smallDiag * orthoTolerance;

  orthoRhombic = 1;
  
  for (i = 0; i < 3; i++ ) {
    for (j = 0 ; j < 3; j++) {
      if (i != j) {
        if (orthoRhombic) {
          if ( fabs(Hmat[i][j]) >= tol) orthoRhombic = 0;
        }        
      }
    }
  }

  if( oldOrtho != orthoRhombic ){
    
    if( orthoRhombic ) {
      sprintf( painCave.errMsg,
	       "OOPSE is switching from the default Non-Orthorhombic\n"
               "\tto the faster Orthorhombic periodic boundary computations.\n"
	       "\tThis is usually a good thing, but if you wan't the\n"
               "\tNon-Orthorhombic computations, make the orthoBoxTolerance\n"
               "\tvariable ( currently set to %G ) smaller.\n",
	       orthoTolerance);
      painCave.severity = OOPSE_INFO;
      simError();
    }
    else {
      sprintf( painCave.errMsg,
	       "OOPSE is switching from the faster Orthorhombic to the more\n"
               "\tflexible Non-Orthorhombic periodic boundary computations.\n"
	       "\tThis is usually because the box has deformed under\n"
               "\tNPTf integration. If you wan't to live on the edge with\n"
               "\tthe Orthorhombic computations, make the orthoBoxTolerance\n"
               "\tvariable ( currently set to %G ) larger.\n",
	       orthoTolerance);
      painCave.severity = OOPSE_WARNING;
      simError();
    }
  }
}

void SimInfo::calcBoxL( void ){

  double dx, dy, dz, dsq;

  // boxVol = Determinant of Hmat

  boxVol = matDet3( Hmat );

  // boxLx
  
  dx = Hmat[0][0]; dy = Hmat[1][0]; dz = Hmat[2][0];
  dsq = dx*dx + dy*dy + dz*dz;
  boxL[0] = sqrt( dsq );
  //maxCutoff = 0.5 * boxL[0];

  // boxLy
  
  dx = Hmat[0][1]; dy = Hmat[1][1]; dz = Hmat[2][1];
  dsq = dx*dx + dy*dy + dz*dz;
  boxL[1] = sqrt( dsq );
  //if( (0.5 * boxL[1]) < maxCutoff ) maxCutoff = 0.5 * boxL[1];


  // boxLz
  
  dx = Hmat[0][2]; dy = Hmat[1][2]; dz = Hmat[2][2];
  dsq = dx*dx + dy*dy + dz*dz;
  boxL[2] = sqrt( dsq );
  //if( (0.5 * boxL[2]) < maxCutoff ) maxCutoff = 0.5 * boxL[2];

  //calculate the max cutoff
  maxCutoff =  calcMaxCutOff(); 
  
  checkCutOffs();

}


double SimInfo::calcMaxCutOff(){

  double ri[3], rj[3], rk[3];
  double rij[3], rjk[3], rki[3];
  double minDist;

  ri[0] = Hmat[0][0];
  ri[1] = Hmat[1][0];
  ri[2] = Hmat[2][0];

  rj[0] = Hmat[0][1];
  rj[1] = Hmat[1][1];
  rj[2] = Hmat[2][1];

  rk[0] = Hmat[0][2];
  rk[1] = Hmat[1][2];
  rk[2] = Hmat[2][2];
    
  crossProduct3(ri, rj, rij);
  distXY = dotProduct3(rk,rij) / norm3(rij);

  crossProduct3(rj,rk, rjk);
  distYZ = dotProduct3(ri,rjk) / norm3(rjk);

  crossProduct3(rk,ri, rki);
  distZX = dotProduct3(rj,rki) / norm3(rki);

  minDist = min(min(distXY, distYZ), distZX);
  return minDist/2;
  
}

void SimInfo::wrapVector( double thePos[3] ){

  int i;
  double scaled[3];

  if( !orthoRhombic ){
    // calc the scaled coordinates.
  

    matVecMul3(HmatInv, thePos, scaled);
    
    for(i=0; i<3; i++)
      scaled[i] -= roundMe(scaled[i]);
    
    // calc the wrapped real coordinates from the wrapped scaled coordinates
    
    matVecMul3(Hmat, scaled, thePos);

  }
  else{
    // calc the scaled coordinates.
    
    for(i=0; i<3; i++)
      scaled[i] = thePos[i]*HmatInv[i][i];
    
    // wrap the scaled coordinates
    
    for(i=0; i<3; i++)
      scaled[i] -= roundMe(scaled[i]);
    
    // calc the wrapped real coordinates from the wrapped scaled coordinates
    
    for(i=0; i<3; i++)
      thePos[i] = scaled[i]*Hmat[i][i];
  }
    
}


int SimInfo::getNDF(){
  int ndf_local;

  ndf_local = 0;
  
  for(int i = 0; i < integrableObjects.size(); i++){
    ndf_local += 3;
    if (integrableObjects[i]->isDirectional()) {
      if (integrableObjects[i]->isLinear())
        ndf_local += 2;
      else
        ndf_local += 3;
    }
  }

  // n_constraints is local, so subtract them on each processor:

  ndf_local -= n_constraints;

#ifdef IS_MPI
  MPI_Allreduce(&ndf_local,&ndf,1,MPI_INT,MPI_SUM, MPI_COMM_WORLD);
#else
  ndf = ndf_local;
#endif

  // nZconstraints is global, as are the 3 COM translations for the 
  // entire system:

  ndf = ndf - 3 - nZconstraints;

  return ndf;
}

int SimInfo::getNDFraw() {
  int ndfRaw_local;

  // Raw degrees of freedom that we have to set
  ndfRaw_local = 0;

  for(int i = 0; i < integrableObjects.size(); i++){
    ndfRaw_local += 3;
    if (integrableObjects[i]->isDirectional()) {
       if (integrableObjects[i]->isLinear())
        ndfRaw_local += 2;
      else
        ndfRaw_local += 3;
    }
  }
    
#ifdef IS_MPI
  MPI_Allreduce(&ndfRaw_local,&ndfRaw,1,MPI_INT,MPI_SUM, MPI_COMM_WORLD);
#else
  ndfRaw = ndfRaw_local;
#endif

  return ndfRaw;
}

int SimInfo::getNDFtranslational() {
  int ndfTrans_local;

  ndfTrans_local = 3 * integrableObjects.size() - n_constraints;


#ifdef IS_MPI
  MPI_Allreduce(&ndfTrans_local,&ndfTrans,1,MPI_INT,MPI_SUM, MPI_COMM_WORLD);
#else
  ndfTrans = ndfTrans_local;
#endif

  ndfTrans = ndfTrans - 3 - nZconstraints;

  return ndfTrans;
}

int SimInfo::getTotIntegrableObjects() {
  int nObjs_local;
  int nObjs;

  nObjs_local =  integrableObjects.size();


#ifdef IS_MPI
  MPI_Allreduce(&nObjs_local,&nObjs,1,MPI_INT,MPI_SUM, MPI_COMM_WORLD);
#else
  nObjs = nObjs_local;
#endif


  return nObjs;
}

void SimInfo::refreshSim(){

  simtype fInfo;
  int isError;
  int n_global;
  int* excl;

  fInfo.dielect = 0.0;

  if( useDipoles ){
    if( useReactionField )fInfo.dielect = dielectric;
  }

  fInfo.SIM_uses_PBC = usePBC;

  if (useSticky || useDipoles || useGayBerne || useShapes) {
    useDirectionalAtoms = 1;
    fInfo.SIM_uses_DirectionalAtoms = useDirectionalAtoms;
  }

  fInfo.SIM_uses_LennardJones = useLennardJones;

  if (useCharges || useDipoles) {
    useElectrostatics = 1;
    fInfo.SIM_uses_Electrostatics = useElectrostatics;
  }

  fInfo.SIM_uses_Charges = useCharges;
  fInfo.SIM_uses_Dipoles = useDipoles;
  fInfo.SIM_uses_Sticky = useSticky;
  fInfo.SIM_uses_GayBerne = useGayBerne;
  fInfo.SIM_uses_EAM = useEAM;
  fInfo.SIM_uses_Shapes = useShapes;
  fInfo.SIM_uses_FLARB = useFLARB;
  fInfo.SIM_uses_RF = useReactionField;

  n_exclude = excludes->getSize();
  excl = excludes->getFortranArray();
  
#ifdef IS_MPI
  n_global = mpiSim->getNAtomsGlobal();
#else
  n_global = n_atoms;
#endif
  
  isError = 0;
  
  getFortranGroupArrays(this, FglobalGroupMembership, mfact);
  //it may not be a good idea to pass the address of first element in vector
  //since c++ standard does not require vector to be stored continuously in meomory
  //Most of the compilers will organize the memory of vector continuously
  setFortranSim( &fInfo, &n_global, &n_atoms, identArray, &n_exclude, excl, 
                  &nGlobalExcludes, globalExcludes, molMembershipArray, 
                  &mfact[0], &ngroup, &FglobalGroupMembership[0], &isError); 

  if( isError ){
    
    sprintf( painCave.errMsg,
             "There was an error setting the simulation information in fortran.\n" );
    painCave.isFatal = 1;
    painCave.severity = OOPSE_ERROR;
    simError();
  }
  
#ifdef IS_MPI
  sprintf( checkPointMsg,
	   "succesfully sent the simulation information to fortran.\n");
  MPIcheckPoint();
#endif // is_mpi
  
  this->ndf = this->getNDF();
  this->ndfRaw = this->getNDFraw();
  this->ndfTrans = this->getNDFtranslational();
}

void SimInfo::setDefaultRcut( double theRcut ){
  
  haveRcut = 1;
  rCut = theRcut;
  rList = rCut + 1.0; 
  
  notifyFortranCutoffs( &rCut, &rSw, &rList );
}

void SimInfo::setDefaultRcut( double theRcut, double theRsw ){

  rSw = theRsw;
  setDefaultRcut( theRcut );
}


void SimInfo::checkCutOffs( void ){
  
  if( boxIsInit ){
    
    //we need to check cutOffs against the box
    
    if( rCut > maxCutoff ){
      sprintf( painCave.errMsg,
	       "cutoffRadius is too large for the current periodic box.\n"
               "\tCurrent Value of cutoffRadius = %G at time %G\n "
               "\tThis is larger than half of at least one of the\n"
               "\tperiodic box vectors.  Right now, the Box matrix is:\n"
	       "\n"
	       "\t[ %G %G %G ]\n"
	       "\t[ %G %G %G ]\n"
	       "\t[ %G %G %G ]\n",
	       rCut, currentTime,
	       Hmat[0][0], Hmat[0][1], Hmat[0][2],
	       Hmat[1][0], Hmat[1][1], Hmat[1][2],
	       Hmat[2][0], Hmat[2][1], Hmat[2][2]);
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();
    }    
  } else {
    // initialize this stuff before using it, OK?
    sprintf( painCave.errMsg,
             "Trying to check cutoffs without a box.\n"
             "\tOOPSE should have better programmers than that.\n" );
    painCave.severity = OOPSE_ERROR;
    painCave.isFatal = 1;
    simError();      
  }
  
}

void SimInfo::addProperty(GenericData* prop){

  map<string, GenericData*>::iterator result;
  result = properties.find(prop->getID());
  
  //we can't simply use  properties[prop->getID()] = prop,
  //it will cause memory leak if we already contain a propery which has the same name of prop
  
  if(result != properties.end()){
    
    delete (*result).second;
    (*result).second = prop;
      
  }
  else{

    properties[prop->getID()] = prop;

  }
    
}

GenericData* SimInfo::getProperty(const string& propName){
 
  map<string, GenericData*>::iterator result;
  
  //string lowerCaseName = ();
  
  result = properties.find(propName);
  
  if(result != properties.end()) 
    return (*result).second;  
  else   
    return NULL;  
}


void SimInfo::getFortranGroupArrays(SimInfo* info, 
                                    vector<int>& FglobalGroupMembership,
                                    vector<double>& mfact){
  
  Molecule* myMols;
  Atom** myAtoms;
  int numAtom;
  double mtot;
  int numMol;
  int numCutoffGroups;
  CutoffGroup* myCutoffGroup;
  vector<CutoffGroup*>::iterator iterCutoff;
  Atom* cutoffAtom;
  vector<Atom*>::iterator iterAtom;
  int atomIndex;
  double totalMass;
  
  mfact.clear();
  FglobalGroupMembership.clear();
  

  // Fix the silly fortran indexing problem
#ifdef IS_MPI
  numAtom = mpiSim->getNAtomsGlobal();
#else
  numAtom = n_atoms;
#endif
  for (int i = 0; i < numAtom; i++) 
    FglobalGroupMembership.push_back(globalGroupMembership[i] + 1);
  

  myMols = info->molecules;
  numMol = info->n_mol;
  for(int i  = 0; i < numMol; i++){
    numCutoffGroups = myMols[i].getNCutoffGroups();
    for(myCutoffGroup =myMols[i].beginCutoffGroup(iterCutoff); 
        myCutoffGroup != NULL; 
        myCutoffGroup =myMols[i].nextCutoffGroup(iterCutoff)){

      totalMass = myCutoffGroup->getMass();
      
      for(cutoffAtom = myCutoffGroup->beginAtom(iterAtom); 
          cutoffAtom != NULL; 
          cutoffAtom = myCutoffGroup->nextAtom(iterAtom)){
        mfact.push_back(cutoffAtom->getMass()/totalMass);
      }  
    }
  }

}
