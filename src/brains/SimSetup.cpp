#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <string>
#include <sprng.h> 
#include "brains/SimSetup.hpp"
#include "io/ReadWrite.hpp"
#include "io/parse_me.h"
#include "integrators/Integrator.hpp"
#include "utils/simError.h"
#include "primitives/RigidBody.hpp"
#include "minimizers/OOPSEMinimizer.hpp"

#ifdef IS_MPI
#include "io/mpiBASS.h"
#include "brains/mpiSimulation.hpp"
#endif

// some defines for ensemble and Forcefield  cases

#define NVE_ENS        0 
#define NVT_ENS        1
#define NPTi_ENS       2
#define NPTf_ENS       3
#define NPTxyz_ENS     4


#define FF_DUFF  0
#define FF_LJ    1
#define FF_EAM   2
#define FF_H2O   3

using namespace std;
using namespace oopse;

/**
 * Check whether dividend is divisble by divisor or not
 */
bool isDivisible(double dividend, double divisor){
  double tolerance = 0.000001;
  double quotient;
  double diff;
  int intQuotient;
  
  quotient = dividend / divisor;

  if (quotient < 0)
    quotient = -quotient;

  intQuotient = int (quotient + tolerance);

  diff = fabs(fabs(dividend) - intQuotient  * fabs(divisor));

  if (diff <= tolerance)
    return true;
  else
    return false;  
}

string getPrefix(const string& str ){
  string prefix;
  string suffix;
  int pos; 

  pos = str.rfind(".");

  if (pos >= 0) {
     prefix = str.substr(0, pos);
     suffix = str.substr(pos, str.size());

     // leave .bass there in case we've reverted to old habits 
     if (LowerCase(suffix) == ".md" || LowerCase(suffix) == ".bass")
       return prefix;
     else 
       return str;
     
  } else 
    return str;
};


SimSetup::SimSetup(){
  
  initSuspend = false;
  isInfoArray = 0;
  nInfo = 1;

  stamps = new MakeStamps();
  globals = new Globals();


#ifdef IS_MPI
  strcpy(checkPointMsg, "SimSetup creation successful");
  MPIcheckPoint();
#endif // IS_MPI
}

SimSetup::~SimSetup(){
  // clean up the forcefield
  the_ff->cleanMe();

  delete stamps;
  delete globals;
}

void SimSetup::setSimInfo(SimInfo* the_info, int theNinfo){
  info = the_info;
  nInfo = theNinfo;
  isInfoArray = 1;
  initSuspend = true;
}


void SimSetup::parseFile(char* fileName){
#ifdef IS_MPI
  if (worldRank == 0){
#endif // is_mpi

    inFileName = fileName;

    globals->initalize();
    set_interface_stamps(stamps, globals);

#ifdef IS_MPI
    mpiEventInit();
#endif

    yacc_BASS(fileName);

#ifdef IS_MPI
    throwMPIEvent(NULL);
  }
  else{
    receiveParse();
  }
#endif

}

#ifdef IS_MPI
void SimSetup::receiveParse(void){
  set_interface_stamps(stamps, globals);
  mpiEventInit();
  MPIcheckPoint();
  mpiEventLoop();
}

#endif // is_mpi

void SimSetup::createSim(void){

  // gather all of the information from the meta-data file

  gatherInfo();

  // creation of complex system objects

  sysObjectsCreation();

  // check on the post processing info

  finalInfoCheck();

  // initialize the system coordinates

  if ( !initSuspend ){
    initSystemCoords();

    if( !(globals->getUseInitTime()) )
      info[0].currentTime = 0.0;
  }  

  // make the output filenames

  makeOutNames();
  
#ifdef IS_MPI
  mpiSim->mpiRefresh();
#endif

  // initialize the Fortran

  initFortran();

  if (globals->haveMinimizer())
    // make minimizer
    makeMinimizer();
  else
    // make the integrator
    makeIntegrator();

}


void SimSetup::makeMolecules(void){
  int i, j, k;
  int exI, exJ, exK, exL, slI, slJ;
  int tempI, tempJ, tempK, tempL;
  int molI, globalID;
  int stampID, atomOffset, rbOffset, groupOffset;
  molInit molInfo;
  DirectionalAtom* dAtom;
  RigidBody* myRB;
  StuntDouble* mySD;
  LinkedAssign* extras;
  LinkedAssign* current_extra;
  AtomStamp* currentAtom;
  BondStamp* currentBond;
  BendStamp* currentBend;
  TorsionStamp* currentTorsion;
  RigidBodyStamp* currentRigidBody;
  CutoffGroupStamp* currentCutoffGroup;
  CutoffGroup* myCutoffGroup;
  int nCutoffGroups;// number of cutoff group of a molecule defined in mdl file
  set<int> cutoffAtomSet; //atoms belong to  cutoffgroup defined at mdl file

  bond_pair* theBonds;
  bend_set* theBends;
  torsion_set* theTorsions;

  set<int> skipList;

  double phi, theta, psi;
  char* molName;
  char rbName[100];

  int whichRigidBody; 
  int consAtomIndex;  //index of constraint atom in rigid body's atom array
  double bondLength2;
  //init the forceField paramters

  the_ff->readParams();

  // init the atoms

  int nMembers, nNew, rb1, rb2;

  for (k = 0; k < nInfo; k++){
    the_ff->setSimInfo(&(info[k]));

#ifdef IS_MPI
    info[k].globalGroupMembership = new int[mpiSim->getNAtomsGlobal()];
    for (i = 0; i < mpiSim->getNAtomsGlobal(); i++) 
      info[k].globalGroupMembership[i] = 0;
#else
    info[k].globalGroupMembership = new int[info[k].n_atoms];
    for (i = 0; i < info[k].n_atoms; i++) 
      info[k].globalGroupMembership[i] = 0;
#endif

    atomOffset = 0;
    groupOffset = 0;

    for (i = 0; i < info[k].n_mol; i++){
      stampID = info[k].molecules[i].getStampID();
      molName = comp_stamps[stampID]->getID();

      molInfo.nAtoms = comp_stamps[stampID]->getNAtoms();
      molInfo.nBonds = comp_stamps[stampID]->getNBonds();
      molInfo.nBends = comp_stamps[stampID]->getNBends();
      molInfo.nTorsions = comp_stamps[stampID]->getNTorsions();
      molInfo.nRigidBodies = comp_stamps[stampID]->getNRigidBodies();

      nCutoffGroups = comp_stamps[stampID]->getNCutoffGroups();
      
      molInfo.myAtoms = &(info[k].atoms[atomOffset]);

      if (molInfo.nBonds > 0) 
        molInfo.myBonds = new Bond*[molInfo.nBonds];
      else 
        molInfo.myBonds = NULL;

      if (molInfo.nBends > 0) 
        molInfo.myBends = new Bend*[molInfo.nBends];
      else 
        molInfo.myBends = NULL;

      if (molInfo.nTorsions > 0) 
        molInfo.myTorsions = new Torsion *[molInfo.nTorsions];
      else 
        molInfo.myTorsions = NULL;

      theBonds = new bond_pair[molInfo.nBonds];
      theBends = new bend_set[molInfo.nBends];
      theTorsions = new torsion_set[molInfo.nTorsions];
      
      // make the Atoms

      for (j = 0; j < molInfo.nAtoms; j++){
        currentAtom = comp_stamps[stampID]->getAtom(j);

        if (currentAtom->haveOrientation()){
          dAtom = new DirectionalAtom((j + atomOffset),
                                      info[k].getConfiguration());
          info[k].n_oriented++;
          molInfo.myAtoms[j] = dAtom;

          // Directional Atoms have standard unit vectors which are oriented
          // in space using the three Euler angles.  We assume the standard
          // unit vector was originally along the z axis below.

          phi = currentAtom->getEulerPhi() * M_PI / 180.0;
          theta = currentAtom->getEulerTheta() * M_PI / 180.0;
          psi = currentAtom->getEulerPsi()* M_PI / 180.0;

          dAtom->setUnitFrameFromEuler(phi, theta, psi);
            
        }
        else{

          molInfo.myAtoms[j] = new Atom((j + atomOffset), info[k].getConfiguration());

        }

        molInfo.myAtoms[j]->setType(currentAtom->getType());
#ifdef IS_MPI
        molInfo.myAtoms[j]->setGlobalIndex(globalAtomIndex[j + atomOffset]);
#endif // is_mpi
      } 

      // make the bonds
      for (j = 0; j < molInfo.nBonds; j++){
        currentBond = comp_stamps[stampID]->getBond(j);
        theBonds[j].a = currentBond->getA() + atomOffset;
        theBonds[j].b = currentBond->getB() + atomOffset;

        tempI = theBonds[j].a;
        tempJ = theBonds[j].b;

#ifdef IS_MPI
        exI = info[k].atoms[tempI]->getGlobalIndex() + 1;
        exJ = info[k].atoms[tempJ]->getGlobalIndex() + 1;
#else
        exI = tempI + 1;
        exJ = tempJ + 1;
#endif

        info[k].excludes->addPair(exI, exJ);
      }

      //make the bends
      for (j = 0; j < molInfo.nBends; j++){
        currentBend = comp_stamps[stampID]->getBend(j);
        theBends[j].a = currentBend->getA() + atomOffset;
        theBends[j].b = currentBend->getB() + atomOffset;
        theBends[j].c = currentBend->getC() + atomOffset;

        if (currentBend->haveExtras()){
          extras = currentBend->getExtras();
          current_extra = extras;

          while (current_extra != NULL){
            if (!strcmp(current_extra->getlhs(), "ghostVectorSource")){
              switch (current_extra->getType()){
                case 0:
                  theBends[j].ghost = current_extra->getInt() + atomOffset;
                  theBends[j].isGhost = 1;
                  break;

                case 1:
                  theBends[j].ghost = (int) current_extra->getDouble() +
                                      atomOffset;
                  theBends[j].isGhost = 1;
                  break;

                default:
                  sprintf(painCave.errMsg,
                          "SimSetup Error: ghostVectorSource was neither a "
                          "double nor an int.\n"
                          "-->Bend[%d] in %s\n",
                          j, comp_stamps[stampID]->getID());
                  painCave.isFatal = 1;
                  simError();
              }
            }
            else{
              sprintf(painCave.errMsg,
                      "SimSetup Error: unhandled bend assignment:\n"
                      "    -->%s in Bend[%d] in %s\n",
                      current_extra->getlhs(), j, comp_stamps[stampID]->getID());
              painCave.isFatal = 1;
              simError();
            }

            current_extra = current_extra->getNext();
          }
        }

        if (theBends[j].isGhost) {
          
          tempI = theBends[j].a;
          tempJ = theBends[j].b;
          
#ifdef IS_MPI
          exI = info[k].atoms[tempI]->getGlobalIndex() + 1;
          exJ = info[k].atoms[tempJ]->getGlobalIndex() + 1;
#else
          exI = tempI + 1;
          exJ = tempJ + 1;
#endif          
          info[k].excludes->addPair(exI, exJ);

        } else {

          tempI = theBends[j].a;
          tempJ = theBends[j].b;
          tempK = theBends[j].c;
          
#ifdef IS_MPI
          exI = info[k].atoms[tempI]->getGlobalIndex() + 1;
          exJ = info[k].atoms[tempJ]->getGlobalIndex() + 1;
          exK = info[k].atoms[tempK]->getGlobalIndex() + 1;
#else
          exI = tempI + 1;
          exJ = tempJ + 1;
          exK = tempK + 1;
#endif
          
          info[k].excludes->addPair(exI, exK);
          info[k].excludes->addPair(exI, exJ);
          info[k].excludes->addPair(exJ, exK);
        }
      }

      for (j = 0; j < molInfo.nTorsions; j++){
        currentTorsion = comp_stamps[stampID]->getTorsion(j);
        theTorsions[j].a = currentTorsion->getA() + atomOffset;
        theTorsions[j].b = currentTorsion->getB() + atomOffset;
        theTorsions[j].c = currentTorsion->getC() + atomOffset;
        theTorsions[j].d = currentTorsion->getD() + atomOffset;

        tempI = theTorsions[j].a;       
        tempJ = theTorsions[j].b;
        tempK = theTorsions[j].c;
        tempL = theTorsions[j].d;

#ifdef IS_MPI
        exI = info[k].atoms[tempI]->getGlobalIndex() + 1;
        exJ = info[k].atoms[tempJ]->getGlobalIndex() + 1;
        exK = info[k].atoms[tempK]->getGlobalIndex() + 1;
        exL = info[k].atoms[tempL]->getGlobalIndex() + 1;
#else
        exI = tempI + 1;
        exJ = tempJ + 1;
        exK = tempK + 1;
        exL = tempL + 1;
#endif

        info[k].excludes->addPair(exI, exJ);
        info[k].excludes->addPair(exI, exK);
        info[k].excludes->addPair(exI, exL);        
        info[k].excludes->addPair(exJ, exK);
        info[k].excludes->addPair(exJ, exL);
        info[k].excludes->addPair(exK, exL);
      }

      
      molInfo.myRigidBodies.clear();
      
      for (j = 0; j < molInfo.nRigidBodies; j++){

        currentRigidBody = comp_stamps[stampID]->getRigidBody(j);
        nMembers = currentRigidBody->getNMembers(); 

        // Create the Rigid Body:

        myRB = new RigidBody();

        sprintf(rbName,"%s_RB_%d", molName, j);
        myRB->setType(rbName);
        
        for (rb1 = 0; rb1 < nMembers; rb1++) {

          // molI is atom numbering inside this molecule
          molI = currentRigidBody->getMember(rb1);     

          // tempI is atom numbering on local processor
          tempI = molI + atomOffset;

          // currentAtom is the AtomStamp (which we need for 
          // rigid body reference positions)
          currentAtom = comp_stamps[stampID]->getAtom(molI);

          // When we add to the rigid body, add the atom itself and 
          // the stamp info:

          myRB->addAtom(info[k].atoms[tempI], currentAtom);
          
          // Add this atom to the Skip List for the integrators
#ifdef IS_MPI
          slI = info[k].atoms[tempI]->getGlobalIndex();
#else
          slI = tempI;
#endif
          skipList.insert(slI);
          
        }
        
        for(rb1 = 0; rb1 < nMembers - 1; rb1++) {
          for(rb2 = rb1+1; rb2 < nMembers; rb2++) {
            
            tempI = currentRigidBody->getMember(rb1);
            tempJ = currentRigidBody->getMember(rb2);
            
            // Some explanation is required here.
            // Fortran indexing starts at 1, while c indexing starts at 0
            // Also, in parallel computations, the GlobalIndex is
            // used for the exclude list:
            
#ifdef IS_MPI
            exI = molInfo.myAtoms[tempI]->getGlobalIndex() + 1;
            exJ = molInfo.myAtoms[tempJ]->getGlobalIndex() + 1;
#else
            exI = molInfo.myAtoms[tempI]->getIndex() + 1;
            exJ = molInfo.myAtoms[tempJ]->getIndex() + 1;
#endif
            
            info[k].excludes->addPair(exI, exJ);
            
          }
        }

        molInfo.myRigidBodies.push_back(myRB);
        info[k].rigidBodies.push_back(myRB);
      }
      

      //create cutoff group for molecule

      cutoffAtomSet.clear();
      molInfo.myCutoffGroups.clear();
      
      for (j = 0; j < nCutoffGroups; j++){

        currentCutoffGroup = comp_stamps[stampID]->getCutoffGroup(j);
        nMembers = currentCutoffGroup->getNMembers(); 

        myCutoffGroup = new CutoffGroup();
        
#ifdef IS_MPI
        myCutoffGroup->setGlobalIndex(globalGroupIndex[groupOffset]);
#else
        myCutoffGroup->setGlobalIndex(groupOffset);
#endif
        
        for (int cg = 0; cg < nMembers; cg++) {

          // molI is atom numbering inside this molecule
          molI = currentCutoffGroup->getMember(cg);     

          // tempI is atom numbering on local processor
          tempI = molI + atomOffset;

#ifdef IS_MPI
          globalID = info[k].atoms[tempI]->getGlobalIndex();
          info[k].globalGroupMembership[globalID] = globalGroupIndex[groupOffset];
#else 
          globalID = info[k].atoms[tempI]->getIndex();
          info[k].globalGroupMembership[globalID] = groupOffset;
#endif                    
          myCutoffGroup->addAtom(info[k].atoms[tempI]);
          cutoffAtomSet.insert(tempI);
        }
        
        molInfo.myCutoffGroups.push_back(myCutoffGroup);
        groupOffset++;

      }//end for (j = 0; j < molInfo.nCutoffGroups; j++)
      
      
      // create a cutoff group for every atom in current molecule which
      // does not belong to cutoffgroup defined at mdl file
      
      for(j = 0; j < molInfo.nAtoms; j++){
        
        if(cutoffAtomSet.find(molInfo.myAtoms[j]->getIndex()) == cutoffAtomSet.end()){
          myCutoffGroup = new CutoffGroup();
          myCutoffGroup->addAtom(molInfo.myAtoms[j]);
          
#ifdef IS_MPI
          myCutoffGroup->setGlobalIndex(globalGroupIndex[groupOffset]);
          globalID = info[k].atoms[atomOffset + j]->getGlobalIndex();
          info[k].globalGroupMembership[globalID] = globalGroupIndex[groupOffset]; 
#else
          myCutoffGroup->setGlobalIndex(groupOffset);
          globalID = info[k].atoms[atomOffset + j]->getIndex();
          info[k].globalGroupMembership[globalID] = groupOffset;
#endif
          molInfo.myCutoffGroups.push_back(myCutoffGroup);
          groupOffset++;
        }          
      }

      // After this is all set up, scan through the atoms to 
      // see if they can be added to the integrableObjects:

      molInfo.myIntegrableObjects.clear();
      

      for (j = 0; j < molInfo.nAtoms; j++){

#ifdef IS_MPI
        slJ = molInfo.myAtoms[j]->getGlobalIndex();
#else
        slJ = j+atomOffset;
#endif

        // if they aren't on the skip list, then they can be integrated

        if (skipList.find(slJ) == skipList.end()) {
          mySD = (StuntDouble *) molInfo.myAtoms[j];
          info[k].integrableObjects.push_back(mySD);
          molInfo.myIntegrableObjects.push_back(mySD);
        }
      }

      // all rigid bodies are integrated:

      for (j = 0; j < molInfo.nRigidBodies; j++) {
        mySD = (StuntDouble *) molInfo.myRigidBodies[j];
        info[k].integrableObjects.push_back(mySD);      
        molInfo.myIntegrableObjects.push_back(mySD);
      }
         
      // send the arrays off to the forceField for init.
      
      the_ff->initializeAtoms(molInfo.nAtoms, molInfo.myAtoms);
      the_ff->initializeBonds(molInfo.nBonds, molInfo.myBonds, theBonds);
      the_ff->initializeBends(molInfo.nBends, molInfo.myBends, theBends);
      the_ff->initializeTorsions(molInfo.nTorsions, molInfo.myTorsions,
                                 theTorsions);

      info[k].molecules[i].initialize(molInfo);
      
      
      atomOffset += molInfo.nAtoms;
      delete[] theBonds;
      delete[] theBends;
      delete[] theTorsions;
    }



#ifdef IS_MPI    
    // Since the globalGroupMembership has been zero filled and we've only
    // poked values into the atoms we know, we can do an Allreduce
    // to get the full globalGroupMembership array (We think).
    // This would be prettier if we could use MPI_IN_PLACE like the MPI-2
    // docs said we could.

    int* ggMjunk = new int[mpiSim->getNAtomsGlobal()];    

    MPI_Allreduce(info[k].globalGroupMembership,
                  ggMjunk,
                  mpiSim->getNAtomsGlobal(),
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    for (i = 0; i < mpiSim->getNAtomsGlobal(); i++) 
      info[k].globalGroupMembership[i] = ggMjunk[i];

    delete[] ggMjunk;
    
#endif



  }

#ifdef IS_MPI
  sprintf(checkPointMsg, "all molecules initialized succesfully");
  MPIcheckPoint();
#endif // is_mpi

}

void SimSetup::gatherInfo(void){
  int i;

  ensembleCase = -1;
  ffCase = -1;

  // set the easy ones first

  for (i = 0; i < nInfo; i++){
    if (globals->haveTargetTemp()) {
      info[i].target_temp = globals->getTargetTemp();
      info[i].have_target_temp = 1;
    } else {
      info[i].have_target_temp = 0;
    }
    if (globals->haveDt()) {
      info[i].dt = globals->getDt();
    }
    if (globals->haveRunTime()) {
      info[i].run_time = globals->getRunTime();
    }
    }
  n_components = globals->getNComponents();


  // get the forceField

  strcpy(force_field, globals->getForceField());

  if (!strcasecmp(force_field, "DUFF")){
    ffCase = FF_DUFF;
  }
  else if (!strcasecmp(force_field, "LJ")){
    ffCase = FF_LJ;
  }
  else if (!strcasecmp(force_field, "EAM")){
    ffCase = FF_EAM;
  }
  else if (!strcasecmp(force_field, "WATER")){
    ffCase = FF_H2O;
  }
  else{
    sprintf(painCave.errMsg, "SimSetup Error. Unrecognized force field -> %s\n",
            force_field);
         painCave.isFatal = 1;
         simError();
  }
  if (globals->haveForceFieldVariant()) {
    strcpy(forcefield_variant, globals->getForceFieldVariant());
    has_forcefield_variant = 1;
  }
  
  // get the ensemble


  if (globals->haveEnsemble()) {
    
    strcpy(ensemble, globals->getEnsemble());
    
    if (!strcasecmp(ensemble, "NVE")){
      ensembleCase = NVE_ENS;
    }
    else if (!strcasecmp(ensemble, "NVT")){
      ensembleCase = NVT_ENS;
    }
    else if (!strcasecmp(ensemble, "NPTi") || !strcasecmp(ensemble, "NPT")){
      ensembleCase = NPTi_ENS;
    }
    else if (!strcasecmp(ensemble, "NPTf")){
      ensembleCase = NPTf_ENS;
    }
    else if (!strcasecmp(ensemble, "NPTxyz")){
      ensembleCase = NPTxyz_ENS;
    }
    else{
      sprintf(painCave.errMsg,
              "SimSetup Warning. Unrecognized Ensemble -> %s \n"
              "\treverting to NVE for this simulation.\n",
              ensemble);
      painCave.isFatal = 0;
      simError();
      strcpy(ensemble, "NVE");
      ensembleCase = NVE_ENS;
    }  
    
    for (i = 0; i < nInfo; i++) 
      strcpy(info[i].ensemble, ensemble);
      
    
    //check whether sample time, status time, thermal time and reset time are divisble by dt
    if (globals->haveSampleTime() && !isDivisible(globals->getSampleTime(), globals->getDt())){
      sprintf(painCave.errMsg,
              "Sample time is not divisible by dt.\n"
              "\tThis will result in samples that are not uniformly\n"
              "\tdistributed in time.  If this is a problem, change\n"
              "\tyour sampleTime variable.\n");
      painCave.isFatal = 0;
      simError();    
    }
    
    if (globals->haveStatusTime() && !isDivisible(globals->getStatusTime(), globals->getDt())){
      sprintf(painCave.errMsg,
              "Status time is not divisible by dt.\n"
              "\tThis will result in status reports that are not uniformly\n"
              "\tdistributed in time.  If this is a problem, change \n"
              "\tyour statusTime variable.\n");
      painCave.isFatal = 0;
      simError();    
    }
    
    if (globals->haveThermalTime() && !isDivisible(globals->getThermalTime(), globals->getDt())){
      sprintf(painCave.errMsg,
              "Thermal time is not divisible by dt.\n"
              "\tThis will result in thermalizations that are not uniformly\n"
              "\tdistributed in time.  If this is a problem, change \n"
              "\tyour thermalTime variable.\n");
      painCave.isFatal = 0;
      simError();    
    }  
    
    if (globals->haveResetTime() && !isDivisible(globals->getResetTime(), globals->getDt())){
      sprintf(painCave.errMsg,
              "Reset time is not divisible by dt.\n"
              "\tThis will result in integrator resets that are not uniformly\n"
              "\tdistributed in time.  If this is a problem, change\n"
              "\tyour resetTime variable.\n");
      painCave.isFatal = 0;
      simError();    
    } 
    
    // set the status, sample, and thermal kick times
    
    for (i = 0; i < nInfo; i++){
      if (globals->haveSampleTime()){
        info[i].sampleTime = globals->getSampleTime();
        info[i].statusTime = info[i].sampleTime;
      }
      else{
        info[i].sampleTime = globals->getRunTime();
        info[i].statusTime = info[i].sampleTime;
      }
      
      if (globals->haveStatusTime()){
        info[i].statusTime = globals->getStatusTime();
      }
      
      if (globals->haveThermalTime()){
        info[i].thermalTime = globals->getThermalTime();
      } else {
        info[i].thermalTime = globals->getRunTime();
      }
      
      info[i].resetIntegrator = 0;
      if( globals->haveResetTime() ){
        info[i].resetTime = globals->getResetTime();
        info[i].resetIntegrator = 1;
      }        
    }

    for (i=0; i < nInfo; i++) {
      
      // check for the temperature set flag
      
      if (globals->haveTempSet())
        info[i].setTemp = globals->getTempSet();
      
      // check for the extended State init
      
      info[i].useInitXSstate = globals->getUseInitXSstate();
      info[i].orthoTolerance = globals->getOrthoBoxTolerance();
      
      // check for thermodynamic integration
      if (globals->getUseSolidThermInt() && !globals->getUseLiquidThermInt()) {
        if (globals->haveThermIntLambda() && globals->haveThermIntK()) {
          info[i].useSolidThermInt = globals->getUseSolidThermInt();
          info[i].thermIntLambda = globals->getThermIntLambda();
          info[i].thermIntK = globals->getThermIntK();
          
          Restraints *myRestraint = new Restraints(info[i].thermIntLambda, info[i].thermIntK);
          info[i].restraint = myRestraint;
        }
        else {
          sprintf(painCave.errMsg,
                  "SimSetup Error:\n"
                  "\tKeyword useSolidThermInt was set to 'true' but\n"
                  "\tthermodynamicIntegrationLambda (and/or\n"
                  "\tthermodynamicIntegrationK) was not specified.\n"
                  "\tPlease provide a lambda value and k value in your meta-data file.\n");
          painCave.isFatal = 1;
          simError();    
        }
      }
      else if(globals->getUseLiquidThermInt()) {
        if (globals->getUseSolidThermInt()) {
          sprintf( painCave.errMsg,
                   "SimSetup Warning: It appears that you have both solid and\n"
                   "\tliquid thermodynamic integration activated in your meta-data\n"
                   "\tfile. To avoid confusion, specify only one technique in\n"
                   "\tyour meta-data file. Liquid-state thermodynamic integration\n"
                   "\twill be assumed for the current simulation. If this is not\n"
                   "\twhat you desire, set useSolidThermInt to 'true' and\n"
                   "\tuseLiquidThermInt to 'false' in your meta-data file.\n");
          painCave.isFatal = 0;
          simError();
        }
        if (globals->haveThermIntLambda() && globals->haveThermIntK()) {
          info[i].useLiquidThermInt = globals->getUseLiquidThermInt();
          info[i].thermIntLambda = globals->getThermIntLambda();
          info[i].thermIntK = globals->getThermIntK();
        }
        else {
          sprintf(painCave.errMsg,
                  "SimSetup Error:\n"
                  "\tKeyword useLiquidThermInt was set to 'true' but\n"
                  "\tthermodynamicIntegrationLambda (and/or\n"
                  "\tthermodynamicIntegrationK) was not specified.\n"
                  "\tPlease provide a lambda value and k value in your meta-data file.\n");
          painCave.isFatal = 1;
          simError();    
        }
      }
      else if(globals->haveThermIntLambda() || globals->haveThermIntK()){
        sprintf(painCave.errMsg,
                "SimSetup Warning: If you want to use Thermodynamic\n"
                "\tIntegration, set useSolidThermInt or useLiquidThermInt to\n"
                "\t'true' in your meta-data file.  These keywords are set to\n"
                "\t'false' by default, so your lambda and/or k values are\n"
                "\tbeing ignored.\n");
        painCave.isFatal = 0;
        simError();   
      }
    }        
  }
  
  for (i = 0; i < nInfo; i++) {
    info[i].usePBC = globals->getPBC();
  }
  
  // get the components and calculate the tot_nMol and indvidual n_mol
  
  the_components = globals->getComponents();
  components_nmol = new int[n_components];
  
  if (!globals->haveNMol()){
    // we don't have the total number of molecules, so we assume it is
    // given in each component

    tot_nmol = 0;
    for (i = 0; i < n_components; i++){
      if (!the_components[i]->haveNMol()){
        // we have a problem
        sprintf(painCave.errMsg,
                "SimSetup Error. No global NMol or component NMol given.\n"
                "\tCannot calculate the number of atoms.\n");
        painCave.isFatal = 1;
        simError();
      }

      tot_nmol += the_components[i]->getNMol();
      components_nmol[i] = the_components[i]->getNMol();
    }
  }
  else{
    sprintf(painCave.errMsg,
            "SimSetup error.\n"
            "\tSorry, the ability to specify total"
            " nMols and then give molfractions in the components\n"
            "\tis not currently supported."
            " Please give nMol in the components.\n");
    painCave.isFatal = 1;
    simError();
  }
  


  
  //setup seed for random number generator
  int seedValue;

  if (globals->haveSeed()){
    seedValue = globals->getSeed();

    if(seedValue / 1E9 == 0){
      sprintf(painCave.errMsg,
              "Seed for sprng library should contain at least 9 digits\n"
              "OOPSE will generate a seed for user\n");
      painCave.isFatal = 0;
      simError();

      //using seed generated by system instead of invalid seed set by user 
#ifndef IS_MPI
      seedValue = make_sprng_seed();
#else 
      if (worldRank == 0){
        seedValue = make_sprng_seed();
      }
      MPI_Bcast(&seedValue, 1, MPI_INT, 0, MPI_COMM_WORLD);   
#endif      
    }
  }//end of if branch of globals->haveSeed()
  else{
    
#ifndef IS_MPI
    seedValue = make_sprng_seed();
#else 
    if (worldRank == 0){
      seedValue = make_sprng_seed();
    }
    MPI_Bcast(&seedValue, 1, MPI_INT, 0, MPI_COMM_WORLD);   
#endif
  }//end of globals->haveSeed()

  for (int i = 0; i < nInfo; i++){
    info[i].setSeed(seedValue);
  }
  
#ifdef IS_MPI
  strcpy(checkPointMsg, "Successfully gathered all information from meta-data file\n");
  MPIcheckPoint();
#endif // is_mpi
}


void SimSetup::finalInfoCheck(void){
  int index;
  int usesDipoles;
  int usesCharges;
  int i;

  for (i = 0; i < nInfo; i++){
    // check electrostatic parameters

    index = 0;
    usesDipoles = 0;
    while ((index < info[i].n_atoms) && !usesDipoles){
      usesDipoles = (info[i].atoms[index])->hasDipole();
      index++;
    }
    index = 0;
    usesCharges = 0;
    while ((index < info[i].n_atoms) && !usesCharges){
      usesCharges= (info[i].atoms[index])->hasCharge();
      index++;
    }
#ifdef IS_MPI
    int myUse = usesDipoles;
    MPI_Allreduce(&myUse, &usesDipoles, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
#endif //is_mpi

    double theRcut, theRsw;

    if (globals->haveRcut()) {
      theRcut = globals->getRcut();

      if (globals->haveRsw()) 
        theRsw = globals->getRsw();
      else 
        theRsw = theRcut;
      
      info[i].setDefaultRcut(theRcut, theRsw);

    } else {
      
      the_ff->calcRcut();
      theRcut = info[i].getRcut();

      if (globals->haveRsw()) 
        theRsw = globals->getRsw();
      else 
        theRsw = theRcut;
      
      info[i].setDefaultRcut(theRcut, theRsw);
    }

    if (globals->getUseRF()){
      info[i].useReactionField = 1;
      
      if (!globals->haveRcut()){
        sprintf(painCave.errMsg,
                "SimSetup Warning: No value was set for the cutoffRadius.\n"
                "\tOOPSE will use a default value of 15.0 angstroms"
                "\tfor the cutoffRadius.\n");
        painCave.isFatal = 0;
        simError();
        theRcut = 15.0;
      }
      else{
        theRcut = globals->getRcut();
      }

      if (!globals->haveRsw()){
        sprintf(painCave.errMsg,
                "SimSetup Warning: No value was set for switchingRadius.\n"
                "\tOOPSE will use a default value of\n"
                "\t0.95 * cutoffRadius for the switchingRadius\n");
        painCave.isFatal = 0;
        simError();
        theRsw = 0.95 * theRcut;
      }
      else{
        theRsw = globals->getRsw();
      }

      info[i].setDefaultRcut(theRcut, theRsw);

      if (!globals->haveDielectric()){
        sprintf(painCave.errMsg,
                "SimSetup Error: No Dielectric constant was set.\n"
                "\tYou are trying to use Reaction Field without"
                "\tsetting a dielectric constant!\n");
        painCave.isFatal = 1;
        simError();
      }
      info[i].dielectric = globals->getDielectric();
    }
    else{
      if (usesDipoles || usesCharges){

        if (!globals->haveRcut()){
          sprintf(painCave.errMsg,
                  "SimSetup Warning: No value was set for the cutoffRadius.\n"
                  "\tOOPSE will use a default value of 15.0 angstroms"
                  "\tfor the cutoffRadius.\n");
          painCave.isFatal = 0;
          simError();
          theRcut = 15.0;
      }
        else{
          theRcut = globals->getRcut();
        }
        
        if (!globals->haveRsw()){
          sprintf(painCave.errMsg,
                  "SimSetup Warning: No value was set for switchingRadius.\n"
                  "\tOOPSE will use a default value of\n"
                  "\t0.95 * cutoffRadius for the switchingRadius\n");
          painCave.isFatal = 0;
          simError();
          theRsw = 0.95 * theRcut;
        }
        else{
          theRsw = globals->getRsw();
        }
        
        info[i].setDefaultRcut(theRcut, theRsw);
	
      }
    }
  }
#ifdef IS_MPI
  strcpy(checkPointMsg, "post processing checks out");
  MPIcheckPoint();
#endif // is_mpi

}
  
void SimSetup::initSystemCoords(void){
  int i;

  char* inName;

  (info[0].getConfiguration())->createArrays(info[0].n_atoms);

  for (i = 0; i < info[0].n_atoms; i++)
    info[0].atoms[i]->setCoords();

  if (globals->haveInitialConfig()){
    InitializeFromFile* fileInit;
#ifdef IS_MPI // is_mpi
    if (worldRank == 0){
#endif //is_mpi
      inName = globals->getInitialConfig();
      fileInit = new InitializeFromFile(inName);
#ifdef IS_MPI
    }
    else
      fileInit = new InitializeFromFile(NULL);
#endif
    fileInit->readInit(info); // default velocities on

    delete fileInit;
  }
  else{
    
    // no init from md file 
    
    sprintf(painCave.errMsg,
            "Cannot intialize a simulation without an initial configuration file.\n");
    painCave.isFatal = 1;;
    simError();
    
  }

#ifdef IS_MPI
  strcpy(checkPointMsg, "Successfully read in the initial configuration");
  MPIcheckPoint();
#endif // is_mpi
}


void SimSetup::makeOutNames(void){
  int k;
  string prefix;

  for (k = 0; k < nInfo; k++){
#ifdef IS_MPI
    if (worldRank == 0){
#endif // is_mpi
      
      if(globals->haveFinalConfig()) 
	prefix = getPrefix(globals->getFinalConfig());   
      else
	prefix = getPrefix(inFileName);

      info[k].finalName = prefix + ".eor";	
      info[k].sampleName = prefix + ".dump";
      info[k].statusName = prefix + ".stat";

#ifdef IS_MPI

    }
#endif // is_mpi
  }
}


void SimSetup::sysObjectsCreation(void){
  int i, k;

  // create the forceField

  createFF();

  // extract componentList

  compList();

  // calc the number of atoms, bond, bends, and torsions

  calcSysValues();

#ifdef IS_MPI
  // divide the molecules among the processors

  mpiMolDivide();
#endif //is_mpi

  // create the atom and SRI arrays. Also initialize Molecule Stamp ID's

  makeSysArrays();

  // make and initialize the molecules (all but atomic coordinates)

  makeMolecules();

  for (k = 0; k < nInfo; k++){
    info[k].identArray = new int[info[k].n_atoms];
    for (i = 0; i < info[k].n_atoms; i++){
      info[k].identArray[i] = info[k].atoms[i]->getIdent();
    }
  }
}


void SimSetup::createFF(void){
  switch (ffCase){
    case FF_DUFF:
        the_ff = new DUFF();
      break;

    case FF_LJ:
      the_ff = new LJFF();
      break;

    case FF_EAM:
      if (has_forcefield_variant) 
        the_ff = new EAM_FF(forcefield_variant);
      else
        the_ff = new EAM_FF();
      break;

    case FF_H2O:
      the_ff = new WATER();
      break;

    default:
      sprintf(painCave.errMsg,
              "SimSetup Error. Unrecognized force field in case statement.\n");
      painCave.isFatal = 1;
      simError();
  }


#ifdef IS_MPI
  strcpy(checkPointMsg, "ForceField creation successful");
  MPIcheckPoint();
#endif // is_mpi
}


void SimSetup::compList(void){
  int i;
  char* id;
  LinkedMolStamp* headStamp = new LinkedMolStamp();
  LinkedMolStamp* currentStamp = NULL;
  comp_stamps = new MoleculeStamp * [n_components];
  bool haveCutoffGroups;

  haveCutoffGroups = false;
  
  // make an array of molecule stamps that match the components used.
  // also extract the used stamps out into a separate linked list

  for (i = 0; i < nInfo; i++){
    info[i].nComponents = n_components;
    info[i].componentsNmol = components_nmol;
    info[i].compStamps = comp_stamps;
    info[i].headStamp = headStamp;
  }


  for (i = 0; i < n_components; i++){
    id = the_components[i]->getType();
    comp_stamps[i] = NULL;

    // check to make sure the component isn't already in the list

    comp_stamps[i] = headStamp->match(id);
    if (comp_stamps[i] == NULL){
      // extract the component from the list;

      currentStamp = stamps->extractMolStamp(id);
      if (currentStamp == NULL){
        sprintf(painCave.errMsg,
                "SimSetup error: Component \"%s\" was not found in the "
                "list of declared molecules\n",
                id);
        painCave.isFatal = 1;
        simError();
      }

      headStamp->add(currentStamp);
      comp_stamps[i] = headStamp->match(id);
    }

    if(comp_stamps[i]->getNCutoffGroups() > 0)
      haveCutoffGroups = true;    
  }
    
  for (i = 0; i < nInfo; i++)
    info[i].haveCutoffGroups = haveCutoffGroups;

#ifdef IS_MPI
  strcpy(checkPointMsg, "Component stamps successfully extracted\n");
  MPIcheckPoint();
#endif // is_mpi
}

void SimSetup::calcSysValues(void){
  int i, j;
  int ncutgroups, atomsingroups, ngroupsinstamp;

  int* molMembershipArray;
  CutoffGroupStamp* cg;

  tot_atoms = 0;
  tot_bonds = 0;
  tot_bends = 0;
  tot_torsions = 0;
  tot_rigid = 0;
  tot_groups = 0;
  for (i = 0; i < n_components; i++){
    tot_atoms += components_nmol[i] * comp_stamps[i]->getNAtoms();
    tot_bonds += components_nmol[i] * comp_stamps[i]->getNBonds();
    tot_bends += components_nmol[i] * comp_stamps[i]->getNBends();
    tot_torsions += components_nmol[i] * comp_stamps[i]->getNTorsions();
    tot_rigid += components_nmol[i] * comp_stamps[i]->getNRigidBodies();

    ncutgroups = comp_stamps[i]->getNCutoffGroups();
    atomsingroups = 0;
    for (j=0; j < ncutgroups; j++) {
      cg = comp_stamps[i]->getCutoffGroup(j);
      atomsingroups += cg->getNMembers();
    }
    ngroupsinstamp = comp_stamps[i]->getNAtoms() - atomsingroups + ncutgroups;
    tot_groups += components_nmol[i] * ngroupsinstamp;    
  }
  
  tot_SRI = tot_bonds + tot_bends + tot_torsions;
  molMembershipArray = new int[tot_atoms];

  for (i = 0; i < nInfo; i++){
    info[i].n_atoms = tot_atoms;
    info[i].n_bonds = tot_bonds;
    info[i].n_bends = tot_bends;
    info[i].n_torsions = tot_torsions;
    info[i].n_SRI = tot_SRI;
    info[i].n_mol = tot_nmol;
    info[i].ngroup = tot_groups;
    info[i].molMembershipArray = molMembershipArray;
  }
}

#ifdef IS_MPI

void SimSetup::mpiMolDivide(void){
  int i, j, k;
  int localMol, allMol;
  int local_atoms, local_bonds, local_bends, local_torsions, local_SRI;
  int local_rigid, local_groups;
  vector<int> globalMolIndex;
  int ncutgroups, atomsingroups, ngroupsinstamp;
  CutoffGroupStamp* cg;

  mpiSim = new mpiSimulation(info);

  mpiSim->divideLabor();
  globalAtomIndex = mpiSim->getGlobalAtomIndex();
  globalGroupIndex = mpiSim->getGlobalGroupIndex();
  //globalMolIndex = mpiSim->getGlobalMolIndex();

  // set up the local variables 

  mol2proc = mpiSim->getMolToProcMap();
  molCompType = mpiSim->getMolComponentType();

  allMol = 0;
  localMol = 0;
  local_atoms = 0;
  local_bonds = 0;
  local_bends = 0;
  local_torsions = 0;
  local_rigid = 0;
  local_groups = 0;
  globalAtomCounter = 0;

  for (i = 0; i < n_components; i++){
    for (j = 0; j < components_nmol[i]; j++){
      if (mol2proc[allMol] == worldRank){
        local_atoms += comp_stamps[i]->getNAtoms();
        local_bonds += comp_stamps[i]->getNBonds();
        local_bends += comp_stamps[i]->getNBends();
        local_torsions += comp_stamps[i]->getNTorsions();
        local_rigid += comp_stamps[i]->getNRigidBodies();

        ncutgroups = comp_stamps[i]->getNCutoffGroups();
        atomsingroups = 0;
        for (k=0; k < ncutgroups; k++) {
          cg = comp_stamps[i]->getCutoffGroup(k);
          atomsingroups += cg->getNMembers();
        }
        ngroupsinstamp = comp_stamps[i]->getNAtoms() - atomsingroups + 
          ncutgroups;
        local_groups += ngroupsinstamp;    

        localMol++;
      }      
      for (k = 0; k < comp_stamps[i]->getNAtoms(); k++){
        info[0].molMembershipArray[globalAtomCounter] = allMol;
        globalAtomCounter++;
      }

      allMol++;
    }
  }
  local_SRI = local_bonds + local_bends + local_torsions;

  info[0].n_atoms = mpiSim->getNAtomsLocal();  
  
  if (local_atoms != info[0].n_atoms){
    sprintf(painCave.errMsg,
            "SimSetup error: mpiSim's localAtom (%d) and SimSetup's\n"
            "\tlocalAtom (%d) are not equal.\n",
            info[0].n_atoms, local_atoms);
    painCave.isFatal = 1;
    simError();
  }

  info[0].ngroup = mpiSim->getNGroupsLocal();   
  if (local_groups != info[0].ngroup){
    sprintf(painCave.errMsg,
            "SimSetup error: mpiSim's localGroups (%d) and SimSetup's\n"
            "\tlocalGroups (%d) are not equal.\n",
            info[0].ngroup, local_groups);
    painCave.isFatal = 1;
    simError();
  }
  
  info[0].n_bonds = local_bonds;
  info[0].n_bends = local_bends;
  info[0].n_torsions = local_torsions;
  info[0].n_SRI = local_SRI;
  info[0].n_mol = localMol;

  strcpy(checkPointMsg, "Passed nlocal consistency check.");
  MPIcheckPoint();
}

#endif // is_mpi


void SimSetup::makeSysArrays(void){
 
#ifndef IS_MPI
  int k, j;
#endif // is_mpi
  int i, l;

  Atom** the_atoms;
  Molecule* the_molecules;

  for (l = 0; l < nInfo; l++){
    // create the atom and short range interaction arrays

    the_atoms = new Atom * [info[l].n_atoms];
    the_molecules = new Molecule[info[l].n_mol];
    int molIndex;

    // initialize the molecule's stampID's

#ifdef IS_MPI


    molIndex = 0;
    for (i = 0; i < mpiSim->getNMolGlobal(); i++){
      if (mol2proc[i] == worldRank){
        the_molecules[molIndex].setStampID(molCompType[i]);
        the_molecules[molIndex].setMyIndex(molIndex);
        the_molecules[molIndex].setGlobalIndex(i);
        molIndex++;
      }
    }

#else // is_mpi

    molIndex = 0;
    globalAtomCounter = 0;
    for (i = 0; i < n_components; i++){
      for (j = 0; j < components_nmol[i]; j++){
        the_molecules[molIndex].setStampID(i);
        the_molecules[molIndex].setMyIndex(molIndex);
        the_molecules[molIndex].setGlobalIndex(molIndex);
        for (k = 0; k < comp_stamps[i]->getNAtoms(); k++){
          info[l].molMembershipArray[globalAtomCounter] = molIndex;
          globalAtomCounter++;
        }
        molIndex++;
      }
    }


#endif // is_mpi

    info[l].globalExcludes = new int;
    info[l].globalExcludes[0] = 0;
    
    // set the arrays into the SimInfo object

    info[l].atoms = the_atoms;
    info[l].molecules = the_molecules;
    info[l].nGlobalExcludes = 0;
    
    the_ff->setSimInfo(info);
  }
}

void SimSetup::makeIntegrator(void){
  int k;

  NVE<Integrator<BaseIntegrator> >* myNVE = NULL;
  NVT<Integrator<BaseIntegrator> >* myNVT = NULL;
  NPTi<NPT<Integrator<BaseIntegrator> > >* myNPTi = NULL;
  NPTf<NPT<Integrator<BaseIntegrator> > >* myNPTf = NULL;
  NPTxyz<NPT<Integrator<BaseIntegrator> > >* myNPTxyz = NULL;
  
  for (k = 0; k < nInfo; k++){
    switch (ensembleCase){
      case NVE_ENS:
        if (globals->haveZconstraints()){
          setupZConstraint(info[k]);
          myNVE = new ZConstraint<NVE<RealIntegrator> >(&(info[k]), the_ff);
        }
        else{
          myNVE = new NVE<RealIntegrator>(&(info[k]), the_ff);
	}
	
	info->the_integrator = myNVE;
        break;

      case NVT_ENS:
        if (globals->haveZconstraints()){
          setupZConstraint(info[k]);
          myNVT = new ZConstraint<NVT<RealIntegrator> >(&(info[k]), the_ff);
        }
        else
          myNVT = new NVT<RealIntegrator>(&(info[k]), the_ff);

        
        if (globals->haveTargetTemp())
          myNVT->setTargetTemp(globals->getTargetTemp());
        else{
          sprintf(painCave.errMsg,
                  "SimSetup error: If you use the NVT\n"
                  "\tensemble, you must set targetTemp.\n");
          painCave.isFatal = 1;
          simError();
        }

        if (globals->haveTauThermostat())
          myNVT->setTauThermostat(globals->getTauThermostat());
        else{
          sprintf(painCave.errMsg,
                  "SimSetup error: If you use the NVT\n"
                  "\tensemble, you must set tauThermostat.\n");
          painCave.isFatal = 1;
          simError();
        }

	info->the_integrator = myNVT;
        break;

      case NPTi_ENS:
        if (globals->haveZconstraints()){
          setupZConstraint(info[k]);
          myNPTi = new ZConstraint<NPTi<NPT <RealIntegrator> > >(&(info[k]), the_ff);
        }
        else
          myNPTi = new NPTi<NPT<RealIntegrator> >(&(info[k]), the_ff);

        if (globals->haveTargetTemp())
          myNPTi->setTargetTemp(globals->getTargetTemp());
        else{
          sprintf(painCave.errMsg,
                  "SimSetup error: If you use a constant pressure\n"
                  "\tensemble, you must set targetTemp.\n");
          painCave.isFatal = 1;
          simError();
        }

        if (globals->haveTargetPressure())
          myNPTi->setTargetPressure(globals->getTargetPressure());
        else{
          sprintf(painCave.errMsg,
                  "SimSetup error: If you use a constant pressure\n"
                  "\tensemble, you must set targetPressure in the meta-data file.\n");
          painCave.isFatal = 1;
          simError();
        }

        if (globals->haveTauThermostat())
          myNPTi->setTauThermostat(globals->getTauThermostat());
        else{
          sprintf(painCave.errMsg,
                  "SimSetup error: If you use an NPT\n"
                  "\tensemble, you must set tauThermostat.\n");
          painCave.isFatal = 1;
          simError();
        }

        if (globals->haveTauBarostat())
          myNPTi->setTauBarostat(globals->getTauBarostat());
        else{
          sprintf(painCave.errMsg,
                  "SimSetup error: If you use an NPT\n"
                  "\tensemble, you must set tauBarostat.\n");
          painCave.isFatal = 1;
          simError();
        }

	info->the_integrator = myNPTi;
        break;

      case NPTf_ENS:
        if (globals->haveZconstraints()){
          setupZConstraint(info[k]);
          myNPTf = new ZConstraint<NPTf<NPT <RealIntegrator> > >(&(info[k]), the_ff);
        }
        else
          myNPTf = new NPTf<NPT <RealIntegrator> >(&(info[k]), the_ff);

        if (globals->haveTargetTemp())
          myNPTf->setTargetTemp(globals->getTargetTemp());
        else{
          sprintf(painCave.errMsg,
                  "SimSetup error: If you use a constant pressure\n"
                  "\tensemble, you must set targetTemp.\n");
          painCave.isFatal = 1;
          simError();
        }

        if (globals->haveTargetPressure())
          myNPTf->setTargetPressure(globals->getTargetPressure());
        else{
          sprintf(painCave.errMsg,
                  "SimSetup error: If you use a constant pressure\n"
                  "\tensemble, you must set targetPressure in the meta-data file.\n");
          painCave.isFatal = 1;
          simError();
        }    

        if (globals->haveTauThermostat())
          myNPTf->setTauThermostat(globals->getTauThermostat());

        else{
          sprintf(painCave.errMsg,
                  "SimSetup error: If you use an NPT\n"
                  "\tensemble, you must set tauThermostat.\n");
          painCave.isFatal = 1;
          simError();
        }

        if (globals->haveTauBarostat())
          myNPTf->setTauBarostat(globals->getTauBarostat());

        else{
          sprintf(painCave.errMsg,
                  "SimSetup error: If you use an NPT\n"
                  "\tensemble, you must set tauBarostat.\n");
          painCave.isFatal = 1;
          simError();
        }

	info->the_integrator = myNPTf;
        break;

      case NPTxyz_ENS:
        if (globals->haveZconstraints()){
          setupZConstraint(info[k]);
          myNPTxyz = new ZConstraint<NPTxyz<NPT <RealIntegrator> > >(&(info[k]), the_ff);
        }
        else
          myNPTxyz = new NPTxyz<NPT <RealIntegrator> >(&(info[k]), the_ff);

        if (globals->haveTargetTemp())
          myNPTxyz->setTargetTemp(globals->getTargetTemp());
        else{
          sprintf(painCave.errMsg,
                  "SimSetup error: If you use a constant pressure\n"
                  "\tensemble, you must set targetTemp.\n");
          painCave.isFatal = 1;
          simError();
        }

        if (globals->haveTargetPressure())
          myNPTxyz->setTargetPressure(globals->getTargetPressure());
        else{
          sprintf(painCave.errMsg,
                  "SimSetup error: If you use a constant pressure\n"
                  "\tensemble, you must set targetPressure in the meta-data file.\n");
          painCave.isFatal = 1;
          simError();
        }    

        if (globals->haveTauThermostat())
          myNPTxyz->setTauThermostat(globals->getTauThermostat());
        else{
          sprintf(painCave.errMsg,
                  "SimSetup error: If you use an NPT\n"
                  "\tensemble, you must set tauThermostat.\n");
          painCave.isFatal = 1;
          simError();
        }

        if (globals->haveTauBarostat())
          myNPTxyz->setTauBarostat(globals->getTauBarostat());
        else{
          sprintf(painCave.errMsg,
                  "SimSetup error: If you use an NPT\n"
                  "\tensemble, you must set tauBarostat.\n");
          painCave.isFatal = 1;
          simError();
        }

	info->the_integrator = myNPTxyz;
        break;

      default:
        sprintf(painCave.errMsg,
                "SimSetup Error. Unrecognized ensemble in case statement.\n");
        painCave.isFatal = 1;
        simError();
    }
  }
}

void SimSetup::initFortran(void){
  info[0].refreshSim();

  the_ff->initForceField();

#ifdef IS_MPI
  strcpy(checkPointMsg, "Successfully intialized the fortran portion of the force field.");
  MPIcheckPoint();
#endif // is_mpi
}

void SimSetup::setupZConstraint(SimInfo& theInfo){
  int nZConstraints;
  ZconStamp** zconStamp;

  if (globals->haveZconstraintTime()){
    //add sample time of z-constraint  into SimInfo's property list                    
    DoubleGenericData* zconsTimeProp = new DoubleGenericData();
    zconsTimeProp->setID(ZCONSTIME_ID);
    zconsTimeProp->setData(globals->getZconsTime());
    theInfo.addProperty(zconsTimeProp);
  }
  else{
    sprintf(painCave.errMsg,
            "ZConstraint error: If you use a ZConstraint,\n"
            "\tyou must set zconsTime.\n");
    painCave.isFatal = 1;
    simError();
  }

  //push zconsTol into siminfo, if user does not specify
  //value for zconsTol, a default value will be used
  DoubleGenericData* zconsTol = new DoubleGenericData();
  zconsTol->setID(ZCONSTOL_ID);
  if (globals->haveZconsTol()){
    zconsTol->setData(globals->getZconsTol());
  }
  else{
    double defaultZConsTol = 0.01;
    sprintf(painCave.errMsg,
            "ZConstraint Warning: Tolerance for z-constraint method is not specified.\n"
            "\tOOPSE will use a default value of %f.\n"
            "\tTo set the tolerance, use the zconsTol variable.\n",
            defaultZConsTol);
    painCave.isFatal = 0;
    simError();      

    zconsTol->setData(defaultZConsTol);
  }
  theInfo.addProperty(zconsTol);

  //set Force Subtraction Policy
  StringGenericData* zconsForcePolicy = new StringGenericData();
  zconsForcePolicy->setID(ZCONSFORCEPOLICY_ID);

  if (globals->haveZconsForcePolicy()){
    zconsForcePolicy->setData(globals->getZconsForcePolicy());
  }
  else{
    sprintf(painCave.errMsg,
            "ZConstraint Warning: No force subtraction policy was set.\n"
            "\tOOPSE will use PolicyByMass.\n"
            "\tTo set the policy, use the zconsForcePolicy variable.\n");
    painCave.isFatal = 0;
    simError(); 
    zconsForcePolicy->setData("BYMASS");
  }

  theInfo.addProperty(zconsForcePolicy);

  //set zcons gap
  DoubleGenericData* zconsGap = new DoubleGenericData();
  zconsGap->setID(ZCONSGAP_ID);

  if (globals->haveZConsGap()){
    zconsGap->setData(globals->getZconsGap());
    theInfo.addProperty(zconsGap);  
  }

  //set zcons fixtime
  DoubleGenericData* zconsFixtime = new DoubleGenericData();
  zconsFixtime->setID(ZCONSFIXTIME_ID);

  if (globals->haveZConsFixTime()){
    zconsFixtime->setData(globals->getZconsFixtime());
    theInfo.addProperty(zconsFixtime);  
  }

  //set zconsUsingSMD
  IntGenericData* zconsUsingSMD = new IntGenericData();
  zconsUsingSMD->setID(ZCONSUSINGSMD_ID);

  if (globals->haveZConsUsingSMD()){
    zconsUsingSMD->setData(globals->getZconsUsingSMD());
    theInfo.addProperty(zconsUsingSMD);  
  }

  //Determine the name of ouput file and add it into SimInfo's property list 
  //Be careful, do not use inFileName, since it is a pointer which
  //point to a string at master node, and slave nodes do not contain that string

  string zconsOutput(theInfo.finalName);

  zconsOutput = zconsOutput.substr(0, zconsOutput.rfind(".")) + ".fz";

  StringGenericData* zconsFilename = new StringGenericData();
  zconsFilename->setID(ZCONSFILENAME_ID);
  zconsFilename->setData(zconsOutput);

  theInfo.addProperty(zconsFilename);

  //setup index, pos and other parameters of z-constraint molecules
  nZConstraints = globals->getNzConstraints();
  theInfo.nZconstraints = nZConstraints;

  zconStamp = globals->getZconStamp();
  ZConsParaItem tempParaItem;

  ZConsParaData* zconsParaData = new ZConsParaData();
  zconsParaData->setID(ZCONSPARADATA_ID);

  for (int i = 0; i < nZConstraints; i++){
    tempParaItem.havingZPos = zconStamp[i]->haveZpos();
    tempParaItem.zPos = zconStamp[i]->getZpos();
    tempParaItem.zconsIndex = zconStamp[i]->getMolIndex();
    tempParaItem.kRatio = zconStamp[i]->getKratio();
    tempParaItem.havingCantVel = zconStamp[i]->haveCantVel();
    tempParaItem.cantVel = zconStamp[i]->getCantVel();    
    zconsParaData->addItem(tempParaItem);
  }

  //check the uniqueness of index  
  if(!zconsParaData->isIndexUnique()){
    sprintf(painCave.errMsg,
            "ZConstraint Error: molIndex is not unique!\n");
    painCave.isFatal = 1;
    simError(); 
  }

  //sort the parameters by index of molecules
  zconsParaData->sortByIndex();
  
  //push data into siminfo, therefore, we can retrieve later
  theInfo.addProperty(zconsParaData);
}

void SimSetup::makeMinimizer(){

  OOPSEMinimizer* myOOPSEMinimizer;
  MinimizerParameterSet* param;
  char minimizerName[100];
  
  for (int i = 0; i < nInfo; i++){
    
    //prepare parameter set for minimizer
    param = new MinimizerParameterSet();
    param->setDefaultParameter();

    if (globals->haveMinimizer()){
      param->setFTol(globals->getMinFTol());
    }

    if (globals->haveMinGTol()){
      param->setGTol(globals->getMinGTol());
    }

    if (globals->haveMinMaxIter()){
      param->setMaxIteration(globals->getMinMaxIter());
    }

    if (globals->haveMinWriteFrq()){
      param->setMaxIteration(globals->getMinMaxIter());
    }

    if (globals->haveMinWriteFrq()){
      param->setWriteFrq(globals->getMinWriteFrq());
    }
    
    if (globals->haveMinStepSize()){
      param->setStepSize(globals->getMinStepSize());
    }

    if (globals->haveMinLSMaxIter()){
      param->setLineSearchMaxIteration(globals->getMinLSMaxIter());
    }    

    if (globals->haveMinLSTol()){
      param->setLineSearchTol(globals->getMinLSTol());
    }     

    strcpy(minimizerName, globals->getMinimizer());

    if (!strcasecmp(minimizerName, "CG")){
      myOOPSEMinimizer = new PRCGMinimizer(&(info[i]), the_ff, param);
    }
    else if (!strcasecmp(minimizerName, "SD")){
    //myOOPSEMinimizer = MinimizerFactory.creatMinimizer("", &(info[i]), the_ff, param);
      myOOPSEMinimizer = new SDMinimizer(&(info[i]), the_ff, param);
    }
    else{
          sprintf(painCave.errMsg,
                  "SimSetup error: Unrecognized Minimizer, use Conjugate Gradient \n");
          painCave.isFatal = 0;
          simError();

      myOOPSEMinimizer = new PRCGMinimizer(&(info[i]), the_ff, param);          
    }
     info[i].the_integrator = myOOPSEMinimizer;

     //store the minimizer into simInfo
     info[i].the_minimizer = myOOPSEMinimizer;
     info[i].has_minimizer = true;
  }

}
