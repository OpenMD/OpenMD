#define _LARGEFILE_SOURCE64
#define _FILE_OFFSET_BITS 64

#include <sys/types.h>
#include <sys/stat.h>

#include <iostream>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "io/ReadWrite.hpp"
#include "utils/simError.h"

#ifdef IS_MPI
#include <mpi.h>
#include "brains/mpiSimulation.hpp"
#define TAKE_THIS_TAG_CHAR 0
#define TAKE_THIS_TAG_INT 1
#define TAKE_THIS_TAG_DOUBLE 2
#endif // is_mpi


RestraintReader :: RestraintReader( SimInfo* the_simnfo ){
  
  simnfo = the_simnfo;

  idealName = "idealCrystal.in";
  
  isScanned = false;

#ifdef IS_MPI
  if (worldRank == 0) {
#endif

  inIdealFile = fopen(idealName, "r");
  if(inIdealFile == NULL){
    sprintf(painCave.errMsg,
	    "Cannot open file: %s\n", idealName);
    painCave.isFatal = 1;
    simError();
  }
  
  inIdealFileName = idealName;
#ifdef IS_MPI
  }
  strcpy( checkPointMsg, "Restraint file opened for reading successfully." );
  MPIcheckPoint();
#endif
  return;  
}

RestraintReader :: ~RestraintReader( ){
#ifdef IS_MPI
  if (worldRank == 0) {
#endif
  vector<fpos_t*>::iterator i;

  int error;
  error = fclose( inIdealFile );
  if( error ){
    sprintf( painCave.errMsg,
	     "Error closing %s\n", inIdealFileName.c_str());
    simError();
  }

  for(i = framePos.begin(); i != framePos.end(); ++i)
    delete *i;
  framePos.clear();
  
#ifdef IS_MPI
  }
  strcpy( checkPointMsg, "Restraint file closed successfully." );
  MPIcheckPoint();
#endif

  return;
}


void RestraintReader :: readIdealCrystal(){

  int i;
  unsigned int j;

#ifdef IS_MPI
  int done, which_node, which_atom; // loop counter
#endif //is_mpi

  const int BUFFERSIZE = 2000; // size of the read buffer
  int nTotObjs; // the number of atoms
  char read_buffer[BUFFERSIZE]; //the line buffer for reading

  char *eof_test; // ptr to see when we reach the end of the file
  char *parseErr;

  vector<StuntDouble*> integrableObjects;


#ifndef IS_MPI

  eof_test = fgets(read_buffer, sizeof(read_buffer), inIdealFile);
  if( eof_test == NULL ){
    sprintf( painCave.errMsg,
	     "RestraintReader error: error reading 1st line of \"%s\"\n",
	     inIdealFileName.c_str() );
    painCave.isFatal = 1;
    simError();
  }

  nTotObjs = atoi( read_buffer );

  if( nTotObjs != simnfo->getTotIntegrableObjects() ){
    sprintf( painCave.errMsg,
	     "RestraintReader error. %s nIntegrable, %d, "
	     "does not match the meta-data file's nIntegrable, %d.\n",
	     inIdealFileName.c_str(), nTotObjs, simnfo->getTotIntegrableObjects());
    painCave.isFatal = 1;
    simError();
  }

  // skip over the comment line
  eof_test = fgets(read_buffer, sizeof(read_buffer), inIdealFile);
  if(eof_test == NULL){
    sprintf( painCave.errMsg,
	     "error in reading commment in %s\n", inIdealFileName.c_str());
    painCave.isFatal = 1;
    simError();
  }

  // parse the ideal crystal lines
  /* 
     Note: we assume that there is a one-to-one correspondence between
     integrable objects and lines in the idealCrystal.in file.  Thermodynamic
     integration is only supported for simple rigid bodies. 
  */
  for( i=0; i < simnfo->n_mol; i++){

    integrableObjects = (simnfo->molecules[i]).getIntegrableObjects();

    for(j = 0; j < integrableObjects.size(); j++){

      eof_test = fgets(read_buffer, sizeof(read_buffer), inIdealFile);
      if(eof_test == NULL){
        sprintf(painCave.errMsg,
  	      "error in reading file %s\n"
  	      "natoms  = %d; index = %d\n"
  	      "error reading the line from the file.\n",
  	      inIdealFileName.c_str(), nTotObjs, i );
        painCave.isFatal = 1;
        simError();
      }
      
      parseErr = parseIdealLine( read_buffer, integrableObjects[j]);
      if( parseErr != NULL ){
        strcpy( painCave.errMsg, parseErr );
        painCave.isFatal = 1;
        simError();
      }
    }
  }

  // MPI Section of code..........
#else //IS_MPI

  // first thing first, suspend fatalities.
  painCave.isEventLoop = 1;

  int myStatus; // 1 = wakeup & success; 0 = error; -1 = AllDone
  int haveError;

  MPI_Status istatus;
  int *MolToProcMap = mpiSim->getMolToProcMap();
  int localIndex;
  int nCurObj;
  int nitems;

  nTotObjs = simnfo->getTotIntegrableObjects();
  haveError = 0;
  if (worldRank == 0) {
    eof_test = fgets(read_buffer, sizeof(read_buffer), inIdealFile);
    if( eof_test == NULL ){
      sprintf( painCave.errMsg,
	       "Error reading 1st line of %s \n ",inIdealFileName.c_str());
      haveError = 1;
      simError();
    }
    
    nitems = atoi( read_buffer );
    
    // Check to see that the number of integrable objects in the
    // intial configuration file is the same as derived from the
    // meta-data file.
    if( nTotObjs != nitems){
      sprintf( painCave.errMsg,
	       "RestraintReader Error. %s nIntegrable, %d, "
	       "does not match the meta-data file's nIntegrable, %d.\n",
	       inIdealFileName.c_str(), nTotObjs, simnfo->getTotIntegrableObjects());
      haveError= 1;
      simError();
    }
    
    // skip over the comment line
    eof_test = fgets(read_buffer, sizeof(read_buffer), inIdealFile);
    if(eof_test == NULL){
      sprintf( painCave.errMsg,
	       "error in reading commment in %s\n", inIdealFileName.c_str());
      haveError = 1;
      simError();
    }
    
    for (i=0 ; i < mpiSim->getNMolGlobal(); i++) {
      which_node = MolToProcMap[i];
      if(which_node == 0){
	//molecules belong to master node
	
	localIndex = mpiSim->getGlobalToLocalMol(i);
	
	if(localIndex == -1) {
	  strcpy(painCave.errMsg, "Molecule not found on node 0!");
	  haveError = 1;
	  simError();
	}
	
	integrableObjects = (simnfo->molecules[localIndex]).getIntegrableObjects();
	for(j=0; j < integrableObjects.size(); j++){
	  
          eof_test = fgets(read_buffer, sizeof(read_buffer), inIdealFile);
          if(eof_test == NULL){
	    sprintf(painCave.errMsg,
		    "error in reading file %s\n"
		    "natoms  = %d; index = %d\n"
		    "error reading the line from the file.\n",
		    inIdealFileName.c_str(), nTotObjs, i );
	    haveError= 1;
	    simError();
          }
          
          if(haveError) nodeZeroError();
	  
          parseIdealLine(read_buffer, integrableObjects[j]);
	}
      }
      else{
	//molecule belongs to slave nodes
	
        MPI_Recv(&nCurObj, 1, MPI_INT, which_node,
		 TAKE_THIS_TAG_INT, MPI_COMM_WORLD, &istatus);
	
	for(j=0; j < nCurObj; j++){
	  
          eof_test = fgets(read_buffer, sizeof(read_buffer), inIdealFile);
          if(eof_test == NULL){
	    sprintf(painCave.errMsg,
		    "error in reading file %s\n"
		    "natoms  = %d; index = %d\n"
		    "error reading the line from the file.\n",
		    inIdealFileName.c_str(), nTotObjs, i );
	    haveError= 1;
	    simError();
          }
          
          if(haveError) nodeZeroError();
	  
	  MPI_Send(read_buffer, BUFFERSIZE, MPI_CHAR, which_node,
		   TAKE_THIS_TAG_CHAR, MPI_COMM_WORLD);
	}
      }
    }
  }
  else{
    //actions taken at slave nodes
    for (i=0 ; i < mpiSim->getNMolGlobal(); i++) {
      which_node = MolToProcMap[i];
      
      if(which_node == worldRank){
	//molecule with global index i belongs to this processor
	
        localIndex = mpiSim->getGlobalToLocalMol(i);
	
        if(localIndex == -1) {
          sprintf(painCave.errMsg, "Molecule not found on node %d\n", worldRank);
          haveError = 1;
          simError();
        }
	
        integrableObjects = (simnfo->molecules[localIndex]).getIntegrableObjects();        
	
        nCurObj = integrableObjects.size();
        
        MPI_Send(&nCurObj, 1, MPI_INT, 0,
		 TAKE_THIS_TAG_INT, MPI_COMM_WORLD);
	
        for(j = 0; j < integrableObjects.size(); j++){
	  
          MPI_Recv(read_buffer, BUFFERSIZE, MPI_CHAR, 0,
		   TAKE_THIS_TAG_CHAR, MPI_COMM_WORLD, &istatus);
	  
          parseErr = parseIdealLine(read_buffer, integrableObjects[j]);
	  
          if( parseErr != NULL ){
	    strcpy( painCave.errMsg, parseErr );
	    simError();
          }
        }
      }
    }
  }
#endif
}

char* RestraintReader::parseIdealLine(char* readLine, StuntDouble* sd){

  char *foo; // the pointer to the current string token

  double pos[3];        // position place holders
  double q[4];          // the quaternions
  double RfromQ[3][3];  // the rotation matrix 
  double normalize;     // to normalize the reference unit vector
  double uX, uY, uZ;    // reference unit vector place holders

  // set the string tokenizer

  foo = strtok(readLine, " ,;\t");

  // check the atom name to the current atom

  if( strcmp( foo, sd->getType() ) ){
    sprintf( painCave.errMsg,
	     "RestraintReader error.  Does not"
	     " match the meta-data atom %s.\n",
	     sd->getType() );
    return strdup( painCave.errMsg );
  }

  // get the positions

  foo = strtok(NULL, " ,;\t");
  if(foo == NULL){
    sprintf( painCave.errMsg,
	     "error in reading position x from %s\n",
	     inIdealFileName.c_str());
    return strdup( painCave.errMsg );
  }
  pos[0] = atof( foo );

  foo = strtok(NULL, " ,;\t");
  if(foo == NULL){
    sprintf( painCave.errMsg,
	     "error in reading position y from %s\n",
	     inIdealFileName.c_str());
    return strdup( painCave.errMsg );
  }
  pos[1] = atof( foo );

  foo = strtok(NULL, " ,;\t");
  if(foo == NULL){
    sprintf( painCave.errMsg,
	     "error in reading position z from %s\n",
	     inIdealFileName.c_str());
    return strdup( painCave.errMsg );
  }
  pos[2] = atof( foo );

  // store the positions in the stuntdouble as generic data doubles
  DoubleGenericData* refPosX = new DoubleGenericData();
  refPosX->setID("refPosX");
  refPosX->setData(pos[0]);
  sd->addProperty(refPosX);
  
  DoubleGenericData* refPosY = new DoubleGenericData();
  refPosY->setID("refPosY");
  refPosY->setData(pos[1]);
  sd->addProperty(refPosY);
  
  DoubleGenericData* refPosZ = new DoubleGenericData();
  refPosZ->setID("refPosZ");
  refPosZ->setData(pos[2]);
  sd->addProperty(refPosZ);

  // we don't need the velocities
  foo = strtok(NULL, " ,;\t");
  foo = strtok(NULL, " ,;\t");
  foo = strtok(NULL, " ,;\t");

  if (!sd->isDirectional())
    return NULL;

  // get the quaternions
  if( sd->isDirectional() ){

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      sprintf( painCave.errMsg,
	             "error in reading quaternion 0 from %s\n",
	              inIdealFileName.c_str() );
      return strdup( painCave.errMsg );
    }
    q[0] = atof( foo );

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      sprintf( painCave.errMsg,
	             "error in reading quaternion 1 from %s\n",
	              inIdealFileName.c_str() );
      return strdup( painCave.errMsg );
    }
    q[1] = atof( foo );

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      sprintf( painCave.errMsg,
	             "error in reading quaternion 2 from %s\n",
	              inIdealFileName.c_str() );
      return strdup( painCave.errMsg );
    }
    q[2] = atof( foo );

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      sprintf( painCave.errMsg,
	             "error in reading quaternion 3 from %s\n",
	              inIdealFileName.c_str() );
      return strdup( painCave.errMsg );
    }
    q[3] = atof( foo );

    // now build the rotation matrix and find the unit vectors
    RfromQ[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
    RfromQ[0][1] = 2*(q[1]*q[2] + q[0]*q[3]);
    RfromQ[0][2] = 2*(q[1]*q[3] - q[0]*q[2]);
    RfromQ[1][0] = 2*(q[1]*q[2] - q[0]*q[3]);
    RfromQ[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
    RfromQ[1][2] = 2*(q[2]*q[3] + q[0]*q[1]);
    RfromQ[2][0] = 2*(q[1]*q[3] + q[0]*q[2]);
    RfromQ[2][1] = 2*(q[2]*q[3] - q[0]*q[1]);
    RfromQ[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
    
    normalize = sqrt(RfromQ[2][0]*RfromQ[2][0] + RfromQ[2][1]*RfromQ[2][1] 
		     + RfromQ[2][2]*RfromQ[2][2]);
    uX = RfromQ[2][0]/normalize;
    uY = RfromQ[2][1]/normalize;
    uZ = RfromQ[2][2]/normalize;

    // store reference unit vectors as generic data in the stuntdouble
    DoubleGenericData* refVectorX = new DoubleGenericData();
    refVectorX->setID("refVectorX");
    refVectorX->setData(uX);
    sd->addProperty(refVectorX);

    DoubleGenericData* refVectorY = new DoubleGenericData();
    refVectorY->setID("refVectorY");
    refVectorY->setData(uY);
    sd->addProperty(refVectorY);

    DoubleGenericData* refVectorZ = new DoubleGenericData();
    refVectorZ->setID("refVectorZ");
    refVectorZ->setData(uZ);
    sd->addProperty(refVectorZ);
  }

  // we don't need the angular velocities, so let's exit the line
  return NULL;
}

#ifdef IS_MPI
void RestraintReader::nodeZeroError( void ){
  int j, myStatus;

  myStatus = 0;
  for (j = 0; j < mpiSim->getNProcessors(); j++) {
    MPI_Send( &myStatus, 1, MPI_INT, j,
	      TAKE_THIS_TAG_INT, MPI_COMM_WORLD);
  }


  MPI_Finalize();
  exit (0);

}

void RestraintReader::anonymousNodeDie( void ){

  MPI_Finalize();
  exit (0);
}
#endif

void RestraintReader::readZangle( const char *in_name ){

  int i;
  unsigned int j;
  int isPresent;

#ifdef IS_MPI
  int done, which_node, which_atom; // loop counter
  int nProc;
  MPI_Status istatus;
#endif //is_mpi

  const int BUFFERSIZE = 2000; // size of the read buffer
  int nTotObjs; // the number of atoms
  char read_buffer[BUFFERSIZE]; //the line buffer for reading

  char *eof_test; // ptr to see when we reach the end of the file
  char *parseErr;

  vector<StuntDouble*> vecParticles;
  vector<double> tempZangs;

  // open the omega value file for reading
#ifdef IS_MPI
  if (worldRank == 0) {
#endif
    isPresent = 1;
    inAngFile = fopen(in_name, "r");
    if(!inAngFile){
      sprintf(painCave.errMsg,
	      "Restraints Warning: %s file is not present\n"
	      "\tAll omega values will be initialized to zero. If the\n"
	      "\tsimulation is starting from the idealCrystal.in reference\n"
	      "\tconfiguration, this is the desired action. If this is not\n"
	      "\tthe case, the energy calculations will be incorrect.\n",
	      in_name);
      painCave.severity = OOPSE_WARNING;
      painCave.isFatal = 0;
      simError();   
      isPresent = 0;
    }

#ifdef IS_MPI
    // let the other nodes know the status of the file search
    MPI_Bcast(&isPresent, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif // is_mpi

    if (!isPresent)
      return;
    
    inAngFileName = in_name;
#ifdef IS_MPI
  }

  // listen to node 0 to see if we should exit this function
  if (worldRank != 0) {
    MPI_Bcast(&isPresent, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!isPresent)
      return;
  }
  
  strcpy( checkPointMsg, "zAngle file opened successfully for reading." );
  MPIcheckPoint();
#endif

#ifndef IS_MPI

  eof_test = fgets(read_buffer, sizeof(read_buffer), inAngFile);
  if( eof_test == NULL ){
    sprintf( painCave.errMsg,
	     "RestraintReader error: error reading 1st line of \"%s\"\n",
	     inAngFileName.c_str() );
    painCave.isFatal = 1;
    simError();
  }

  eof_test = fgets(read_buffer, sizeof(read_buffer), inAngFile);
  while ( eof_test != NULL ) {
    // check for and ignore blank lines
    if ( read_buffer != NULL )
      tempZangs.push_back( atof(read_buffer) );
    eof_test = fgets(read_buffer, sizeof(read_buffer), inAngFile);
  }

  nTotObjs = simnfo->getTotIntegrableObjects();

  if( nTotObjs != tempZangs.size() ){
    sprintf( painCave.errMsg,
	     "RestraintReader zAngle reading error. %s nIntegrable, %d, "
	     "does not match the meta-data file's nIntegrable, %d.\n",
	     inAngFileName.c_str(), tempZangs.size(), nTotObjs );
    painCave.isFatal = 1;
    simError();
  }

  // load the zAngles into the integrable objects
  vecParticles = simnfo->integrableObjects;

  for ( i=0; i<vecParticles.size(); i++ ) {
    vecParticles[i]->setZangle(tempZangs[i]);
  }

  // MPI Section of code..........
#else //IS_MPI

  // first thing first, suspend fatalities.
  painCave.isEventLoop = 1;

  int myStatus; // 1 = wakeup & success; 0 = error; -1 = AllDone
  int haveError, index;

  int *MolToProcMap = mpiSim->getMolToProcMap();
  int localIndex;
  int nCurObj;
  double angleTranfer;

  nTotObjs = simnfo->getTotIntegrableObjects();
  haveError = 0;
  if (worldRank == 0) {

    eof_test = fgets(read_buffer, sizeof(read_buffer), inAngFile);
    if( eof_test == NULL ){
      sprintf( painCave.errMsg,
	       "Error reading 1st line of %s \n ",inAngFileName.c_str());
      haveError = 1;
      simError();
    }

    // let node 0 load the temporary angle vector
    eof_test = fgets(read_buffer, sizeof(read_buffer), inAngFile);
    while ( eof_test != NULL ) {
      // check for and ignore blank lines
      if ( read_buffer != NULL )
	tempZangs.push_back( atof(read_buffer) );
      eof_test = fgets(read_buffer, sizeof(read_buffer), inAngFile);
    }
 
    // Check to see that the number of integrable objects in the
    // intial configuration file is the same as derived from the
    // meta-data file.
    if( nTotObjs != tempZangs.size() ){
      sprintf( painCave.errMsg,
	       "RestraintReader zAngle reading Error. %s nIntegrable, %d, "
	       "does not match the meta-data file's nIntegrable, %d.\n",
	       inAngFileName.c_str(), tempZangs.size(), nTotObjs);
      haveError= 1;
      simError();
    }

  }
  // At this point, node 0 has a tempZangs vector completed, and 
  // everyone else has nada
  index = 0;

  for (i=0 ; i < mpiSim->getNMolGlobal(); i++) {
    // Get the Node number which has this atom
    which_node = MolToProcMap[i];
    
    if (worldRank == 0) {
      if (which_node == 0) {
	localIndex = mpiSim->getGlobalToLocalMol(i);
	
	if(localIndex == -1) {
	  strcpy(painCave.errMsg, "Molecule not found on node 0!");
	  haveError = 1;
	  simError();
	}

        vecParticles = (simnfo->molecules[localIndex]).getIntegrableObjects();	
        for(j = 0; j < vecParticles.size(); j++){	  
	  vecParticles[j]->setZangle(tempZangs[index]);
	  index++;
	}	
	
// 	// restraints is limited to a single zAngle per molecule
// 	vecParticles = (simnfo->molecules[localIndex]).getIntegrableObjects();
// 	for(j=0; j < vecParticles.size(); j++)
// 	  vecParticles[j]->setZangle(tempZangs[i]);
	
      } else {
	// I am MASTER OF THE UNIVERSE, but I don't own this molecule

	MPI_Recv(&nCurObj, 1, MPI_INT, which_node,
		 TAKE_THIS_TAG_INT, MPI_COMM_WORLD, &istatus);

	for(j=0; j < nCurObj; j++){	 	 
	  angleTransfer = tempZangs[index];
	  MPI_Send(&angleTransfer, 1, MPI_DOUBLE, which_node, 
		   TAKE_THIS_TAG_DOUBLE, MPI_COMM_WORLD);      	  
	  index++;
	}

      }
      
    } else {
      // I am SLAVE TO THE MASTER

      if (which_node == worldRank) {

	// BUT I OWN THIS MOLECULE!!!

	localIndex = mpiSim->getGlobalToLocalMol(i);
        vecParticles = (simnfo->molecules[localIndex]).getIntegrableObjects();	
        nCurObj = vecParticles.size();
        
        MPI_Send(&nCurObj, 1, MPI_INT, 0,
                        TAKE_THIS_TAG_INT, MPI_COMM_WORLD);

        for(j = 0; j < vecParticles.size(); j++){
	  
	  MPI_Recv(&angleTransfer, 1, MPI_DOUBLE, 0,
		   TAKE_THIS_TAG_DOUBLE, MPI_COMM_WORLD, &istatus);
	  vecParticles[j]->setZangle(angleTransfer);
	}	
      }
    }
  }
#endif
}

void RestraintReader :: zeroZangle(){

  int i;
  unsigned int j;
  int nTotObjs; // the number of atoms

  vector<StuntDouble*> vecParticles;

#ifndef IS_MPI
  // set all zAngles to 0.0
  vecParticles = simnfo->integrableObjects;

  for ( i=0; i<vecParticles.size(); i++ ) {
    vecParticles[i]->setZangle( 0.0 );
  }

  // MPI Section of code..........
#else //IS_MPI

  // first thing first, suspend fatalities.
  painCave.isEventLoop = 1;

  int myStatus; // 1 = wakeup & success; 0 = error; -1 = AllDone
  int haveError, index;
  int which_node;

  MPI_Status istatus;
  int *MolToProcMap = mpiSim->getMolToProcMap();
  int localIndex;
  int nCurObj;
  double angleTranfer;

  nTotObjs = simnfo->getTotIntegrableObjects();
  haveError = 0;

  for (i=0 ; i < mpiSim->getNMolGlobal(); i++) {
    // Get the Node number which has this atom
    which_node = MolToProcMap[i];
    
    // let's let node 0 pass out constant values to all the processors
    if (worldRank == 0) {
      if (which_node == 0) {
	localIndex = mpiSim->getGlobalToLocalMol(i);
	
	if(localIndex == -1) {
	  strcpy(painCave.errMsg, "Molecule not found on node 0!");
	  haveError = 1;
	  simError();
	}
	
	vecParticles = (simnfo->molecules[localIndex]).getIntegrableObjects();	
	for(j = 0; j < vecParticles.size(); j++){	  
	  vecParticles[j]->setZangle( 0.0 );
	}	
	
      } else {
	// I am MASTER OF THE UNIVERSE, but I don't own this molecule
	
	MPI_Recv(&nCurObj, 1, MPI_INT, which_node,
		 TAKE_THIS_TAG_INT, MPI_COMM_WORLD, &istatus);
	
	for(j=0; j < nCurObj; j++){	 	 
	  angleTransfer = 0.0;
	  MPI_Send(&angleTransfer, 1, MPI_DOUBLE, which_node, 
		   TAKE_THIS_TAG_DOUBLE, MPI_COMM_WORLD);      	  
	  index++;
	}
      }
    } else {
      // I am SLAVE TO THE MASTER
      
      if (which_node == worldRank) {
	
	// BUT I OWN THIS MOLECULE!!!
	
	localIndex = mpiSim->getGlobalToLocalMol(i);
	vecParticles = (simnfo->molecules[localIndex]).getIntegrableObjects();	
	nCurObj = vecParticles.size();
	
	MPI_Send(&nCurObj, 1, MPI_INT, 0,
		 TAKE_THIS_TAG_INT, MPI_COMM_WORLD);
	
	for(j = 0; j < vecParticles.size(); j++){
	  
	  MPI_Recv(&angleTransfer, 1, MPI_DOUBLE, 0,
		   TAKE_THIS_TAG_DOUBLE, MPI_COMM_WORLD, &istatus);
	  vecParticles[j]->setZangle(angleTransfer);
	}	
      }
    }
  }
#endif
}
