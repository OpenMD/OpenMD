#define _LARGEFILE_SOURCE64
#define _FILE_OFFSET_BITS 64

#include <sys/types.h>
#include <sys/stat.h>

#include <iostream>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "ReadWrite.hpp"
#include "simError.h"

#ifdef IS_MPI
#include <mpi.h>
#include "mpiSimulation.hpp"
#define TAKE_THIS_TAG_CHAR 0
#define TAKE_THIS_TAG_INT 1
#endif // is_mpi


DumpReader :: DumpReader(const char *in_name ){

  isScanned = false;

#ifdef IS_MPI
  if (worldRank == 0) {
#endif

  inFile = fopen(in_name, "r");
  if(inFile == NULL){
    sprintf(painCave.errMsg,
	    "Cannot open file: %s\n", in_name);
    painCave.isFatal = 1;
    simError();
  }
  
  inFileName = in_name;
#ifdef IS_MPI
  }
  strcpy( checkPointMsg, "Dump file opened for reading successfully." );
  MPIcheckPoint();
#endif
  return;  
}

DumpReader :: ~DumpReader( ){
#ifdef IS_MPI
  if (worldRank == 0) {
#endif
  vector<fpos_t*>::iterator i;

  int error;
  error = fclose( inFile );
  if( error ){
    sprintf( painCave.errMsg,
	     "Error closing %s\n", inFileName.c_str());
    simError();
  }

  for(i = framePos.begin(); i != framePos.end(); ++i)
    delete *i;
  framePos.clear();
  
#ifdef IS_MPI
  }
  strcpy( checkPointMsg, "Dump file closed successfully." );
  MPIcheckPoint();
#endif

  return;
}

int DumpReader::getNframes( void ){

  if( !isScanned )
    scanFile();
  return framePos.size();
}

void DumpReader::scanFile( void ){

  int i, j;
  int lineNum = 0;
  char readBuffer[2000];
  fpos_t *currPos;

#ifdef IS_MPI
  if( worldRank == 0 ){
#endif // is_mpi
    
    rewind( inFile );
    
    currPos = new fpos_t;
    fgetpos( inFile, currPos );
    fgets( readBuffer, sizeof( readBuffer ), inFile );
    lineNum++;
    if( feof( inFile ) ){
      sprintf( painCave.errMsg,
	       "File \"%s\" ended unexpectedly at line %d\n", 
	       inFileName.c_str(),
	       lineNum );
      painCave.isFatal = 1;
      simError();
    }

    while( !feof( inFile ) ){
      
      framePos.push_back(currPos);

      i = atoi(readBuffer);
      
      fgets( readBuffer, sizeof( readBuffer ), inFile );
      lineNum++;    
      if( feof( inFile ) ){
	sprintf( painCave.errMsg,
		 "File \"%s\" ended unexpectedly at line %d\n",
		 inFileName.c_str(),
		 lineNum );
	painCave.isFatal = 1;
	simError();
      }
            
      for(j=0; j<i; j++){
	
	fgets( readBuffer, sizeof( readBuffer ), inFile );
	lineNum++;    
	if( feof( inFile ) ){
	  sprintf( painCave.errMsg,
		   "File \"%s\" ended unexpectedly at line %d,"
		   " with atom %d\n", 
		   inFileName.c_str(),
		   lineNum, 
		   j );
	  painCave.isFatal = 1;
	  simError();
	}
	
      }
      
      currPos = new fpos_t;
      fgetpos( inFile, currPos );
      fgets( readBuffer, sizeof( readBuffer ), inFile );
      lineNum++;
    }
    
    delete currPos;
    rewind( inFile );
     
    isScanned = true;

#ifdef IS_MPI
  }
  strcpy( checkPointMsg, "Successfully scanned DumpFile\n" );
  MPIcheckPoint();
#endif // is_mpi
}

void DumpReader :: readFrame( SimInfo* the_simnfo, int whichFrame){

  simnfo = the_simnfo;

  this->readSet( whichFrame );
}



void DumpReader :: readSet( int whichFrame ){

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

  fsetpos(inFile, framePos[whichFrame]);
  eof_test = fgets(read_buffer, sizeof(read_buffer), inFile);
  if( eof_test == NULL ){
    sprintf( painCave.errMsg,
	     "DumpReader error: error reading 1st line of \"%s\"\n",
	     inFileName.c_str() );
    painCave.isFatal = 1;
    simError();
  }

  nTotObjs = atoi( read_buffer );

  if( nTotObjs != simnfo->getTotIntegrableObjects() ){
    sprintf( painCave.errMsg,
	     "DumpReader error. %s nIntegrable, %d, "
	     "does not match the meta-data file's nIntegrable, %d.\n",
	     inFileName.c_str(), nTotObjs, simnfo->getTotIntegrableObjects());
    painCave.isFatal = 1;
    simError();
  }

  //read the box mat from the comment line

  eof_test = fgets(read_buffer, sizeof(read_buffer), inFile);
  if(eof_test == NULL){
    sprintf( painCave.errMsg,
	     "error in reading commment in %s\n", inFileName.c_str());
    painCave.isFatal = 1;
    simError();
  }

  parseErr = parseCommentLine( read_buffer, simnfo);
  if( parseErr != NULL ){
    strcpy( painCave.errMsg, parseErr );
    painCave.isFatal = 1;
    simError();
  }

  //parse dump lines

  for( i=0; i < simnfo->n_mol; i++){

    integrableObjects = (simnfo->molecules[i]).getIntegrableObjects();

    for(j = 0; j < integrableObjects.size(); j++){

      eof_test = fgets(read_buffer, sizeof(read_buffer), inFile);
      if(eof_test == NULL){
        sprintf(painCave.errMsg,
  	      "error in reading file %s\n"
  	      "natoms  = %d; index = %d\n"
  	      "error reading the line from the file.\n",
  	      inFileName.c_str(), nTotObjs, i );
        painCave.isFatal = 1;
        simError();
      }
      
      parseErr = parseDumpLine( read_buffer, integrableObjects[j]);
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
     fsetpos(inFile,  framePos[whichFrame]);

    eof_test = fgets(read_buffer, sizeof(read_buffer), inFile);
    if( eof_test == NULL ){
      sprintf( painCave.errMsg,
	       "Error reading 1st line of %s \n ",inFileName.c_str());
      haveError = 1;
      simError();
    }

    nitems = atoi( read_buffer );

    // Check to see that the number of integrable objects in the
    // intial configuration file is the same as derived from the
    // meta-data file.

    if( nTotObjs != nitems){
      sprintf( painCave.errMsg,
	       "DumpReader Error. %s nIntegrable, %d, "
	       "does not match the meta-data file's nIntegrable, %d.\n",
	       inFileName.c_str(), nTotObjs, simnfo->getTotIntegrableObjects());
      haveError= 1;
      simError();
    }

    //read the boxMat from the comment line

    eof_test = fgets(read_buffer, sizeof(read_buffer), inFile);
    if(eof_test == NULL){
      sprintf( painCave.errMsg,
	       "error in reading commment in %s\n", inFileName.c_str());
      haveError = 1;
      simError();
    }

    //Every single processor will parse the comment line by itself
    //By using this way, we might lose some efficiency, but if we want to add
    //more parameters into comment line, we only need to modify function
    //parseCommentLine

    MPI_Bcast(read_buffer, BUFFERSIZE, MPI_CHAR, 0, MPI_COMM_WORLD);

    parseErr = parseCommentLine( read_buffer, simnfo);

    if( parseErr != NULL ){
      strcpy( painCave.errMsg, parseErr );
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
        
          eof_test = fgets(read_buffer, sizeof(read_buffer), inFile);
          if(eof_test == NULL){
        	sprintf(painCave.errMsg,
		    "error in reading file %s\n"
		    "natoms  = %d; index = %d\n"
		    "error reading the line from the file.\n",
		    inFileName.c_str(), nTotObjs, i );
	        haveError= 1;
	        simError();
          }
          
          if(haveError) nodeZeroError();

          parseDumpLine(read_buffer, integrableObjects[j]);
           
       }


      }
      else{
      //molecule belongs to slave nodes

        MPI_Recv(&nCurObj, 1, MPI_INT, which_node,
	       TAKE_THIS_TAG_INT, MPI_COMM_WORLD, &istatus);
      
       for(j=0; j < nCurObj; j++){
        
          eof_test = fgets(read_buffer, sizeof(read_buffer), inFile);
          if(eof_test == NULL){
        	sprintf(painCave.errMsg,
		    "error in reading file %s\n"
		    "natoms  = %d; index = %d\n"
		    "error reading the line from the file.\n",
		    inFileName.c_str(), nTotObjs, i );
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
    MPI_Bcast(read_buffer, BUFFERSIZE, MPI_CHAR, 0, MPI_COMM_WORLD);

    parseErr = parseCommentLine( read_buffer, simnfo);

    if( parseErr != NULL ){
      strcpy( painCave.errMsg, parseErr );
      haveError = 1;
      simError();
    }
  
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

          parseErr = parseDumpLine(read_buffer, integrableObjects[j]);

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

char* DumpReader::parseDumpLine(char* readLine, StuntDouble* sd){

  char *foo; // the pointer to the current string token

  double pos[3]; // position place holders
  double vel[3]; // velocity placeholders
  double q[4]; // the quaternions
  double ji[3]; // angular velocity placeholders;
  double qSqr, qLength; // needed to normalize the quaternion vector.


  // set the string tokenizer

  foo = strtok(readLine, " ,;\t");

  // check the atom name to the current atom

  if( strcmp( foo, sd->getType() ) ){
    sprintf( painCave.errMsg,
	     "DumpReader error.  Does not"
	     " match the meta-data atom %s.\n",
	     sd->getType() );
    return strdup( painCave.errMsg );
  }

  // get the positions

  foo = strtok(NULL, " ,;\t");
  if(foo == NULL){
    sprintf( painCave.errMsg,
	     "error in reading postition x from %s\n",
	     inFileName.c_str());
    return strdup( painCave.errMsg );
  }
  pos[0] = atof( foo );

  foo = strtok(NULL, " ,;\t");
  if(foo == NULL){
    sprintf( painCave.errMsg,
	     "error in reading postition y from %s\n",
	     inFileName.c_str());
    return strdup( painCave.errMsg );
  }
  pos[1] = atof( foo );

  foo = strtok(NULL, " ,;\t");
  if(foo == NULL){
    sprintf( painCave.errMsg,
	     "error in reading postition z from %s\n",
	     inFileName.c_str());
    return strdup( painCave.errMsg );
  }
  pos[2] = atof( foo );


  // get the velocities

  foo = strtok(NULL, " ,;\t");
  if(foo == NULL){
    sprintf( painCave.errMsg,
	     "error in reading velocity x from %s\n",
	     inFileName.c_str() );
    return strdup( painCave.errMsg );
  }
  vel[0] = atof( foo );

  foo = strtok(NULL, " ,;\t");
  if(foo == NULL){
    sprintf( painCave.errMsg,
	     "error in reading velocity x from %s\n",
	     inFileName.c_str() );
    return strdup( painCave.errMsg );
  }
  vel[1] = atof( foo );

  foo = strtok(NULL, " ,;\t");
  if(foo == NULL){
    sprintf( painCave.errMsg,
	     "error in reading velocity x from %s\n",
	     inFileName.c_str() );
    return strdup( painCave.errMsg );
  }
  vel[2] = atof( foo );


  // add the positions and velocities to the atom

  sd->setPos( pos );
  sd->setVel( vel );

  if (!sd->isDirectional())
    return NULL;

  // get the quaternions

  if( sd->isDirectional() ){

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      sprintf( painCave.errMsg,
	             "error in reading velocity x from %s\n",
	              inFileName.c_str() );
      return strdup( painCave.errMsg );
    }
    q[0] = atof( foo );

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      sprintf( painCave.errMsg,
	             "error in reading velocity x from %s\n",
	              inFileName.c_str() );
      return strdup( painCave.errMsg );
    }
    q[1] = atof( foo );

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      sprintf( painCave.errMsg,
	             "error in reading velocity x from %s\n",
	              inFileName.c_str() );
      return strdup( painCave.errMsg );
    }
    q[2] = atof( foo );

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      sprintf( painCave.errMsg,
	             "error in reading velocity x from %s\n",
	              inFileName.c_str() );
      return strdup( painCave.errMsg );
    }
    q[3] = atof( foo );

    // get the angular velocities

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      sprintf( painCave.errMsg,
	             "error in reading velocity x from %s\n",
	              inFileName.c_str() );
      return strdup( painCave.errMsg );
    }
    ji[0] = atof( foo );

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      sprintf( painCave.errMsg,
	             "error in reading velocity x from %s\n",
	              inFileName.c_str() );
      return strdup( painCave.errMsg );
    }
    ji[1] = atof(foo );

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      sprintf( painCave.errMsg,
	             "error in reading velocity x from %s\n",
	              inFileName.c_str() );
      return strdup( painCave.errMsg );
    }
    ji[2] = atof( foo );


    // check that the quaternion vector is normalized

    qSqr = (q[0] * q[0]) + (q[1] * q[1]) + (q[2] * q[2]) + (q[3] * q[3]);

    if (fabs(qSqr) < 1e-6) {
      sprintf(painCave.errMsg,
          "initial quaternion error (q0^2 + q1^2 + q2^2 + q3^2 ~ 0).\n");
       return strdup(painCave.errMsg);
    }

    qLength = sqrt( qSqr );
    q[0] = q[0] / qLength;
    q[1] = q[1] / qLength;
    q[2] = q[2] / qLength;
    q[3] = q[3] / qLength;

    // add quaternion and angular velocities

    sd->setQ( q );
    sd->setJ( ji );
  }



  return NULL;
}


char* DumpReader::parseCommentLine(char* readLine, SimInfo* entry_plug){

  double currTime;
  double boxMat[9];
  double theBoxMat3[3][3];
  double chi;
  double integralOfChidt;
  double eta[9];

  char *foo; // the pointer to the current string token

  // set the string tokenizer

  foo = strtok(readLine, " ,;\t");
  // set the timeToken.

  if(foo == NULL){
    sprintf( painCave.errMsg,
	     "error in reading Time from %s\n",
	     inFileName.c_str() );
    return strdup( painCave.errMsg );
  }

  currTime = atof( foo );
  entry_plug->setTime( currTime );

  //get H-Matrix

  for(int i = 0 ; i < 9; i++){
    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      sprintf( painCave.errMsg,
               "error in reading H[%d] from %s\n", i, inFileName.c_str() );
      return strdup( painCave.errMsg );
    }
    boxMat[i] = atof( foo );
  }

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++) theBoxMat3[i][j] = boxMat[3*j+i];

  //set H-Matrix
  entry_plug->setBoxM( theBoxMat3 );

  //get chi and integralOfChidt, they should appear by pair

  if( entry_plug->useInitXSstate ){
    foo = strtok(NULL, " ,;\t\n");
    if(foo != NULL){
      chi = atof(foo);
      
      foo = strtok(NULL, " ,;\t\n");
      if(foo == NULL){
	sprintf( painCave.errMsg,
		 "chi and integralOfChidt should appear by pair in %s\n", inFileName.c_str() );
	return strdup( painCave.errMsg );
      }
      integralOfChidt = atof( foo );
      
      //push chi and integralOfChidt into SimInfo::properties which can be
      //retrieved by integrator later
      DoubleData* chiValue = new DoubleData();
      chiValue->setID(CHIVALUE_ID);
      chiValue->setData(chi);
      entry_plug->addProperty(chiValue);
      
      DoubleData* integralOfChidtValue = new DoubleData();
      integralOfChidtValue->setID(INTEGRALOFCHIDT_ID);
      integralOfChidtValue->setData(integralOfChidt);
      entry_plug->addProperty(integralOfChidtValue);
      
    }
    else
      return NULL;
    
    //get eta
    foo = strtok(NULL, " ,;\t\n");
    if(foo != NULL ){
  
      for(int i = 0 ; i < 9; i++){
	
	if(foo == NULL){
	  sprintf( painCave.errMsg,
		   "error in reading eta[%d] from %s\n", i, inFileName.c_str() );
	  return strdup( painCave.errMsg );
	}
	eta[i] = atof( foo );
	foo = strtok(NULL, " ,;\t\n");
      }
    }
    else
      return NULL;
    
    //push eta into SimInfo::properties which can be
    //retrieved by integrator later
    //entry_plug->setBoxM( theBoxMat3 );
    DoubleArrayData* etaValue = new DoubleArrayData();
    etaValue->setID(ETAVALUE_ID);
    etaValue->setData(eta, 9);
    entry_plug->addProperty(etaValue);
  }

  return NULL;
}

#ifdef IS_MPI
void DumpReader::nodeZeroError( void ){
  int j, myStatus;

  myStatus = 0;
  for (j = 0; j < mpiSim->getNProcessors(); j++) {
    MPI_Send( &myStatus, 1, MPI_INT, j,
	      TAKE_THIS_TAG_INT, MPI_COMM_WORLD);
  }


  MPI_Finalize();
  exit (0);

}

void DumpReader::anonymousNodeDie( void ){

  MPI_Finalize();
  exit (0);
}
#endif
