/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */

#define _LARGEFILE_SOURCE64
#define _FILE_OFFSET_BITS 64

#include <sys/types.h>
#include <sys/stat.h>

#include <iostream>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "primitives/Molecule.hpp"
#include "utils/MemoryUtils.hpp"
#include "utils/StringTokenizer.hpp"
#include "io/RestReader.hpp"
#include "utils/simError.h"

#ifdef IS_MPI
#include <mpi.h>
#define TAKE_THIS_TAG_CHAR 0
#define TAKE_THIS_TAG_INT 1
#define TAKE_THIS_TAG_DOUBLE 2
#endif // is_mpi

namespace oopse {
  
  RestReader::RestReader( SimInfo* info ) : info_(info){
        
    idealName = "idealCrystal.in";
    
    isScanned = false;
    
#ifdef IS_MPI
    if (worldRank == 0) {
#endif
      
      inIdealFile = fopen(idealName, "r");
      if(inIdealFile == NULL){
        sprintf(painCave.errMsg,
                "RestReader: Cannot open file: %s\n", idealName);
        painCave.isFatal = 1;
        simError();
      }
      
      inIdealFileName = idealName;
#ifdef IS_MPI
    }
    strcpy( checkPointMsg, 
            "File \"idealCrystal.in\" opened successfully for reading." );
    MPIcheckPoint();
#endif
    return;  
  }
  
  RestReader :: ~RestReader( ){
#ifdef IS_MPI
    if (worldRank == 0) {
#endif
      int error;
      error = fclose( inIdealFile );
      
      if( error ){
        sprintf( painCave.errMsg,
                 "Error closing %s\n", inIdealFileName.c_str());
        simError();
      }
      
      MemoryUtils::deletePointers(framePos);
      
#ifdef IS_MPI
    }
    strcpy( checkPointMsg, "Restraint file closed successfully." );
    MPIcheckPoint();
#endif
    
    return;
  }
  
  
  void RestReader :: readIdealCrystal(){
    
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
    
    std::vector<StuntDouble*> integrableObjects;
    
    Molecule* mol;
    StuntDouble* integrableObject;
    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ii;
    
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
    
    if( nTotObjs != info_->getNGlobalIntegrableObjects() ){
      sprintf( painCave.errMsg,
               "RestraintReader error. %s nIntegrable, %d, "
               "does not match the meta-data file's nIntegrable, %d.\n",
               inIdealFileName.c_str(), nTotObjs, 
               info_->getNGlobalIntegrableObjects());
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
     * Note: we assume that there is a one-to-one correspondence between
     * integrable objects and lines in the idealCrystal.in file.  Thermodynamic
     * integration is only supported for simple rigid bodies. 
     */
    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) {
      
      for (integrableObject = mol->beginIntegrableObject(ii); 
           integrableObject != NULL; 
           integrableObject = mol->nextIntegrableObject(ii)) {    
        
        eof_test = fgets(read_buffer, sizeof(read_buffer), inIdealFile);
        if(eof_test == NULL){
          sprintf(painCave.errMsg,
                  "RestReader Error: error in reading file %s\n"
                  "natoms  = %d; index = %d\n"
                  "error reading the line from the file.\n",
                  inIdealFileName.c_str(), nTotObjs, 
                  integrableObject->getGlobalIndex() );
          painCave.isFatal = 1;
          simError();
        }
        
        parseErr = parseIdealLine( read_buffer, integrableObject);
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
    
    int masterNode = 0;
    int myStatus; // 1 = wakeup & success; 0 = error; -1 = AllDone
    
    MPI_Status istatus;
    int nCurObj;
    int nitems;
    int haveError;

    nTotObjs = info_->getNGlobalIntegrableObjects();
    haveError = 0;
    
    if (worldRank == masterNode) {
      eof_test = fgets(read_buffer, sizeof(read_buffer), inIdealFile);
      if( eof_test == NULL ){
        sprintf( painCave.errMsg,
                 "Error reading 1st line of %s \n ",inIdealFileName.c_str());
        painCave.isFatal = 1;
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
                 inIdealFileName.c_str(), nTotObjs, 
                 info_->getNGlobalIntegrableObjects());
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
      
      for (i=0 ; i < info_->getNGlobalMolecules(); i++) {
        int which_node = info_->getMolToProc(i);
        
        if(which_node == masterNode){
          //molecules belong to master node
          
          mol = info_->getMoleculeByGlobalIndex(i);
          
          if(mol == NULL) {
	    sprintf(painCave.errMsg, 
                   "RestReader Error: Molecule not found on node %d!\n",
                   worldRank);
            painCave.isFatal = 1;
            simError();
          }
          
          for (integrableObject = mol->beginIntegrableObject(ii); 
               integrableObject != NULL; 
               integrableObject = mol->nextIntegrableObject(ii)){
            
            eof_test = fgets(read_buffer, sizeof(read_buffer), inIdealFile);
            
            if(eof_test == NULL){
              sprintf(painCave.errMsg,
                      "RestReader Error: error in reading file %s\n"
                      "natoms  = %d; index = %d\n"
                      "error reading the line from the file.\n",
                      inIdealFileName.c_str(), nTotObjs, i );
              painCave.isFatal = 1;
              simError();
            }
            
            parseIdealLine(read_buffer, integrableObjects[j]);
          }
        } else {
          //molecule belongs to slave nodes
          
          MPI_Recv(&nCurObj, 1, MPI_INT, which_node,
                   TAKE_THIS_TAG_INT, MPI_COMM_WORLD, &istatus);
          
          for(j=0; j < nCurObj; j++){
            
            eof_test = fgets(read_buffer, sizeof(read_buffer), inIdealFile);
            if(eof_test == NULL){
              sprintf(painCave.errMsg,
                      "RestReader Error: error in reading file %s\n"
                      "natoms  = %d; index = %d\n"
                      "error reading the line from the file.\n",
                      inIdealFileName.c_str(), nTotObjs, i );
              painCave.isFatal = 1;
              simError();
            }
            
	    MPI_Send(read_buffer, BUFFERSIZE, MPI_CHAR, which_node,
                     TAKE_THIS_TAG_CHAR, MPI_COMM_WORLD);
          }
        }
      }
    } else {
      //actions taken at slave nodes
      for (i=0 ; i < info_->getNGlobalMolecules(); i++) {
        int which_node = info_->getMolToProc(i);
        
        if(which_node == worldRank){
          //molecule with global index i belongs to this processor
          
          mol = info_->getMoleculeByGlobalIndex(i);
          
          if(mol == NULL) {
            sprintf(painCave.errMsg, 
                    "RestReader Error: molecule not found on node %d\n", 
                    worldRank);
            painCave.isFatal = 1;
            simError();
          }
          
          nCurObj = mol->getNIntegrableObjects();
          
          MPI_Send(&nCurObj, 1, MPI_INT, masterNode,
                   TAKE_THIS_TAG_INT, MPI_COMM_WORLD);
          
          for (integrableObject = mol->beginIntegrableObject(ii); 
               integrableObject != NULL; 
               integrableObject = mol->nextIntegrableObject(ii)){
            
            MPI_Recv(read_buffer, BUFFERSIZE, MPI_CHAR, masterNode,
                     TAKE_THIS_TAG_CHAR, MPI_COMM_WORLD, &istatus);
            
            parseErr = parseIdealLine(read_buffer, integrableObject);
            
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
  
  char* RestReader::parseIdealLine(char* readLine, StuntDouble* sd){
    
    char *foo; // the pointer to the current string token
    
    double pos[3];        // position place holders
    double q[4];          // the quaternions
    double RfromQ[3][3];  // the rotation matrix 
    double normalize;     // to normalize the reference unit vector
    double uX, uY, uZ;    // reference unit vector place holders
    double uselessToken;
    StringTokenizer tokenizer(readLine);
    int nTokens;
    
    nTokens = tokenizer.countTokens();
    
    if (nTokens < 14) {
      sprintf(painCave.errMsg,
              "RestReader Error: Not enough Tokens.\n");
      painCave.isFatal = 1;
      simError();
    }
    
    std::string name = tokenizer.nextToken();
    
    if (name != sd->getType()) {
      
      sprintf(painCave.errMsg,
              "RestReader Error: Atom type [%s] in %s does not "
              "match Atom Type [%s] in .md file.\n",
              name.c_str(), inIdealFileName.c_str(), 
              sd->getType().c_str());
      painCave.isFatal = 1;
      simError();        
    }
    
    pos[0] = tokenizer.nextTokenAsDouble();
    pos[1] = tokenizer.nextTokenAsDouble();
    pos[2] = tokenizer.nextTokenAsDouble();
    
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
    uselessToken = tokenizer.nextTokenAsDouble();
    uselessToken = tokenizer.nextTokenAsDouble();
    uselessToken = tokenizer.nextTokenAsDouble();
    
    if (sd->isDirectional()) {
      
      q[0] = tokenizer.nextTokenAsDouble();
      q[1] = tokenizer.nextTokenAsDouble();
      q[2] = tokenizer.nextTokenAsDouble();
      q[3] = tokenizer.nextTokenAsDouble();
      
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
  
  void RestReader::readZangle(){
    
    int i;
    unsigned int j;
    int isPresent;
    
    Molecule* mol;
    StuntDouble* integrableObject;
    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ii;
    
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
    
    std::vector<StuntDouble*> vecParticles;
    std::vector<double> tempZangs;
      
    inAngFileName = info_->getRestFileName();
    
    inAngFileName += "0";
    
    // open the omega value file for reading
#ifdef IS_MPI
    if (worldRank == 0) {
#endif
      isPresent = 1;
      inAngFile = fopen(inAngFileName.c_str(), "r");
      if(!inAngFile){
        sprintf(painCave.errMsg,
                "Restraints Warning: %s file is not present\n"
                "\tAll omega values will be initialized to zero. If the\n"
                "\tsimulation is starting from the idealCrystal.in reference\n"
                "\tconfiguration, this is the desired action. If this is not\n"
                "\tthe case, the energy calculations will be incorrect.\n",
                inAngFileName.c_str());
        painCave.severity = OOPSE_WARNING;
        painCave.isFatal = 0;
        simError();   
        isPresent = 0;
      }
      
#ifdef IS_MPI
      // let the other nodes know the status of the file search
      MPI_Bcast(&isPresent, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif // is_mpi
      
      if (!isPresent) {
        zeroZangle();
        return;
      }
      
#ifdef IS_MPI
    }
    
    // listen to node 0 to see if we should exit this function
    if (worldRank != 0) {
      MPI_Bcast(&isPresent, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (!isPresent) {
        zeroZangle();
        return;
      }
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
    
    nTotObjs = info_->getNGlobalIntegrableObjects();
    
    if( nTotObjs != tempZangs.size() ){
      sprintf( painCave.errMsg,
               "RestraintReader zAngle reading error. %s nIntegrable, %d, "
               "does not match the meta-data file's nIntegrable, %d.\n",
               inAngFileName.c_str(), tempZangs.size(), nTotObjs );
      painCave.isFatal = 1;
      simError();
    }
    
    // load the zAngles into the integrable objects
    i = 0;
    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) {
      
      for (integrableObject = mol->beginIntegrableObject(ii); 
           integrableObject != NULL; 
           integrableObject = mol->nextIntegrableObject(ii)) {    
        
        integrableObject->setZangle(tempZangs[i]);
        i++;
      }
    }
    
    // MPI Section of code..........
#else //IS_MPI
    
    // first thing first, suspend fatalities.
    painCave.isEventLoop = 1;

    int masterNode = 0;
    int myStatus; // 1 = wakeup & success; 0 = error; -1 = AllDone
    int haveError;
    int index;    

    int nCurObj;
    double angleTranfer;
    
    nTotObjs = info_->getNGlobalIntegrableObjects();
    haveError = 0;

    if (worldRank == masterNode) {
      
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
      
      // At this point, node 0 has a tempZangs vector completed, and 
      // everyone else has nada
      index = 0;
      
      for (i=0 ; i < info_->getNGlobalMolecules(); i++) {
	// Get the Node number which has this atom
	which_node = info_->getMolToProc(i);
	
	if (worldRank == masterNode) {
	  mol = info_->getMoleculeByGlobalIndex(i);
	  
	  if(mol == NULL) {
	    strcpy(painCave.errMsg, "Molecule not found on node 0!");
	    haveError = 1;
	    simError();
	  }
	  
	  for (integrableObject = mol->beginIntegrableObject(ii); 
	       integrableObject != NULL; 
	       integrableObject = mol->nextIntegrableObject(ii)){
	    
	    integrableObject->setZangle(tempZangs[index]);
	    index++;
	  }	
	  
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
      }
    } else {
      // I am SLAVE TO THE MASTER
      for (i=0 ; i < info_->getNGlobalMolecules(); i++) {
        int which_node = info_->getMolToProc(i);

	if (which_node == worldRank) {
	  
	  // BUT I OWN THIS MOLECULE!!!
	  
	  mol = info_->getMoleculeByGlobalIndex(i);

          if(mol == NULL) {
            sprintf(painCave.errMsg, 
                    "RestReader Error: molecule not found on node %d\n", 
                    worldRank);
            painCave.isFatal = 1;
            simError();
          }

          nCurObj = mol->getNIntegrableObjects();
	
	  MPI_Send(&nCurObj, 1, MPI_INT, 0,
		   TAKE_THIS_TAG_INT, MPI_COMM_WORLD);
	  
          for (integrableObject = mol->beginIntegrableObject(ii); 
               integrableObject != NULL; 
               integrableObject = mol->nextIntegrableObject(ii)){
	    
	    MPI_Recv(&angleTransfer, 1, MPI_DOUBLE, 0,
		     TAKE_THIS_TAG_DOUBLE, MPI_COMM_WORLD, &istatus);

	    integrableObject->setZangle(angleTransfer);
	  }	
	}
      }
    } 
#endif
  }
  
  void RestReader :: zeroZangle(){
    
    int i;
    unsigned int j;
    int nTotObjs; // the number of atoms
    
    Molecule* mol;
    StuntDouble* integrableObject;
    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ii;
    
    std::vector<StuntDouble*> vecParticles;
    
#ifndef IS_MPI
    // set all zAngles to 0.0
    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) 
      
      for (integrableObject = mol->beginIntegrableObject(ii); 
           integrableObject != NULL; 
           integrableObject = mol->nextIntegrableObject(ii))    
        integrableObject->setZangle( 0.0 );
    
    
    // MPI Section of code..........
#else //IS_MPI
    
    // first thing first, suspend fatalities.
    painCave.isEventLoop = 1;
    
    int masterNode = 0;
    int myStatus; // 1 = wakeup & success; 0 = error; -1 = AllDone
    int haveError;
    int which_node;
    
    MPI_Status istatus;
    
    int nCurObj;
    double angleTranfer;
    
    nTotObjs = info_->getNGlobalIntegrableObjects();
    haveError = 0;
    if (worldRank == masterNode) {

      for (i=0 ; i < info_->getNGlobalMolecules(); i++) {
	// Get the Node number which has this atom
	which_node = info_->getMolToProc(i);
	
	// let's let node 0 pass out constant values to all the processors
	if (worldRank == masterNode) {
	  mol = info_->getMoleculeByGlobalIndex(i);
	  
	  if(mol == NULL) {
	    strcpy(painCave.errMsg, "Molecule not found on node 0!");
	    haveError = 1;
	    simError();
	  }
	  
	  for (integrableObject = mol->beginIntegrableObject(ii); 
	       integrableObject != NULL; 
	       integrableObject = mol->nextIntegrableObject(ii)){
	    
	    integrableObject->setZangle( 0.0 );
	    
	  }
	  
	} else {
	  // I am MASTER OF THE UNIVERSE, but I don't own this molecule
	  
	  MPI_Recv(&nCurObj, 1, MPI_INT, which_node,
		   TAKE_THIS_TAG_INT, MPI_COMM_WORLD, &istatus);
	  
	  for(j=0; j < nCurObj; j++){	 	 
	    angleTransfer = 0.0;
	    MPI_Send(&angleTransfer, 1, MPI_DOUBLE, which_node, 
		     TAKE_THIS_TAG_DOUBLE, MPI_COMM_WORLD);      	  
	    
	  }
	}
      }
    } else {
      // I am SLAVE TO THE MASTER
      for (i=0 ; i < info_->getNGlobalMolecules(); i++) {
	int which_node = info_->getMolToProc(i);
	
	if (which_node == worldRank) {
	  
	  // BUT I OWN THIS MOLECULE!!!
	  mol = info_->getMoleculeByGlobalIndex(i);
	  
	  if(mol == NULL) {
	    sprintf(painCave.errMsg, 
		    "RestReader Error: molecule not found on node %d\n", 
		    worldRank);
	    painCave.isFatal = 1;
	    simError();
	  }
	  
	  nCurObj = mol->getNIntegrableObjects();
	  
	  MPI_Send(&nCurObj, 1, MPI_INT, 0,
		   TAKE_THIS_TAG_INT, MPI_COMM_WORLD);
	  
	  for (integrableObject = mol->beginIntegrableObject(ii); 
               integrableObject != NULL; 
               integrableObject = mol->nextIntegrableObject(ii)){
	    
            MPI_Recv(&angleTransfer, 1, MPI_DOUBLE, 0,
                     TAKE_THIS_TAG_DOUBLE, MPI_COMM_WORLD, &istatus);
            vecParticles[j]->setZangle(angleTransfer);
          }	
        }
      }
    }
#endif
  }
  
} // end namespace oopse
