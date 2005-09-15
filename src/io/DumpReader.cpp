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

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"
#include "utils/MemoryUtils.hpp"
#include "utils/StringTokenizer.hpp"

#ifdef IS_MPI

#include <mpi.h>
#define TAKE_THIS_TAG_CHAR 0
#define TAKE_THIS_TAG_INT 1

#endif // is_mpi


namespace oopse {
  
  DumpReader::DumpReader(SimInfo* info, const std::string& filename)
    : info_(info), filename_(filename), isScanned_(false), nframes_(0) {
    
#ifdef IS_MPI
    
      if (worldRank == 0) {
#endif
      
	inFile_ = fopen(filename_.c_str(), "r");
      
	if (inFile_ == NULL) {
	  sprintf(painCave.errMsg, "DumpReader: Cannot open file: %s\n", filename_.c_str());
	  painCave.isFatal = 1;
	  simError();
	}
      
#ifdef IS_MPI
      
      }
    
      strcpy(checkPointMsg, "Dump file opened for reading successfully.");
      MPIcheckPoint();
    
#endif
    
      return;
    }
  
  DumpReader::~DumpReader() {
    
#ifdef IS_MPI
    
    if (worldRank == 0) {
#endif
      
      int error;
      error = fclose(inFile_);
      
      if (error) {
        sprintf(painCave.errMsg, "DumpReader Error: Error closing %s\n", filename_.c_str());
        painCave.isFatal = 1;            
        simError();
      }
      
      MemoryUtils::deletePointers(framePos_);
      
#ifdef IS_MPI
      
    }
    
    strcpy(checkPointMsg, "Dump file closed successfully.");
    MPIcheckPoint();
    
#endif
    
    return;
  }
  
  int DumpReader::getNFrames(void) {
    
    if (!isScanned_)
      scanFile();
    
    return nframes_;
  }
  
  void DumpReader::scanFile(void) {
    int i, j;
    int lineNum = 0;
    char readBuffer[maxBufferSize];
    fpos_t * currPos;
    
#ifdef IS_MPI
    
    if (worldRank == 0) {
#endif // is_mpi
      
      rewind(inFile_);
      
      currPos = new fpos_t;
      fgetpos(inFile_, currPos);
      fgets(readBuffer, sizeof(readBuffer), inFile_);
      lineNum++;
      
      if (feof(inFile_)) {
        sprintf(painCave.errMsg,
                "DumpReader Error: File \"%s\" ended unexpectedly at line %d\n",
                filename_.c_str(),
                lineNum);
        painCave.isFatal = 1;
        simError();
      }
      
      while (!feof(inFile_)) {
        framePos_.push_back(currPos);
        
        i = atoi(readBuffer);
        
        fgets(readBuffer, sizeof(readBuffer), inFile_);
        lineNum++;
        
        if (feof(inFile_)) {
          sprintf(painCave.errMsg,
                  "DumpReader Error: File \"%s\" ended unexpectedly at line %d\n",
                  filename_.c_str(),
                  lineNum);
          painCave.isFatal = 1;
          simError();
        }
        
        for(j = 0; j < i; j++) {
          fgets(readBuffer, sizeof(readBuffer), inFile_);
          lineNum++;
          
          if (feof(inFile_)) {
            sprintf(painCave.errMsg,
                    "DumpReader Error: File \"%s\" ended unexpectedly at line %d,"
                    " with atom %d\n", filename_.c_str(),
                    lineNum,
                    j);
            
            painCave.isFatal = 1;
            simError();
          }
        }
        
        currPos = new fpos_t;
        fgetpos(inFile_, currPos);
        fgets(readBuffer, sizeof(readBuffer), inFile_);
        lineNum++;
      }
      
      delete currPos;
      rewind(inFile_);
      
      nframes_ = framePos_.size();
#ifdef IS_MPI
    }
    
    MPI_Bcast(&nframes_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    strcpy(checkPointMsg, "Successfully scanned DumpFile\n");
    MPIcheckPoint();
    
#endif // is_mpi
    
    isScanned_ = true;
  }
  
  void DumpReader::readFrame(int whichFrame) {
    if (!isScanned_)
      scanFile();
        
    int storageLayout = info_->getSnapshotManager()->getStorageLayout();
    
    if (storageLayout & DataStorage::dslPosition) {
      needPos_ = true;
    } else {
      needPos_ = false;
    }
    
    if (storageLayout & DataStorage::dslVelocity) {
      needVel_ = true;
    } else {
      needVel_ = false;
    }
    
    if (storageLayout & DataStorage::dslAmat || storageLayout & DataStorage::dslElectroFrame) {
      needQuaternion_ = true;
    } else {
      needQuaternion_ = false;
    }
    
    if (storageLayout & DataStorage::dslAngularMomentum) {
      needAngMom_ = true;
    } else {
      needAngMom_ = false;    
    }
    
    readSet(whichFrame);
  }
  
  void DumpReader::readSet(int whichFrame) {
    int i;
    int nTotObjs;                  // the number of atoms
    char read_buffer[maxBufferSize];  //the line buffer for reading
    char * eof_test;               // ptr to see when we reach the end of the file
    
    Molecule* mol;
    StuntDouble* integrableObject;
    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ii;
    
#ifndef IS_MPI
    
    fsetpos(inFile_, framePos_[whichFrame]);
    eof_test = fgets(read_buffer, sizeof(read_buffer), inFile_);
    
    if (eof_test == NULL) {
      sprintf(painCave.errMsg,
              "DumpReader error: error reading 1st line of \"%s\"\n",
              filename_.c_str());
      painCave.isFatal = 1;
      simError();
    }
    
    nTotObjs = atoi(read_buffer);
    
    if (nTotObjs != info_->getNGlobalIntegrableObjects()) {
      sprintf(painCave.errMsg,
              "DumpReader error. %s nIntegrable, %d, "
              "does not match the meta-data file's nIntegrable, %d.\n",
              filename_.c_str(),
              nTotObjs,
              info_->getNGlobalIntegrableObjects());
      
      painCave.isFatal = 1;
      simError();
    }
    
    //read the box mat from the comment line
    
    eof_test = fgets(read_buffer, sizeof(read_buffer), inFile_);
    
    if (eof_test == NULL) {
      sprintf(painCave.errMsg, "DumpReader Error: error in reading commment in %s\n",
              filename_.c_str());
      painCave.isFatal = 1;
      simError();
    }
    
    parseCommentLine(read_buffer, info_->getSnapshotManager()->getCurrentSnapshot());
    
    //parse dump lines
    
    i = 0;
    for (mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)) {
      
      for (integrableObject = mol->beginIntegrableObject(ii); integrableObject != NULL; 
           integrableObject = mol->nextIntegrableObject(ii)) {           
        
        eof_test = fgets(read_buffer, sizeof(read_buffer), inFile_);
        
        if (eof_test == NULL) {
          sprintf(painCave.errMsg,
                  "DumpReader Error: error in reading file %s\n"
                  "natoms  = %d; index = %d\n"
                  "error reading the line from the file.\n",
                  filename_.c_str(),
                  nTotObjs,
                  i);
          
          painCave.isFatal = 1;
          simError();
        }
        
        parseDumpLine(read_buffer, integrableObject);
        i++;
      }
    }
    
    // MPI Section of code..........
    
#else //IS_MPI
    
    // first thing first, suspend fatalities.
    int masterNode = 0;
    int nCurObj;
    painCave.isEventLoop = 1;
    
    int myStatus; // 1 = wakeup & success; 0 = error; -1 = AllDone
    int haveError;
    
    MPI_Status istatus;
    int nitems;
    
    nTotObjs = info_->getNGlobalIntegrableObjects();
    haveError = 0;
    
    if (worldRank == masterNode) {
      fsetpos(inFile_, framePos_[whichFrame]);
      
      eof_test = fgets(read_buffer, sizeof(read_buffer), inFile_);
      
      if (eof_test == NULL) {
        sprintf(painCave.errMsg, "DumpReader Error: Error reading 1st line of %s \n ",
                filename_.c_str());
        painCave.isFatal = 1;
        simError();
      }
      
      nitems = atoi(read_buffer);
      
      // Check to see that the number of integrable objects in the
      // intial configuration file is the same as derived from the
      // meta-data file.
      
      if (nTotObjs != nitems) {
        sprintf(painCave.errMsg,
                "DumpReader Error. %s nIntegrable, %d, "
                "does not match the meta-data file's nIntegrable, %d.\n",
                filename_.c_str(),
                nTotObjs,
                info_->getNGlobalIntegrableObjects());
        
        painCave.isFatal = 1;
        simError();
      }
      
      //read the boxMat from the comment line
      
      eof_test = fgets(read_buffer, sizeof(read_buffer), inFile_);
      
      if (eof_test == NULL) {
        sprintf(painCave.errMsg, "DumpReader Error: error in reading commment in %s\n",
                filename_.c_str());
        painCave.isFatal = 1;
        simError();
      }
      
      //Every single processor will parse the comment line by itself
      //By using this way, we might lose some efficiency, but if we want to add
      //more parameters into comment line, we only need to modify function
      //parseCommentLine
      
      MPI_Bcast(read_buffer, maxBufferSize, MPI_CHAR, masterNode, MPI_COMM_WORLD);
      parseCommentLine(read_buffer, info_->getSnapshotManager()->getCurrentSnapshot());
      
      for(i = 0; i < info_->getNGlobalMolecules(); i++) {
        int which_node = info_->getMolToProc(i);
        
        if (which_node == masterNode) {
          //molecules belong to master node
          
          mol = info_->getMoleculeByGlobalIndex(i);
          
          if (mol == NULL) {
            sprintf(painCave.errMsg, "DumpReader Error: Molecule not found on node %d!", worldRank);
            painCave.isFatal = 1;
            simError();
          }
          
          for (integrableObject = mol->beginIntegrableObject(ii); integrableObject != NULL; 
               integrableObject = mol->nextIntegrableObject(ii)){
            
            eof_test = fgets(read_buffer, sizeof(read_buffer), inFile_);
            
            if (eof_test == NULL) {
              sprintf(painCave.errMsg,
                      "DumpReader Error: error in reading file %s\n"
                      "natoms  = %d; index = %d\n"
                      "error reading the line from the file.\n",
                      filename_.c_str(),
                      nTotObjs,
                      i);
              
              painCave.isFatal = 1;
              simError();
            }
            
            parseDumpLine(read_buffer, integrableObject);
          }
        } else {
          //molecule belongs to slave nodes
          
          MPI_Recv(&nCurObj, 1, MPI_INT, which_node, TAKE_THIS_TAG_INT,
                   MPI_COMM_WORLD, &istatus);
          
          for(int j = 0; j < nCurObj; j++) {
            eof_test = fgets(read_buffer, sizeof(read_buffer), inFile_);
            
            if (eof_test == NULL) {
              sprintf(painCave.errMsg,
                      "DumpReader Error: error in reading file %s\n"
                      "natoms  = %d; index = %d\n"
                      "error reading the line from the file.\n",
                      filename_.c_str(),
                      nTotObjs,
                      i);
              
              painCave.isFatal = 1;
              simError();
            }
            
            MPI_Send(read_buffer, maxBufferSize, MPI_CHAR, which_node,
                     TAKE_THIS_TAG_CHAR, MPI_COMM_WORLD);
          }
        }
      }
    } else {
      //actions taken at slave nodes
      MPI_Bcast(read_buffer, maxBufferSize, MPI_CHAR, masterNode, MPI_COMM_WORLD);
      
      /**@todo*/
      parseCommentLine(read_buffer, info_->getSnapshotManager()->getCurrentSnapshot());
      
      for(i = 0; i < info_->getNGlobalMolecules(); i++) {
        int which_node = info_->getMolToProc(i);
        
        if (which_node == worldRank) {
          //molecule with global index i belongs to this processor
          
          mol = info_->getMoleculeByGlobalIndex(i);
          if (mol == NULL) {
            sprintf(painCave.errMsg, "DumpReader Error: Molecule not found on node %d!", worldRank);
            painCave.isFatal = 1;
            simError();
          }
          
          nCurObj = mol->getNIntegrableObjects();
          
          MPI_Send(&nCurObj, 1, MPI_INT, masterNode, TAKE_THIS_TAG_INT,
                   MPI_COMM_WORLD);
          
          for (integrableObject = mol->beginIntegrableObject(ii); integrableObject != NULL; 
               integrableObject = mol->nextIntegrableObject(ii)){
            
            MPI_Recv(read_buffer, maxBufferSize, MPI_CHAR, masterNode,
                     TAKE_THIS_TAG_CHAR, MPI_COMM_WORLD, &istatus);
            
            parseDumpLine(read_buffer, integrableObject);
          }
          
        }
        
      }
      
    }
    
#endif
    
  }
  
  void DumpReader::parseDumpLine(char *line, StuntDouble *integrableObject) {
    
    Vector3d pos;  // position place holders
    Vector3d vel;  // velocity placeholders
    Quat4d q;    // the quaternions
    Vector3d ji;   // angular velocity placeholders;
    StringTokenizer tokenizer(line);
    int nTokens;
    
    nTokens = tokenizer.countTokens();
    
    if (nTokens < 14) {
      sprintf(painCave.errMsg,
              "DumpReader Error: Not enough Tokens.\n%s\n", line);
      painCave.isFatal = 1;
      simError();
    }
    
    std::string name = tokenizer.nextToken();
    
    if (name != integrableObject->getType()) {
      
      sprintf(painCave.errMsg,
              "DumpReader Error: Atom type [%s] in %s does not match Atom Type [%s] in .md file.\n",
              name.c_str(), filename_.c_str(), integrableObject->getType().c_str());
      painCave.isFatal = 1;
      simError();        
    }
    
    pos[0] = tokenizer.nextTokenAsDouble();
    pos[1] = tokenizer.nextTokenAsDouble();
    pos[2] = tokenizer.nextTokenAsDouble();
    if (needPos_) {
      integrableObject->setPos(pos);
    }
    
    vel[0] = tokenizer.nextTokenAsDouble();
    vel[1] = tokenizer.nextTokenAsDouble();
    vel[2] = tokenizer.nextTokenAsDouble();
    if (needVel_) {
      integrableObject->setVel(vel);
    }
    
    if (integrableObject->isDirectional()) {
      
      q[0] = tokenizer.nextTokenAsDouble();
      q[1] = tokenizer.nextTokenAsDouble();
      q[2] = tokenizer.nextTokenAsDouble();
      q[3] = tokenizer.nextTokenAsDouble();
      
      double qlen = q.length();
      if (qlen < oopse::epsilon) { //check quaternion is not equal to 0
        
        sprintf(painCave.errMsg,
                "DumpReader Error: initial quaternion error (q0^2 + q1^2 + q2^2 + q3^2 ~ 0).\n");
        painCave.isFatal = 1;
        simError();
        
      } 
      
      q.normalize();
      if (needQuaternion_) {           
        integrableObject->setQ(q);
      }
      
      ji[0] = tokenizer.nextTokenAsDouble();
      ji[1] = tokenizer.nextTokenAsDouble();
      ji[2] = tokenizer.nextTokenAsDouble();
      if (needAngMom_) {
        integrableObject->setJ(ji);
      }
    }
    
  }
  
  
  void DumpReader::parseCommentLine(char* line, Snapshot* s) {
    double currTime;
    Mat3x3d hmat;
    double chi;
    double integralOfChiDt;
    Mat3x3d eta;
    
    StringTokenizer tokenizer(line);
    int nTokens;
    
    nTokens = tokenizer.countTokens();
    
    //comment line should at least contain 10 tokens: current time(1 token) and  h-matrix(9 tokens)
    if (nTokens < 10) {
      sprintf(painCave.errMsg,
              "DumpReader Error: Not enough tokens in comment line: %s", line);
      painCave.isFatal = 1;
      simError();   
    }
    
    //read current time
    currTime = tokenizer.nextTokenAsDouble();
    s->setTime(currTime);
    
    //read h-matrix
    hmat(0, 0) = tokenizer.nextTokenAsDouble();
    hmat(0, 1) = tokenizer.nextTokenAsDouble();
    hmat(0, 2) = tokenizer.nextTokenAsDouble();
    hmat(1, 0) = tokenizer.nextTokenAsDouble();
    hmat(1, 1) = tokenizer.nextTokenAsDouble();
    hmat(1, 2) = tokenizer.nextTokenAsDouble();
    hmat(2, 0) = tokenizer.nextTokenAsDouble();
    hmat(2, 1) = tokenizer.nextTokenAsDouble();
    hmat(2, 2) = tokenizer.nextTokenAsDouble();
    s->setHmat(hmat);
    
    //read chi and integralOfChidt, they should apprear in pair
    if (tokenizer.countTokens() >= 2) {
      chi = tokenizer.nextTokenAsDouble();
      integralOfChiDt = tokenizer.nextTokenAsDouble();            
      
      s->setChi(chi);
      s->setIntegralOfChiDt(integralOfChiDt);
    }
    
    //read eta (eta is 3x3 matrix)
    if (tokenizer.countTokens() >= 9) {
      eta(0, 0) = tokenizer.nextTokenAsDouble();
      eta(0, 1) = tokenizer.nextTokenAsDouble();
      eta(0, 2) = tokenizer.nextTokenAsDouble();
      eta(1, 0) = tokenizer.nextTokenAsDouble();
      eta(1, 1) = tokenizer.nextTokenAsDouble();
      eta(1, 2) = tokenizer.nextTokenAsDouble();
      eta(2, 0) = tokenizer.nextTokenAsDouble();
      eta(2, 1) = tokenizer.nextTokenAsDouble();
      eta(2, 2) = tokenizer.nextTokenAsDouble();      
      
      s->setEta(eta);
    }
    
    
  }
  
}//end namespace oopse
