/* 
 * Copyright (c) 2009 The University of Notre Dame. All Rights Reserved. 
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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */ 

 
#include <sys/types.h> 
#include <sys/stat.h> 
 
#include <math.h> 
#include <string>
#include <sstream>
#include <iostream> 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
 
#include "io/RestReader.hpp" 
#include "primitives/Molecule.hpp" 
#include "utils/simError.h" 
#include "utils/MemoryUtils.hpp" 
#include "utils/StringTokenizer.hpp"
#include "restraints/ObjectRestraint.hpp"
#include "restraints/MolecularRestraint.hpp" 
 
#ifdef IS_MPI 
 
#include <mpi.h> 
#endif

namespace OpenMD { 

  void RestReader::scanFile(){
    int lineNo = 0; 
    std::streampos prevPos;
    std::streampos  currPos;
    
#ifdef IS_MPI 
    
    if (worldRank == 0) { 
#endif // is_mpi 

      inFile_->clear();
      currPos = inFile_->tellg();
      prevPos = currPos;
      
      bool foundOpenSnapshotTag = false;
      
      while(!foundOpenSnapshotTag && inFile_->getline(buffer, bufferSize)) {
        ++lineNo;
        
        std::string line = buffer;
        currPos = inFile_->tellg(); 
        if (line.find("<Snapshot>")!= std::string::npos) {
          foundOpenSnapshotTag = true;
          framePos_ = prevPos;
        }
        prevPos = currPos;
      }
      
#ifdef IS_MPI 
    }
    MPI_Bcast(&framePos_, 1, MPI_INT, 0, MPI_COMM_WORLD); 
#endif // is_mpi 
  } 


  void RestReader::readSet(){
    std::string line;
    
#ifndef IS_MPI
    
    inFile_->clear();
    inFile_->seekg(framePos_);
    
    std::istream& inputStream = *inFile_;
#else
    
    int masterNode = 0;
    std::stringstream sstream;
    if (worldRank == masterNode) {
      std::string sendBuffer;
      
      inFile_->clear();  
      inFile_->seekg(framePos_); 
      
      while (inFile_->getline(buffer, bufferSize)) {
        
        line = buffer;
        sendBuffer += line;
        sendBuffer += '\n';
        if (line.find("</Snapshot>") != std::string::npos) {
          break;
        }        
      }
      
      int sendBufferSize = sendBuffer.size();
      MPI_Bcast(&sendBufferSize, 1, MPI_INT, masterNode, MPI_COMM_WORLD);     
      MPI_Bcast((void *)sendBuffer.c_str(), sendBufferSize, MPI_CHAR, masterNode, MPI_COMM_WORLD);     
      
      sstream.str(sendBuffer);
    } else {
      int sendBufferSize;
      MPI_Bcast(&sendBufferSize, 1, MPI_INT, masterNode, MPI_COMM_WORLD);     
      char * recvBuffer = new char[sendBufferSize+1];
      assert(recvBuffer);
      recvBuffer[sendBufferSize] = '\0';
      MPI_Bcast(recvBuffer, sendBufferSize, MPI_CHAR, masterNode, MPI_COMM_WORLD);     
      sstream.str(recvBuffer);
      delete [] recvBuffer;
    }      
    
    std::istream& inputStream = sstream;  
#endif
    
    inputStream.getline(buffer, bufferSize);

    line = buffer;
    if (line.find("<Snapshot>") == std::string::npos) {
      sprintf(painCave.errMsg, 
              "RestReader Error: can not find <Snapshot>\n"); 
      painCave.isFatal = 1; 
      simError(); 
    } 
    
    //read frameData
    readFrameProperties(inputStream);
    
    //read StuntDoubles
    readStuntDoubles(inputStream);
    
    inputStream.getline(buffer, bufferSize);
    line = buffer;
    if (line.find("</Snapshot>") == std::string::npos) {
      sprintf(painCave.errMsg, 
              "RestReader Error: can not find </Snapshot>\n"); 
      painCave.isFatal = 1; 
      simError(); 
    }            
  }
    
  void RestReader::readReferenceStructure() {

    // We need temporary storage to keep track of all StuntDouble positions
    // in case some of the restraints are molecular (i.e. if they use
    // multiple SD positions to determine restrained orientations or positions:

    all_pos_.clear();
    all_pos_.resize(info_->getNGlobalIntegrableObjects()) ;
    
    // Restraint files are just standard dump files, but with the reference
    // structure stored in the first frame (frame 0).
    // RestReader overloads readSet and explicitly handles all of the
    // ObjectRestraints in that method:

    scanFile();

    readSet();


    // all ObjectRestraints have been handled, now we have to worry about
    // molecular restraints:

    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule* mol;
    StuntDouble* sd;
    
    // no need to worry about parallel molecules, as molecules are not
    // split across processor boundaries.  Just loop over all molecules
    // we know about:
    
    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {          
      
      // is this molecule restrained?    
      GenericData* data = mol->getPropertyByName("Restraint");
      
      if (data != NULL) {
        
        // make sure we can reinterpret the generic data as restraint data:
        
        RestraintData* restData= dynamic_cast<RestraintData*>(data);        
        
        if (restData != NULL) {
          
          // make sure we can reinterpet the restraint data as a
          // pointer to a MolecularRestraint:
          
          MolecularRestraint* mRest = dynamic_cast<MolecularRestraint*>(restData->getData());
          
          if (mRest != NULL) {          
            
            // now we need to pack the stunt doubles for the reference
            // structure:
            
            std::vector<Vector3d> ref;
            int count = 0;
            RealType mass, mTot;
            Vector3d COM(0.0);
            
            mTot = 0.0;
            
            // loop over the stunt doubles in this molecule in the order we
            // will be looping them in the restraint code:
            
            for (sd = mol->beginIntegrableObject(j); sd != NULL; 
                 sd = mol->nextIntegrableObject(j)) {
              
              // push back the reference positions of the stunt
              // doubles from the *globally* sorted array of
              // positions:
              
              ref.push_back( all_pos_[sd->getGlobalIntegrableObjectIndex()] );
              mass = sd->getMass();              
              COM = COM + mass * ref[count];
              mTot = mTot + mass;
              count = count + 1;
            }
            COM /= mTot;
            mRest->setReferenceStructure(ref, COM);         
          }
        }
      }
    }
  }

   
   
  void RestReader::parseDumpLine(const std::string& line) {        

    StringTokenizer tokenizer(line); 
    int nTokens; 
     
    nTokens = tokenizer.countTokens(); 
     
    if (nTokens < 2) { 
      sprintf(painCave.errMsg, 
              "RestReader Error: Not enough Tokens.\n%s\n", line.c_str()); 
      painCave.isFatal = 1; 
      simError(); 
    } 

    int index = tokenizer.nextTokenAsInt();

    StuntDouble* sd = info_->getIOIndexToIntegrableObject(index);

    if (sd == NULL) {
      return;
    }
  
    std::string type = tokenizer.nextToken(); 
    int size = type.size();
    
    Vector3d pos;
    Quat4d q;

    for(int i = 0; i < size; ++i) {
      switch(type[i]) {

      case 'p': {
        pos[0] = tokenizer.nextTokenAsDouble(); 
        pos[1] = tokenizer.nextTokenAsDouble(); 
        pos[2] = tokenizer.nextTokenAsDouble(); 
        break;
      }
      case 'v' : {
        Vector3d vel;
        vel[0] = tokenizer.nextTokenAsDouble(); 
        vel[1] = tokenizer.nextTokenAsDouble(); 
        vel[2] = tokenizer.nextTokenAsDouble(); 
        break;
      }

      case 'q' : {
        if (sd->isDirectional()) { 
          
          q[0] = tokenizer.nextTokenAsDouble(); 
          q[1] = tokenizer.nextTokenAsDouble(); 
          q[2] = tokenizer.nextTokenAsDouble(); 
          q[3] = tokenizer.nextTokenAsDouble(); 
          
          RealType qlen = q.length(); 
          if (qlen < OpenMD::epsilon) { //check quaternion is not equal to 0 
            
            sprintf(painCave.errMsg, 
                    "RestReader Error: initial quaternion error (q0^2 + q1^2 + q2^2 + q3^2) ~ 0\n"); 
            painCave.isFatal = 1; 
            simError();             
          }  
              
          q.normalize(); 
        }          
        break;
      }  
      case 'j' : {
        Vector3d ji;
        if (sd->isDirectional()) {
          ji[0] = tokenizer.nextTokenAsDouble(); 
          ji[1] = tokenizer.nextTokenAsDouble(); 
          ji[2] = tokenizer.nextTokenAsDouble(); 
        }
        break;
      }  
      case 'f': {        
        Vector3d force;
        force[0] = tokenizer.nextTokenAsDouble(); 
        force[1] = tokenizer.nextTokenAsDouble(); 
        force[2] = tokenizer.nextTokenAsDouble();           
        break;
      }
      case 't' : {        
        Vector3d torque;
        torque[0] = tokenizer.nextTokenAsDouble(); 
        torque[1] = tokenizer.nextTokenAsDouble(); 
        torque[2] = tokenizer.nextTokenAsDouble();           
        break;
      }
      default: {
        sprintf(painCave.errMsg, 
                "RestReader Error: %s is an unrecognized type\n", type.c_str()); 
        painCave.isFatal = 1; 
        simError();     
        break;   
      }
      }
      // keep the position in case we need it for a molecular restraint:

      all_pos_[index] = pos;      
        
      // is this io restrained?
      GenericData* data = sd->getPropertyByName("Restraint");
      ObjectRestraint* oRest;
      
      if (data != NULL) {
        // make sure we can reinterpret the generic data as restraint data:
        RestraintData* restData= dynamic_cast<RestraintData*>(data);        
        if (restData != NULL) {
          // make sure we can reinterpet the restraint data as a pointer to
            // an ObjectRestraint:
          oRest = dynamic_cast<ObjectRestraint*>(restData->getData());
          if (oRest != NULL) {          
            if (sd->isDirectional()) {
              oRest->setReferenceStructure(pos, q.toRotationMatrix3());
            } else {                           
              oRest->setReferenceStructure(pos);
            }
          }
        }
      }
    }
  }
  
  void  RestReader::readStuntDoubles(std::istream& inputStream) {
    
    inputStream.getline(buffer, bufferSize);
    std::string line(buffer);
    
    if (line.find("<StuntDoubles>") == std::string::npos) {
      sprintf(painCave.errMsg, 
              "RestReader Error: Missing <StuntDoubles>\n"); 
      painCave.isFatal = 1; 
      simError(); 
    }
    
    while(inputStream.getline(buffer, bufferSize)) {
      line = buffer;
      
      if(line.find("</StuntDoubles>") != std::string::npos) {
        break;
      }
      
      parseDumpLine(line);
    }
    
  }

  
  void RestReader::readFrameProperties(std::istream& inputStream) {
    inputStream.getline(buffer, bufferSize);
    std::string line(buffer);

    if (line.find("<FrameData>") == std::string::npos) {
      sprintf(painCave.errMsg, 
              "RestReader Error: Missing <FrameData>\n"); 
      painCave.isFatal = 1; 
      simError(); 
    }

    // restraints don't care about frame data (unless we need to wrap
    // coordinates, but we'll worry about that later), so 
    // we'll just scan ahead until the end of the frame data:

    while(inputStream.getline(buffer, bufferSize)) {
      line = buffer;
      
      if(line.find("</FrameData>") != std::string::npos) {
        break;
      }
      
    }
    
  }

   
}//end namespace OpenMD 
