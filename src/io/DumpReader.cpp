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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
  
#define _LARGEFILE_SOURCE64 
#define _FILE_OFFSET_BITS 64 
 
#ifdef IS_MPI 
#include <mpi.h> 
#endif
 
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
#include "brains/Thermo.hpp"
 
 
namespace OpenMD { 
   
  DumpReader::DumpReader(SimInfo* info, const std::string& filename) 
    : info_(info), filename_(filename), isScanned_(false), nframes_(0), needCOMprops_(false) { 
    
#ifdef IS_MPI 
    
    if (worldRank == 0) { 
#endif 
      
      inFile_ = new std::ifstream(filename_.c_str(),   
                                  ifstream::in | ifstream::binary); 
      
      if (inFile_->fail()) { 
	sprintf(painCave.errMsg, 
		"DumpReader: Cannot open file: %s\n", 
		filename_.c_str()); 
	painCave.isFatal = 1; 
	simError(); 
      } 
      
#ifdef IS_MPI 
      
    } 
    
    strcpy(checkPointMsg, "Dump file opened for reading successfully."); 
    errorCheckPoint(); 
    
#endif 
    
    return; 
  } 
  
  DumpReader::~DumpReader() { 
    
#ifdef IS_MPI 
    
    if (worldRank == 0) { 
#endif 
      
      delete inFile_; 
      
#ifdef IS_MPI 
      
    } 
    
    strcpy(checkPointMsg, "Dump file closed successfully."); 
    errorCheckPoint(); 
    
#endif 
    
    return; 
  } 
  
  int DumpReader::getNFrames(void) { 
     
    if (!isScanned_) 
      scanFile(); 
     
    return nframes_; 
  } 
   
  void DumpReader::scanFile(void) { 

    std::streampos prevPos;
    std::streampos  currPos; 
    
#ifdef IS_MPI 
    
    if (worldRank == 0) { 
#endif // is_mpi 
      
      currPos = inFile_->tellg();
      prevPos = currPos;
      bool foundOpenSnapshotTag = false;
      bool foundClosedSnapshotTag = false;

      int lineNo = 0; 
      while(inFile_->getline(buffer, bufferSize)) {
        ++lineNo;
        
        std::string line = buffer;
        currPos = inFile_->tellg(); 
        if (line.find("<Snapshot>")!= std::string::npos) {
          if (foundOpenSnapshotTag) {
            sprintf(painCave.errMsg, 
                    "DumpReader:<Snapshot> is multiply nested at line %d in %s \n", lineNo, 
                    filename_.c_str()); 
            painCave.isFatal = 1; 
            simError();           
          }
          foundOpenSnapshotTag = true;
          foundClosedSnapshotTag = false;
          framePos_.push_back(prevPos);
          
        } else if (line.find("</Snapshot>") != std::string::npos){
          if (!foundOpenSnapshotTag) {
            sprintf(painCave.errMsg, 
                    "DumpReader:</Snapshot> appears before <Snapshot> at line %d in %s \n", lineNo, 
                    filename_.c_str()); 
            painCave.isFatal = 1; 
            simError(); 
          }
          
          if (foundClosedSnapshotTag) {
            sprintf(painCave.errMsg, 
                    "DumpReader:</Snapshot> appears multiply nested at line %d in %s \n", lineNo, 
                    filename_.c_str()); 
            painCave.isFatal = 1; 
            simError(); 
          }
          foundClosedSnapshotTag = true;
          foundOpenSnapshotTag = false;
        }
        prevPos = currPos;
      }
      
      // only found <Snapshot> for the last frame means the file is corrupted, we should discard
      // it and give a warning message
      if (foundOpenSnapshotTag) {
        sprintf(painCave.errMsg, 
                "DumpReader: last frame in %s is invalid\n", filename_.c_str()); 
        painCave.isFatal = 0; 
        simError();       
        framePos_.pop_back();
      }
      
      nframes_ = framePos_.size(); 
      
      if (nframes_ == 0) {
        sprintf(painCave.errMsg, 
                "DumpReader: %s does not contain a valid frame\n", filename_.c_str()); 
        painCave.isFatal = 1; 
        simError();      
      }
#ifdef IS_MPI 
    } 
     
    MPI_Bcast(&nframes_, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    
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
     
    if (storageLayout & DataStorage::dslAmat || 
        storageLayout & DataStorage::dslDipole || 
        storageLayout & DataStorage::dslQuadrupole) { 
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

    if (needCOMprops_) {
      Snapshot* s = info_->getSnapshotManager()->getCurrentSnapshot();
      Thermo thermo(info_);
      Vector3d com;

      if (needPos_ && needVel_) {
        Vector3d comvel;
        Vector3d comw;
        thermo.getComAll(com, comvel);
        comw = thermo.getAngularMomentum();
      } else {
        com = thermo.getCom();
      }                    
    }
  } 
   
  void DumpReader::readSet(int whichFrame) {     
    std::string line;

#ifndef IS_MPI 
    inFile_->clear();  
    inFile_->seekg(framePos_[whichFrame]); 

    std::istream& inputStream = *inFile_;     

#else 
    int masterNode = 0;
    std::stringstream sstream;
    if (worldRank == masterNode) {
      std::string sendBuffer;

      inFile_->clear();  
      inFile_->seekg(framePos_[whichFrame]); 
      
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
      MPI_Bcast((void *)sendBuffer.c_str(), sendBufferSize, 
                MPI_CHAR, masterNode, MPI_COMM_WORLD);     
      
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
              "DumpReader Error: can not find <Snapshot>\n"); 
      painCave.isFatal = 1; 
      simError(); 
    } 
    
    //read frameData
    readFrameProperties(inputStream);

    //read StuntDoubles
    readStuntDoubles(inputStream);     

    inputStream.getline(buffer, bufferSize);
    line = buffer;

    if (line.find("<SiteData>") != std::string::npos) {
      //read SiteData
      readSiteData(inputStream);         
    } else {
      if (line.find("</Snapshot>") == std::string::npos) {
        sprintf(painCave.errMsg, 
                "DumpReader Error: can not find </Snapshot>\n"); 
        painCave.isFatal = 1; 
        simError(); 
      }        
    }
  } 
   
  void DumpReader::parseDumpLine(const std::string& line) { 

       
    StringTokenizer tokenizer(line); 
    int nTokens; 
     
    nTokens = tokenizer.countTokens(); 
     
    if (nTokens < 2) {  
      sprintf(painCave.errMsg, 
              "DumpReader Error: Not enough Tokens.\n%s\n", line.c_str()); 
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

    size_t found;
    
    if (needPos_) {
      found = type.find("p");      
      if (found == std::string::npos) {
        sprintf(painCave.errMsg, 
                "DumpReader Error: StuntDouble %d has no Position\n"
                "\tField (\"p\") specified.\n%s\n", index, 
                line.c_str());  
        painCave.isFatal = 1; 
        simError(); 
      }
    }
    
    if (sd->isDirectional()) {
      if (needQuaternion_) {
        found = type.find("q");      
        if (found == std::string::npos) {
          sprintf(painCave.errMsg, 
                  "DumpReader Error: Directional StuntDouble %d has no\n"
                  "\tQuaternion Field (\"q\") specified.\n%s\n", index, 
                  line.c_str());  
          painCave.isFatal = 1; 
          simError(); 
        }
      }      
    }

    for(int i = 0; i < size; ++i) {
      switch(type[i]) {
        
        case 'p': {
            Vector3d pos;
            pos[0] = tokenizer.nextTokenAsDouble(); 
            pos[1] = tokenizer.nextTokenAsDouble(); 
            pos[2] = tokenizer.nextTokenAsDouble(); 
            if (needPos_) { 
              sd->setPos(pos); 
            }             
            break;
        }
        case 'v' : {
            Vector3d vel;
            vel[0] = tokenizer.nextTokenAsDouble(); 
            vel[1] = tokenizer.nextTokenAsDouble(); 
            vel[2] = tokenizer.nextTokenAsDouble(); 
            if (needVel_) { 
              sd->setVel(vel); 
            } 
            break;
        }

        case 'q' : {
           Quat4d q;
           if (sd->isDirectional()) { 
              
             q[0] = tokenizer.nextTokenAsDouble(); 
             q[1] = tokenizer.nextTokenAsDouble(); 
             q[2] = tokenizer.nextTokenAsDouble(); 
             q[3] = tokenizer.nextTokenAsDouble(); 
              
             RealType qlen = q.length(); 
             if (qlen < OpenMD::epsilon) { //check quaternion is not equal to 0 
                
               sprintf(painCave.errMsg, 
                       "DumpReader Error: initial quaternion error (q0^2 + q1^2 + q2^2 + q3^2) ~ 0\n"); 
               painCave.isFatal = 1; 
               simError(); 
                
             }  
              
             q.normalize(); 
             if (needQuaternion_) {            
               sd->setQ(q); 
             }               
           }            
           break;
        }  
        case 'j' : {
          Vector3d ji;
          if (sd->isDirectional()) {
             ji[0] = tokenizer.nextTokenAsDouble(); 
             ji[1] = tokenizer.nextTokenAsDouble(); 
             ji[2] = tokenizer.nextTokenAsDouble(); 
             if (needAngMom_) { 
               sd->setJ(ji); 
             } 
          }
          break;
        }  
        case 'f': {

          Vector3d force;
          force[0] = tokenizer.nextTokenAsDouble(); 
          force[1] = tokenizer.nextTokenAsDouble(); 
          force[2] = tokenizer.nextTokenAsDouble();           
          sd->setFrc(force); 
          break;
        }
        case 't' : {

           Vector3d torque;
           torque[0] = tokenizer.nextTokenAsDouble(); 
           torque[1] = tokenizer.nextTokenAsDouble(); 
           torque[2] = tokenizer.nextTokenAsDouble();           
           sd->setTrq(torque);          
           break;
        }
        case 'u' : {

           RealType particlePot;
           particlePot = tokenizer.nextTokenAsDouble(); 
           sd->setParticlePot(particlePot);          
           break;
        }
        case 'c' : {

           RealType flucQPos;
           flucQPos = tokenizer.nextTokenAsDouble(); 
           sd->setFlucQPos(flucQPos);          
           break;
        }
        case 'w' : {

           RealType flucQVel;
           flucQVel = tokenizer.nextTokenAsDouble(); 
           sd->setFlucQVel(flucQVel);          
           break;
        }
        case 'g' : {

           RealType flucQFrc;
           flucQFrc = tokenizer.nextTokenAsDouble(); 
           sd->setFlucQFrc(flucQFrc);          
           break;
        }
        case 'e' : {

           Vector3d eField;
           eField[0] = tokenizer.nextTokenAsDouble(); 
           eField[1] = tokenizer.nextTokenAsDouble(); 
           eField[2] = tokenizer.nextTokenAsDouble();           
           sd->setElectricField(eField);          
           break;
        }
        case 's' : {

           RealType sPot;
           sPot = tokenizer.nextTokenAsDouble(); 
           sd->setSitePotential(sPot);          
           break;
        }
        default: {
               sprintf(painCave.errMsg, 
                       "DumpReader Error: %s is an unrecognized type\n", type.c_str()); 
               painCave.isFatal = 1; 
               simError(); 
          break;   
        }

      }
    }
  } 
   

  void DumpReader::parseSiteLine(const std::string& line) { 

    StringTokenizer tokenizer(line); 
    int nTokens; 
         
    nTokens = tokenizer.countTokens(); 
     
    if (nTokens < 2) {  
      sprintf(painCave.errMsg, 
              "DumpReader Error: Not enough Tokens.\n%s\n", line.c_str()); 
      painCave.isFatal = 1; 
      simError(); 
    } 

    /**
     * The first token is the global integrable object index.
     */

    int index = tokenizer.nextTokenAsInt();
    StuntDouble* sd = info_->getIOIndexToIntegrableObject(index);
    if (sd == NULL) {
      return;
    }

    /**
     * Test to see if the next token is an integer or not.  If not,
     * we've got data on the integrable object itself.  If there is an
     * integer, we're parsing data for a site on a rigid body.
     */

    std::string indexTest = tokenizer.peekNextToken();
    std::istringstream i(indexTest);
    int siteIndex;
    if (i >> siteIndex) {
      // chew up this token and parse as an int:
      siteIndex = tokenizer.nextTokenAsInt();

      if (sd->isRigidBody()) {
        RigidBody* rb = static_cast<RigidBody*>(sd);
        sd = rb->getAtoms()[siteIndex];
      }
    }

    /**
     * The next token contains information on what follows.
     */
    std::string type = tokenizer.nextToken(); 
    int size = type.size();
    
    for(int i = 0; i < size; ++i) {
      switch(type[i]) {
        
      case 'u' : {
        
        RealType particlePot;
        particlePot = tokenizer.nextTokenAsDouble(); 
        sd->setParticlePot(particlePot);
        break;
      }
      case 'c' : {
        
        RealType flucQPos;
        flucQPos = tokenizer.nextTokenAsDouble(); 
        sd->setFlucQPos(flucQPos);
        break;
      }
      case 'w' : {
        
        RealType flucQVel;
        flucQVel = tokenizer.nextTokenAsDouble(); 
        sd->setFlucQVel(flucQVel);
        break;
      }
      case 'g' : {
        
        RealType flucQFrc;
        flucQFrc = tokenizer.nextTokenAsDouble(); 
        sd->setFlucQFrc(flucQFrc);
        break;
      }
      case 'e' : {
        
        Vector3d eField;
        eField[0] = tokenizer.nextTokenAsDouble(); 
        eField[1] = tokenizer.nextTokenAsDouble(); 
        eField[2] = tokenizer.nextTokenAsDouble();  
        sd->setElectricField(eField);          
        break;
      }
      case 's' : {
        
        RealType sPot;
        sPot = tokenizer.nextTokenAsDouble(); 
        sd->setSitePotential(sPot);          
        break;
      }
      default: {
        sprintf(painCave.errMsg, 
                "DumpReader Error: %s is an unrecognized type\n", type.c_str()); 
        painCave.isFatal = 1; 
        simError(); 
        break;   
      }
      }
    }    
  } 
  
  
    void  DumpReader::readStuntDoubles(std::istream& inputStream) {
    
    inputStream.getline(buffer, bufferSize);
    std::string line(buffer);
    
    if (line.find("<StuntDoubles>") == std::string::npos) {
      sprintf(painCave.errMsg, 
              "DumpReader Error: Missing <StuntDoubles>\n"); 
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

  void  DumpReader::readSiteData(std::istream& inputStream) {

    std::string line(buffer);

    // We already found the starting <SiteData> tag or we wouldn't be
    // here, so just start parsing until we get to the ending
    // </SiteData> tag:
    
    while(inputStream.getline(buffer, bufferSize)) {
      line = buffer;
      
      if(line.find("</SiteData>") != std::string::npos) {
        break;
      }

      parseSiteLine(line);
    }
  
  }

  void DumpReader::readFrameProperties(std::istream& inputStream) {

    Snapshot* s = info_->getSnapshotManager()->getCurrentSnapshot();
    inputStream.getline(buffer, bufferSize);
    std::string line(buffer);

    if (line.find("<FrameData>") == std::string::npos) {
      sprintf(painCave.errMsg, 
              "DumpReader Error: Missing <FrameData>\n"); 
      painCave.isFatal = 1; 
      simError(); 
    }

    while(inputStream.getline(buffer, bufferSize)) {
      line = buffer;
      
      if(line.find("</FrameData>") != std::string::npos) {
        break;
      }
      
      StringTokenizer tokenizer(line, " ;\t\n\r{}:,");
      if (!tokenizer.hasMoreTokens()) {
        sprintf(painCave.errMsg, 
                "DumpReader Error: Not enough Tokens.\n%s\n", line.c_str()); 
        painCave.isFatal = 1; 
        simError();      
      }

      std::string propertyName = tokenizer.nextToken();
      if (propertyName == "Time") {
        RealType currTime = tokenizer.nextTokenAsDouble(); 
        s->setTime(currTime); 
      } else if (propertyName == "Hmat"){
        Mat3x3d hmat;
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
      } else if (propertyName == "Thermostat") {
        pair<RealType, RealType> thermostat;
        thermostat.first = tokenizer.nextTokenAsDouble();
        thermostat.second = tokenizer.nextTokenAsDouble();
        s->setThermostat(thermostat); 
     } else if (propertyName == "Barostat") {
        Mat3x3d eta;
        eta(0, 0) = tokenizer.nextTokenAsDouble(); 
        eta(0, 1) = tokenizer.nextTokenAsDouble(); 
        eta(0, 2) = tokenizer.nextTokenAsDouble(); 
        eta(1, 0) = tokenizer.nextTokenAsDouble(); 
        eta(1, 1) = tokenizer.nextTokenAsDouble(); 
        eta(1, 2) = tokenizer.nextTokenAsDouble(); 
        eta(2, 0) = tokenizer.nextTokenAsDouble(); 
        eta(2, 1) = tokenizer.nextTokenAsDouble(); 
        eta(2, 2) = tokenizer.nextTokenAsDouble(); 
        s->setBarostat(eta); 
      } else {
        sprintf(painCave.errMsg, 
                "DumpReader Error: %s is an invalid property in <FrameData>\n", propertyName.c_str()); 
        painCave.isFatal = 0; 
        simError();         
      }
      
    }

  }

   
}//end namespace OpenMD 
