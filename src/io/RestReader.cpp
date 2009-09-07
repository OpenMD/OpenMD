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
  
#include "io/DumpReader.hpp" 
#include "io/RestReader.hpp"
#include "primitives/Molecule.hpp"
#include "restraints/ObjectRestraint.hpp"
#include "restraints/MolecularRestraint.hpp"

namespace oopse { 
      
  void RestReader::readReferenceStructure() {

    // some of this comes directly from DumpReader.

    if (!isScanned_) 
      scanFile(); 
         
    int storageLayout = info_->getSnapshotManager()->getStorageLayout(); 
     
    if (storageLayout & DataStorage::dslPosition) { 
      needPos_ = true; 
    } else {
      needPos_ = false; 
    } 
     
    needVel_ = false;
     
    if (storageLayout & DataStorage::dslAmat || storageLayout & DataStorage::dslElectroFrame) { 
      needQuaternion_ = true; 
    } else { 
      needQuaternion_ = false; 
    } 

    needAngMom_ = false;

    // We need temporary storage to keep track of all StuntDouble positions
    // in case some of the restraints are molecular (i.e. if they use
    // multiple SD positions to determine restrained orientations or positions:

    all_pos_.clear();
    all_pos_.resize(info_->getNGlobalIntegrableObjects() );


    // Restraint files are just standard dump files, but with the reference
    // structure stored in the first frame (frame 0).
    // RestReader overloads readSet and explicitly handles all of the
    // ObjectRestraints in that method:

    readSet(0); 
    
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
              
              ref.push_back( all_pos_[sd->getGlobalIndex()] );
              
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
              "DumpReader Error: Not enough Tokens.\n%s\n", line.c_str()); 
      painCave.isFatal = 1; 
      simError(); 
    } 

    int index = tokenizer.nextTokenAsInt();
 
    StuntDouble* integrableObject = info_->getIOIndexToIntegrableObject(index);

    if (integrableObject == NULL) {
      return;
    }   

    std::string type = tokenizer.nextToken(); 
    int size = type.size();
    Vector3d pos;
    Vector3d vel;
    Quat4d q;
    Vector3d ji;
    Vector3d force;
    Vector3d torque;
          
    for(int i = 0; i < size; ++i) {
      switch(type[i]) {
        
        case 'p': {
            pos[0] = tokenizer.nextTokenAsDouble(); 
            pos[1] = tokenizer.nextTokenAsDouble(); 
            pos[2] = tokenizer.nextTokenAsDouble(); 
            break;
        }
        case 'v' : {
            vel[0] = tokenizer.nextTokenAsDouble(); 
            vel[1] = tokenizer.nextTokenAsDouble(); 
            vel[2] = tokenizer.nextTokenAsDouble(); 
            break;
        }

        case 'q' : {
           if (integrableObject->isDirectional()) { 
              
             q[0] = tokenizer.nextTokenAsDouble(); 
             q[1] = tokenizer.nextTokenAsDouble(); 
             q[2] = tokenizer.nextTokenAsDouble(); 
             q[3] = tokenizer.nextTokenAsDouble(); 
              
             RealType qlen = q.length(); 
             if (qlen < oopse::epsilon) { //check quaternion is not equal to 0 
                
               sprintf(painCave.errMsg, 
                       "DumpReader Error: initial quaternion error (q0^2 + q1^2 + q2^2 + q3^2) ~ 0\n"); 
               painCave.isFatal = 1; 
               simError(); 
                
             }  
              
             q.normalize(); 
           }            
           break;
        }  
        case 'j' : {
          if (integrableObject->isDirectional()) {
             ji[0] = tokenizer.nextTokenAsDouble(); 
             ji[1] = tokenizer.nextTokenAsDouble(); 
             ji[2] = tokenizer.nextTokenAsDouble(); 
          }
          break;
        }  
        case 'f': {
          force[0] = tokenizer.nextTokenAsDouble(); 
          force[1] = tokenizer.nextTokenAsDouble(); 
          force[2] = tokenizer.nextTokenAsDouble();           

          break;
        }
        case 't' : {
           torque[0] = tokenizer.nextTokenAsDouble(); 
           torque[1] = tokenizer.nextTokenAsDouble(); 
           torque[2] = tokenizer.nextTokenAsDouble();           

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
     
    // keep the position in case we need it for a molecular restraint:
    
    all_pos_[index] = pos;

    // is this io restrained?    
    GenericData* data = integrableObject->getPropertyByName("Restraint");
    ObjectRestraint* oRest;

    if (data != NULL) {
      // make sure we can reinterpret the generic data as restraint data:
      RestraintData* restData= dynamic_cast<RestraintData*>(data);        
      if (restData != NULL) {
        // make sure we can reinterpet the restraint data as a pointer to
        // an ObjectRestraint:
        oRest = dynamic_cast<ObjectRestraint*>(restData->getData());
        if (oRest != NULL) {          
          if (integrableObject->isDirectional()) {
            oRest->setReferenceStructure(pos, q.toRotationMatrix3());
          } else {                           
            oRest->setReferenceStructure(pos);
          }
        }
      }
    }

  }
   


  void RestReader::readFrameProperties(std::istream& inputStream) {
    inputStream.getline(buffer, bufferSize);
    std::string line(buffer);

    if (line.find("<FrameData>") == std::string::npos) {
      sprintf(painCave.errMsg, 
              "DumpReader Error: Missing <FrameData>\n"); 
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

   
}//end namespace oopse 
