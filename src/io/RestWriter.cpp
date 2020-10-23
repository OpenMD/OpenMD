/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#ifdef IS_MPI
#include <mpi.h>
#endif

#include <string>
#include <sstream>
#include <iostream>

#include "io/RestWriter.hpp"
#include "utils/simError.h"
#include "brains/SnapshotManager.hpp"

namespace OpenMD {
  RestWriter::RestWriter(SimInfo* info, const std::string& filename, 
                         std::vector<Restraint*> restraints ) : 
    info_(info){

    std::vector<Restraint*>::const_iterator resti;

    createRestFile_ = false;  

#ifdef IS_MPI    
    MPI_Status* istatus = NULL;
#endif
    
    int printAny = 0;
    for(resti=restraints.begin(); resti != restraints.end(); ++resti){
      if ((*resti)->getPrintRestraint()) {
        printAny = 1;
      }
    }
    
#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &printAny, 1, MPI_INT, MPI_SUM, 
                  MPI_COMM_WORLD);
#endif

    if (printAny) createRestFile_ = true;

#ifdef IS_MPI
    if(worldRank == 0){
#endif
  
      if (createRestFile_) {
        output_ = new std::ofstream(filename.c_str());
        
        if(!output_){
          sprintf( painCave.errMsg,
                   "Could not open %s for restraint output.\n", 
                   filename.c_str());
          painCave.isFatal = 1;
          simError();
        }
      }
        
#ifdef IS_MPI
    }
#endif // is_mpi


#ifndef IS_MPI
         
    if (createRestFile_) (*output_) << "#time\t";

    for(resti=restraints.begin(); resti != restraints.end(); ++resti){
      if ((*resti)->getPrintRestraint()) {
        std::string myName = (*resti)->getRestraintName();
        int myType = (*resti)->getRestraintType();
        
        (*output_) << myName << ":";
        
        if (myType & Restraint::rtDisplacement)
          (*output_) << "\tPosition(angstroms)\tEnergy(kcal/mol)";

        if (myType & Restraint::rtAbsoluteZ)
          (*output_) << "\tPosition(angstroms)\tEnergy(kcal/mol)";
        
        if (myType & Restraint::rtTwist)
          (*output_) << "\tTwistAngle(radians)\tEnergy(kcal/mol)";
        
        if (myType & Restraint::rtSwingX)
          (*output_) << "\tSwingXAngle(radians)\tEnergy(kcal/mol)";
          
        if (myType & Restraint::rtSwingY)
          (*output_) << "\tSwingYAngle(radians)\tEnergy(kcal/mol)";
      }
    }

    if (createRestFile_) (*output_) << "\n";
    if (createRestFile_) (*output_).flush();
    
#else
    
    std::string buffer;

    for(resti=restraints.begin(); resti != restraints.end(); ++resti){
      if ((*resti)->getPrintRestraint()) {
        std::string myName = (*resti)->getRestraintName();
        int myType = (*resti)->getRestraintType();

        buffer += (myName + ":");
        
        if (myType & Restraint::rtDisplacement)
          buffer += "\tPosition(angstroms)\tEnergy(kcal/mol)";

        if (myType & Restraint::rtAbsoluteZ)
          buffer += "\tPosition(angstroms)\tEnergy(kcal/mol)";
        
        if (myType & Restraint::rtTwist)
          buffer += "\tTwistAngle(radians)\tEnergy(kcal/mol)";
        
        if (myType & Restraint::rtSwingX)
          buffer += "\tSwingXAngle(radians)\tEnergy(kcal/mol)";
        
        if (myType & Restraint::rtSwingY)
          buffer += "\tSwingYAngle(radians)\tEnergy(kcal/mol)";
        
        buffer += "\n";
      }
    }
    
    const int primaryNode = 0;
    
    if (worldRank == primaryNode) {
      if (createRestFile_) (*output_) << "#time\t";
      if (createRestFile_) (*output_) << buffer;
      
      int nProc;
      MPI_Comm_size( MPI_COMM_WORLD, &nProc);

      for (int i = 1; i < nProc; ++i) {
        
        // receive the length of the string buffer that was
        // prepared by processor i
        
        int recvLength;
        MPI_Recv(&recvLength, 1, MPI_INT, i, 0, MPI_COMM_WORLD, istatus);
        char* recvBuffer = new char[recvLength];
        if (recvBuffer == NULL) {
        } else {
          MPI_Recv(recvBuffer, recvLength, MPI_CHAR, i, 0, MPI_COMM_WORLD,
                   istatus);
          if (createRestFile_) (*output_) << recvBuffer;
          delete [] recvBuffer;
        }
      }	
       if (createRestFile_) (*output_).flush();
    } else {
      int sendBufferLength = buffer.size() + 1;
      MPI_Send(&sendBufferLength, 1, MPI_INT, primaryNode, 0, MPI_COMM_WORLD);
      MPI_Send((void *)buffer.c_str(), sendBufferLength, MPI_CHAR,
               primaryNode, 0, MPI_COMM_WORLD);
    }
    
#endif // is_mpi    
    
  }    
  
  void RestWriter::writeRest(std::vector<std::map<int, Restraint::RealPair> > restInfo) {
    
#ifdef IS_MPI
    MPI_Status* istatus = NULL;
#endif
    
#ifndef IS_MPI
    if (createRestFile_)  (*output_) << info_->getSnapshotManager()->getCurrentSnapshot()->getTime();
    
    // output some information about the molecules
    std::vector<std::map<int, Restraint::RealPair> >::const_iterator i;
    std::map<int, Restraint::RealPair>::const_iterator j;
    
    if ( createRestFile_ ) {
      
      for( i = restInfo.begin(); i != restInfo.end(); ++i){        
        for(j = (*i).begin(); j != (*i).end(); ++j){                
          (*output_) << "\t" << (j->second).first << "\t" << (j->second).second;
        }
        (*output_) << std::endl;
      }      
      (*output_).flush();
    }
#else
    std::string buffer, first, second;
    std::stringstream ss;
    
    std::vector<std::map<int, Restraint::RealPair> >::const_iterator i;
    std::map<int, Restraint::RealPair>::const_iterator j;
    
    if ( createRestFile_ ) {
      for( i = restInfo.begin(); i != restInfo.end(); ++i){
        
        for(j = (*i).begin(); j != (*i).end(); ++j){
          ss.clear(); 
          ss << (j->second).first;
          ss >> first;
          ss.clear();
          ss << (j->second).second;
          ss >> second;
          buffer += ("\t" + first + "\t" + second);       
        }
        buffer += "\n";     
      }
    }
    
    const int primaryNode = 0;
    
    if (createRestFile_) {
      if (worldRank == primaryNode) {
        
        (*output_) << info_->getSnapshotManager()->getCurrentSnapshot()->getTime();
        (*output_) << buffer;
      
        int nProc;
        MPI_Comm_size( MPI_COMM_WORLD, &nProc);
        for (int i = 1; i < nProc; ++i) {
          
          // receive the length of the string buffer that was
          // prepared by processor i
          
          int recvLength;
          MPI_Recv(&recvLength, 1, MPI_INT, i, 0, MPI_COMM_WORLD, istatus);
          char* recvBuffer = new char[recvLength];
          if (recvBuffer == NULL) {
          } else {
            MPI_Recv(recvBuffer, recvLength, MPI_CHAR, i, 0, MPI_COMM_WORLD,
                     istatus);
            if (createRestFile_) (*output_) << recvBuffer;
            
            delete [] recvBuffer;
          }
        }	
        (*output_).flush();
      } else {
        int sendBufferLength = buffer.size() + 1;
        MPI_Send(&sendBufferLength, 1, MPI_INT, primaryNode, 0, MPI_COMM_WORLD);
        MPI_Send((void *)buffer.c_str(), sendBufferLength, 
                 MPI_CHAR, primaryNode, 0, MPI_COMM_WORLD);
      }
    }
#endif // is_mpi
  }
  
  
  RestWriter::~RestWriter() {
    
#ifdef IS_MPI
    
    if (worldRank == 0) {
#endif // is_mpi
      if (createRestFile_){
        writeClosing(*output_);
        delete output_;
      }
#ifdef IS_MPI 
    }
#endif // is_mpi
  }
  
  void RestWriter::writeClosing(std::ostream& os) {
    os.flush();
  }
  
}// end namespace OpenMD

