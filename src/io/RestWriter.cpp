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
 
#include <algorithm>
#include <iostream>
#include <map>

#include "primitives/Molecule.hpp"
#include "io/RestWriter.hpp"
#include "utils/simError.h"


namespace oopse {
  RestWriter::RestWriter(SimInfo* info) : 
  info_(info) {
    
    //we use master - slave mode, only master node writes to disk
    outName = info_->getRestFileName();
  }
  
  RestWriter::~RestWriter() {}
  
  void RestWriter::writeZangle(){
    const int BUFFERSIZE = 2000;
    char tempBuffer[BUFFERSIZE];
    char writeLine[BUFFERSIZE];
    
    std::ofstream finalOut;
    
    Molecule* mol;
    StuntDouble* integrableObject;
    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ii;
    
#ifdef IS_MPI
    if(worldRank == 0 ){
#endif    
      finalOut.open( outName.c_str(), std::ios::out | std::ios::trunc );
      if( !finalOut ){
        sprintf( painCave.errMsg,
                 "Could not open \"%s\" for zAngle output.\n",
                 outName.c_str() );
        painCave.isFatal = 1;
        simError();
      }
#ifdef IS_MPI
    }
#endif // is_mpi
    
#ifndef IS_MPI
    // first we do output for the single processor version
    finalOut
      << info_->getSnapshotManager()->getCurrentSnapshot()->getTime()
      << " : omega values at this time\n";
    
    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) {
      
      for (integrableObject = mol->beginIntegrableObject(ii); 
           integrableObject != NULL; 
           integrableObject = mol->nextIntegrableObject(ii)) {    
        
        sprintf( tempBuffer,
                 "%14.10lf\n",
                 integrableObject->getZangle());
        strcpy( writeLine, tempBuffer );
        
        finalOut << writeLine;      
      }
      
    }
    
#else
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    const int masterNode = 0;
    int myNode = worldRank;
    std::vector<int> tmpNIntObjects(nproc, 0);
    std::vector<int> nIntObjectsInProc(nproc, 0);
    tmpNIntObjects[myNode] = info_->getNGlobalIntegrableObjects();
    
    //do MPI_ALLREDUCE to exchange the total number of atoms, rigidbodies and cutoff groups
    MPI_Allreduce(&tmpNIntObjects[0], &nIntObjectsInProc[0], nproc, MPI_INT,
                  MPI_SUM, MPI_COMM_WORLD);
    
    MPI_Status ierr;
    int intObIndex;
    double zAngle;
   
    if (masterNode == 0) {
      std::map<int, double> zAngData;
      for(int i = 0 ; i < nproc; ++i) {
        if (i == masterNode) {
          for (mol = info_->beginMolecule(mi); mol != NULL; 
               mol = info_->nextMolecule(mi)) {
            
            for (integrableObject = mol->beginIntegrableObject(ii); 
                 integrableObject != NULL; 
                 integrableObject = mol->nextIntegrableObject(ii)) { 
              
              intObIndex = integrableObject->getGlobalIndex() ;
              zAngle = integrableObject->getZangle();
              zAngData.insert(std::pair<int, double>(intObIndex, zAngle));
            }      
          }
          
        } else {
          for(int k = 0; k < nIntObjectsInProc[i]; ++k) {
            MPI_Recv(&intObIndex, 1, MPI_INT, i, 0, MPI_COMM_WORLD,&ierr);
            MPI_Recv(&zAngle, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,&ierr);
            zAngData.insert(std::pair<int, double>(intObIndex, zAngle));
          }
        }
        
      }
      
      finalOut
        << info_->getSnapshotManager()->getCurrentSnapshot()->getTime()
        << " : omega values at this time\n";
      
      std::map<int, double>::iterator l;
      for (l = zAngData.begin(); l != zAngData.end(); ++l) {
        finalOut << l->second << "\n";
      }
      
    } else {
      
      for (mol = info_->beginMolecule(mi); mol != NULL; 
           mol = info_->nextMolecule(mi)) {
        
        for (integrableObject = mol->beginIntegrableObject(ii); 
             integrableObject != NULL; 
             integrableObject = mol->nextIntegrableObject(ii)) { 
          intObIndex = integrableObject->getGlobalIndex();            
          zAngle = integrableObject->getZangle();
          MPI_Send(&intObIndex, 1, MPI_INT, masterNode, 0, MPI_COMM_WORLD);
          MPI_Send(&zAngle, 1, MPI_DOUBLE, masterNode, 0, MPI_COMM_WORLD);
        }
      }
    }
#endif
    
#ifdef IS_MPI
    finalOut.close();
#endif
    
  }
  
}
