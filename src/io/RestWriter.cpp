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

#include "io/RestWriter.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"
#include "io/basic_teebuf.hpp"

#ifdef IS_MPI
#include <mpi.h>
#define TAKE_THIS_TAG_INT 1
#define TAKE_THIS_TAG_REAL 2
#endif //is_mpi

namespace oopse {
  RestWriter::RestWriter(SimInfo* info) : 
    info_(info), outName_(info_->getRestFileName()) {
  }
  
  RestWriter::~RestWriter() {}
  
  void RestWriter::writeZAngFile() {
    std::ostream* zangStream;
    
#ifdef IS_MPI
    if (worldRank == 0) {
#endif // is_mpi
      
      zangStream = new std::ofstream(outName_.c_str());
      
#ifdef IS_MPI
    }
#endif // is_mpi    
    
    writeZangle(*zangStream);
    
#ifdef IS_MPI
    if (worldRank == 0) {
#endif // is_mpi
      delete zangStream;
      
#ifdef IS_MPI
    }
#endif // is_mpi  
    
  }

  void RestWriter::writeZangle(std::ostream& finalOut){
    const int BUFFERSIZE = 2000;
    char tempBuffer[BUFFERSIZE];
    char writeLine[BUFFERSIZE];
    
    Molecule* mol;
    StuntDouble* integrableObject;
    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ii;
    
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
    
    MPI_Status ierr;
    int intObIndex;
    int vecLength;
    RealType zAngle;
    std::vector<int> gIndex;
    std::vector<RealType> zValues;

    if (worldRank == masterNode) {
      std::map<int, RealType> zAngData;
      for(int i = 0 ; i < nproc; ++i) {
        if (i == masterNode) {
          for (mol = info_->beginMolecule(mi); mol != NULL; 
               mol = info_->nextMolecule(mi)) {
            
            for (integrableObject = mol->beginIntegrableObject(ii); 
                 integrableObject != NULL; 
                 integrableObject = mol->nextIntegrableObject(ii)) { 
              
              intObIndex = integrableObject->getGlobalIntegrableObjectIndex();

              zAngle = integrableObject->getZangle();
              zAngData.insert(std::pair<int, RealType>(intObIndex, zAngle));
            }      
          }
        } else {
	  MPI_Recv(&vecLength, 1, MPI_INT, i, 
		   TAKE_THIS_TAG_INT, MPI_COMM_WORLD, &ierr);
	  // make sure the vectors are the right size for the incoming data
	  gIndex.resize(vecLength);
	  zValues.resize(vecLength);

	  MPI_Recv(&gIndex[0], vecLength, MPI_INT, i, 
		   TAKE_THIS_TAG_INT, MPI_COMM_WORLD, &ierr);
	  MPI_Recv(&zValues[0], vecLength, MPI_REALTYPE, i, 
		   TAKE_THIS_TAG_REAL, MPI_COMM_WORLD, &ierr);
	  
          for (int k = 0; k < vecLength; k++){
	    zAngData.insert(std::pair<int, RealType>(gIndex[k], zValues[k]));
	  }
	  gIndex.clear();
	  zValues.clear();
        }
      }
      
      finalOut << info_->getSnapshotManager()->getCurrentSnapshot()->getTime()
	       << " : omega values at this time\n";
      
      std::map<int, RealType>::iterator l;
      for (l = zAngData.begin(); l != zAngData.end(); ++l) {

        sprintf( tempBuffer,
                 "%14.10lf\n",
                 l->second);
        strcpy( writeLine, tempBuffer );
        
        finalOut << writeLine;      
      }
      
    } else {
      // pack up and send the appropriate info to the master node
      for(int j = 1; j < nproc; ++j) {
	if (worldRank == j) {
	  for (mol = info_->beginMolecule(mi); mol != NULL; 
	       mol = info_->nextMolecule(mi)) {
	    
	    for (integrableObject = mol->beginIntegrableObject(ii); 
		 integrableObject != NULL; 
		 integrableObject = mol->nextIntegrableObject(ii)) { 
	      
	      // build a vector of the indicies 
	      intObIndex = integrableObject->getGlobalIntegrableObjectIndex();
	      gIndex.push_back(intObIndex);
	      	     
	      // build a vector of the zAngle values
	      zAngle = integrableObject->getZangle();
	      zValues.push_back(zAngle);

	    }      
	  }

	  // let's send these vectors to the master node so that it
	  // can sort them and write to the disk
	  vecLength = gIndex.size();

	  MPI_Send(&vecLength, 1, MPI_INT, masterNode, 
		   TAKE_THIS_TAG_INT, MPI_COMM_WORLD);
	  MPI_Send(&gIndex[0], vecLength, MPI_INT, masterNode, 
		   TAKE_THIS_TAG_INT, MPI_COMM_WORLD);
	  MPI_Send(&zValues[0], vecLength, MPI_REALTYPE, masterNode, 
		   TAKE_THIS_TAG_REAL, MPI_COMM_WORLD);
	
	}
      }
    }

#endif
  }
  
}
