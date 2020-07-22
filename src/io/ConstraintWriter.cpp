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

#include <algorithm>
#include <iostream>
#include <vector>

#include "io/ConstraintWriter.hpp"
#include "utils/simError.h"

namespace OpenMD {
  ConstraintWriter::ConstraintWriter(SimInfo* info, 
                                     const std::string& filename): info_(info) {
    // use a primary - secondary model, only the primary node writes
    // to disk
#ifdef IS_MPI
    if(worldRank == 0){
#endif
      output_.open(filename.c_str());
      
      if(!output_){
	sprintf( painCave.errMsg,
		 "Could not open %s for Constraint output\n", 
                 filename.c_str());
	painCave.isFatal = 1;
	simError();
      }
      
      output_ << "#time(fs)\t"
              << "Index of atom 1\t"
              << "Index of atom 2\tconstraint force" << std::endl; 
      
#ifdef IS_MPI
    }
#endif      
  }
  
  ConstraintWriter::~ConstraintWriter() {    
#ifdef IS_MPI
    if(worldRank == 0 ){
#endif  
      output_.close();  
#ifdef IS_MPI  
    }
#endif
  }
  
  void ConstraintWriter::writeConstraintForces(const std::list<ConstraintPair*>& constraints){
#ifndef IS_MPI
    std::list<ConstraintPair*>::const_iterator i;
    for ( i = constraints.begin(); i != constraints.end(); ++i) {

      if ((*i)->getPrintForce()) {
        output_ << info_->getSnapshotManager()->getCurrentSnapshot()->getTime()  << "\t"
                << (*i)->getConsElem1()->getGlobalIndex() <<"\t" 
                << (*i)->getConsElem2()->getGlobalIndex() <<"\t" 
                << (*i)->getConstraintForce() <<std::endl;
      }
    }
#else
    
    const int primaryNode = 0;
    int nproc;
    int myNode;      
    MPI_Comm_size( MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank( MPI_COMM_WORLD, &myNode);

    std::vector<int> nConstraints(nproc, 0);
    nConstraints[myNode] = constraints.size();
    
    //do MPI_ALLREDUCE to exchange the total number of constraints:
    MPI_Allreduce(MPI_IN_PLACE, &nConstraints[0], nproc, MPI_INT, MPI_SUM, 
                  MPI_COMM_WORLD);
    
    MPI_Status ierr;
    int atom1, atom2, doPrint;
    RealType force;
    
    if (myNode == primaryNode) {
      std::vector<ConstraintData> constraintData;
      ConstraintData tmpData;       
      for(int i = 0 ; i < nproc; ++i) {
        if (i == primaryNode) {
          std::list<ConstraintPair*>::const_iterator j;
          for ( j = constraints.begin(); j != constraints.end(); ++j) {
            tmpData.atom1 = (*j)->getConsElem1()->getGlobalIndex();
            tmpData.atom2 = (*j)->getConsElem2()->getGlobalIndex();         
            tmpData.constraintForce= (*j)->getConstraintForce();
            tmpData.printForce = (*j)->getPrintForce();
            constraintData.push_back(tmpData);
          }                
          
	} else {
	  for(int k = 0; k < nConstraints[i]; ++k) {
	    MPI_Recv(&atom1, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &ierr);
	    MPI_Recv(&atom2, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &ierr);
            MPI_Recv(&force, 1, MPI_REALTYPE, i, 0, MPI_COMM_WORLD, &ierr);
            MPI_Recv(&doPrint, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &ierr);
            
	    tmpData.atom1 = atom1;
	    tmpData.atom2 = atom2;
            tmpData.constraintForce = force;
            tmpData.printForce = (doPrint == 0) ? false : true;
	    constraintData.push_back(tmpData);                                        
	  }
	}            
      }
            
      std::vector<ConstraintData>::iterator l;
      for (l = constraintData.begin(); l != constraintData.end(); ++l) {
        if (l->printForce) {
          output_ << info_->getSnapshotManager()->getCurrentSnapshot()->getTime() << "\t"
                  << l->atom1 <<"\t" 
                  << l->atom2 <<"\t" 
                  << l->constraintForce << std::endl;
        }
      }     
    } else {

      std::list<ConstraintPair*>::const_iterator j;
      for (j = constraints.begin(); j != constraints.end(); ++j) {
        int atom1 = (*j)->getConsElem1()->getGlobalIndex();
        int atom2 = (*j)->getConsElem2()->getGlobalIndex();         
        RealType constraintForce= (*j)->getConstraintForce();
        int printForce = (int) (*j)->getPrintForce();
        
        MPI_Send(&atom1, 1, MPI_INT, primaryNode, 0, MPI_COMM_WORLD);
        MPI_Send(&atom2, 1, MPI_INT, primaryNode, 0, MPI_COMM_WORLD);
        MPI_Send(&constraintForce, 1, MPI_REALTYPE, primaryNode, 0, 
                 MPI_COMM_WORLD);
        MPI_Send(&printForce, 1, MPI_INT, primaryNode, 0, MPI_COMM_WORLD);            
      }
    }
#endif
  }
}
