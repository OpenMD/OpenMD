/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 
#include <algorithm>
#include <iostream>
#include <vector>


#include "io/ZConsWriter.hpp"
#include "utils/simError.h"
#ifdef IS_MPI
#include <mpi.h>
#endif

namespace OpenMD {
  ZConsWriter::ZConsWriter(SimInfo* info, const std::string& filename) : info_(info) {
    //use master - slave mode, only master node writes to disk
#ifdef IS_MPI
    if(worldRank == 0){
#endif

      output_.open(filename.c_str());

      if(!output_){
	sprintf( painCave.errMsg,
		 "Could not open %s for z constrain output_ \n", filename.c_str());
	painCave.isFatal = 1;
	simError();
      }

      output_ << "//time(fs)" << std::endl;
      output_ << "//number of fixed z-constrain molecules" << std::endl;
      output_ << "//global Index of molecule\tzconstrain force\tcurrentZPos" << std::endl;

#ifdef IS_MPI
    }
#endif  

  }

  ZConsWriter::~ZConsWriter()
  {

#ifdef IS_MPI
    if(worldRank == 0 ){
#endif  
      output_.close();  
#ifdef IS_MPI  
    }
#endif
  }

  void ZConsWriter::writeFZ(const std::list<ZconstraintMol>& fixedZmols){
#ifndef IS_MPI
    output_ << info_->getSnapshotManager()->getCurrentSnapshot()->getTime() << std::endl;
    output_ << fixedZmols.size() << std::endl;

    std::list<ZconstraintMol>::const_iterator i;
    for ( i = fixedZmols.begin(); i != fixedZmols.end(); ++i) {
      output_ << i->mol->getGlobalIndex() <<"\t" << i->fz << "\t" << i->zpos << "\t" << i->param.zTargetPos <<std::endl;
    }
#else
    int nproc = MPI::COMM_WORLD.Get_size();
    const int masterNode = 0;
    int myNode = worldRank;
    std::vector<int> tmpNFixedZmols(nproc, 0);
    std::vector<int> nFixedZmolsInProc(nproc, 0);
    tmpNFixedZmols[myNode] = fixedZmols.size();
    
    //do MPI_ALLREDUCE to exchange the total number of atoms,
    //rigidbodies and cutoff groups
    MPI::COMM_WORLD.Allreduce(&tmpNFixedZmols[0], &nFixedZmolsInProc[0], 
                              nproc, MPI::INT, MPI::SUM);

    MPI::Status ierr;
    int zmolIndex;
    RealType data[3];
    
    if (masterNode == 0) {

      std::vector<ZconsData> zconsData;
      ZconsData tmpData;       
      for(int i =0 ; i < nproc; ++i) {
	if (i == masterNode) {
	  std::list<ZconstraintMol>::const_iterator j;
	  for ( j = fixedZmols.begin(); j != fixedZmols.end(); ++j) {
	    tmpData.zmolIndex = j->mol->getGlobalIndex() ;
	    tmpData.zforce= j->fz;
	    tmpData.zpos = j->zpos;
	    tmpData.zconsPos = j->param.zTargetPos;
	    zconsData.push_back(tmpData);
	  }                

	} else {
	  for(int k =0 ; k < nFixedZmolsInProc[i]; ++k) {
            MPI::COMM_WORLD.Recv(&zmolIndex, 1, MPI::INT, i, 0, ierr);
            MPI::COMM_WORLD.Recv(data, 3, MPI::REALTYPE, i, 0, ierr);
	    tmpData.zmolIndex = zmolIndex;
	    tmpData.zforce= data[0];
	    tmpData.zpos = data[1];
	    tmpData.zconsPos = data[2];
	    zconsData.push_back(tmpData);                                        
	  }
	}
            
      }


      output_ << info_->getSnapshotManager()->getCurrentSnapshot()->getTime() << std::endl;
      output_ << zconsData.size() << std::endl;

      std::vector<ZconsData>::iterator l;
      for (l = zconsData.begin(); l != zconsData.end(); ++l) {
	output_ << l->zmolIndex << "\t" << l->zforce << "\t" << l->zpos << "\t" <<  l->zconsPos << std::endl;
      }
        
    } else {

      std::list<ZconstraintMol>::const_iterator j;
      for (j = fixedZmols.begin(); j != fixedZmols.end(); ++j) {
	zmolIndex = j->mol->getGlobalIndex();            
	data[0] = j->fz;
	data[1] = j->zpos;
	data[2] = j->param.zTargetPos;
        MPI::COMM_WORLD.Send(&zmolIndex, 1, MPI::INT, masterNode, 0);
        MPI::COMM_WORLD.Send(data, 3, MPI::REALTYPE, masterNode, 0);
            
      }
    }
#endif
  }

}
