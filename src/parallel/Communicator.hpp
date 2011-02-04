/**
 * @file Communicator.hpp
 * @author Charles Vardeman <cvardema.at.nd.edu>
 * @date 08/18/2010
 * @time 11:56am
 * @version 1.0
 *
 * @section LICENSE
 * Copyright (c) 2010 The University of Notre Dame. All Rights Reserved.
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
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */

#ifndef FORCEDECOMPOSITION_COMMUNICATOR_HPP
#define FORCEDECOMPOSITION_COMMUNICATOR_HPP

#include <config.h>
#include <mpi.h>
#include "math/SquareMatrix3.hpp"

namespace OpenMD{
  
#ifdef IS_MPI

  enum direction {
    Row = 0,
    Column = 1
  };
    
  template<typename T> 
  struct MPITraits
  {
    static const MPI::Datatype datatype;
    static const int dim;
  };
  
  template<> const MPI::Datatype MPITraits<int>::datatype = MPI_INT;
  template<> const int MPITraits<int>::dim = 1;
  template<> const MPI::Datatype MPITraits<RealType>::datatype = MPI_REALTYPE;
  template<> const int MPITraits<RealType>::dim = 1;
  template<> const MPI::Datatype MPITraits<Vector3d>::datatype = MPI_REALTYPE;
  template<> const int MPITraits<Vector3d>::dim = 3;
  template<> const MPI::Datatype MPITraits<Mat3x3d>::datatype = MPI_REALTYPE;
  template<> const int MPITraits<Mat3x3d>::dim = 9;
  
  template<direction D, typename T>
  class Communicator { 
  public: 
    
    Communicator<D, T>(int nObjects) {
      
      int nProc = MPI::COMM_WORLD.Get_size();
      int myRank = MPI::COMM_WORLD.Get_rank();

      int nColumnsMax = (int) sqrt(RealType(nProc));

      int nColumns;
      for (int i = 1; i < nColumnsMax + 1; i++) {
        if (nProc % i == 0) nColumns = i;        
      }
        
      int nRows = nProc / nColumns;
      rowIndex_ = myRank / nColumns;      
      columnIndex_ = myRank % nColumns;

      if (D == Row) {
        myComm = MPI::COMM_WORLD.Split(rowIndex_, 0);
      } else {
        myComm = MPI::COMM_WORLD.Split(columnIndex_, 0);
      }
         
      int nCommProcs = myComm.Get_size();

      counts.reserve(nCommProcs);
      displacements.reserve(nCommProcs);

      planSize_ = MPITraits<T>::dim * nObjects; 

      myComm.Allgather(&planSize_, 1, MPI::INT, &counts[0], 1, MPI::INT);


      displacements[0] = 0;
      for (int i = 1; i < nCommProcs; i++) {
        displacements[i] = displacements[i-1] + counts[i-1];
        size_ += count[i-1];
      }

      size_ = 0;
      for (int i = 0; i < nCommProcs; i++) {
        size_ += counts[i];
      }
    }


    void gather(std::vector<T>& v1, std::vector<T>& v2) {
      
      myComm.Allgatherv(&v1[0], 
                        planSize_, 
                        MPITraits<T>::datatype, 
                        &v2[0], 
                        &counts[0], 
                        &displacements[0], 
                        MPITraits<T>::datatype);      
    }

    
   
    void scatter(std::vector<T>& v1, std::vector<T>& v2) {

      myComm.Reduce_scatter(&v1[0], &v2[0], &counts[0], 
                            MPITraits<T>::datatype, MPI::SUM);
    }

    int getSize() {
      return size_;
    }
    
  private:
    int planSize_;     ///< how many are on local proc
    int rowIndex_;
    int columnIndex_;
    int size_;
    std::vector<int> counts;
    std::vector<int> displacements;
    MPI::Intracomm myComm;
  }; 

#endif
}
#endif


