/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#ifndef PARALLEL_COMMUNICATOR_HPP
#define PARALLEL_COMMUNICATOR_HPP

#include <config.h>
#include <mpi.h>

#include "math/SquareMatrix3.hpp"

using namespace std;
namespace OpenMD {

#ifdef IS_MPI

  enum communicatorType { Global = 0, Row = 1, Column = 2 };

  template<class T>
  class MPITraits {
  public:
    static MPI_Datatype Type();
    static int Length() { return 1; };
  };

  template<>
  inline MPI_Datatype MPITraits<int>::Type() {
    return MPI_INT;
  }
  template<>
  inline MPI_Datatype MPITraits<RealType>::Type() {
    return MPI_REALTYPE;
  }

  template<class T, unsigned int Dim>
  class MPITraits<Vector<T, Dim>> {
  public:
    static MPI_Datatype Type() { return MPITraits<T>::Type(); }
    static int Length() { return Dim; }
  };

  template<class T>
  class MPITraits<Vector3<T>> {
  public:
    static MPI_Datatype Type() { return MPITraits<T>::Type(); }
    static int Length() { return 3; }
  };

  template<class T, unsigned int R, unsigned int C>
  class MPITraits<RectMatrix<T, R, C>> {
  public:
    static MPI_Datatype Type() { return MPITraits<T>::Type(); }
    static int Length() { return R * C; }
  };

  template<class T>
  class MPITraits<SquareMatrix3<T>> {
  public:
    static MPI_Datatype Type() { return MPITraits<T>::Type(); }
    static int Length() { return 9; }
  };

  template<communicatorType D>
  class Communicator {
  public:
    Communicator<D>() {
      int nProc;
      int myRank;

      MPI_Comm_size(MPI_COMM_WORLD, &nProc);
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

      int nColumnsMax = (int)sqrt(RealType(nProc));

      int nColumns(0);
      for (int i = 1; i < nColumnsMax + 1; i++) {
        if (nProc % i == 0) nColumns = i;
      }

      // int nRows = nProc / nColumns;
      rowIndex_    = myRank / nColumns;
      columnIndex_ = myRank % nColumns;

      switch (D) {
      case Row:
        MPI_Comm_split(MPI_COMM_WORLD, rowIndex_, 0, &myComm);
        break;
      case Column:
        MPI_Comm_split(MPI_COMM_WORLD, columnIndex_, 0, &myComm);
        break;
      case Global:
        MPI_Comm_split(MPI_COMM_WORLD, myRank, 0, &myComm);
      }
    }

    MPI_Comm getComm() { return myComm; }

  private:
    int rowIndex_;
    int columnIndex_;
    MPI_Comm myComm;
  };

  template<typename T>
  class Plan {
  public:
    Plan<T>(MPI_Comm comm, int nObjects) : myComm(comm) {
      int nCommProcs;
      MPI_Comm_size(myComm, &nCommProcs);

      counts.resize(nCommProcs, 0);
      displacements.resize(nCommProcs, 0);

      planSize_ = MPITraits<T>::Length() * nObjects;

      MPI_Allgather(&planSize_, 1, MPI_INT, &counts[0], 1, MPI_INT, myComm);

      displacements[0] = 0;
      for (int i = 1; i < nCommProcs; i++) {
        displacements[i] = displacements[i - 1] + counts[i - 1];
      }

      size_ = 0;
      for (int i = 0; i < nCommProcs; i++) {
        size_ += counts[i];
      }
    }

    void gather(vector<T>& v1, vector<T>& v2) {
      // an assert would be helpful here to make sure the vectors are the
      // correct geometry

      MPI_Allgatherv(&v1[0], planSize_, MPITraits<T>::Type(), &v2[0],
                     &counts[0], &displacements[0], MPITraits<T>::Type(),
                     myComm);
    }

    void scatter(vector<T>& v1, vector<T>& v2) {
      // an assert would be helpful here to make sure the vectors are the
      // correct geometry

      MPI_Reduce_scatter(&v1[0], &v2[0], &counts[0], MPITraits<T>::Type(),
                         MPI_SUM, myComm);
    }

    int getSize() { return size_; }

  private:
    int planSize_;  ///< how many are on local proc
    int size_;
    vector<int> counts;
    vector<int> displacements;
    MPI_Comm myComm;
  };

#endif
}  // namespace OpenMD

#endif
