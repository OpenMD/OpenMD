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

#include "io/ConstraintWriter.hpp"

#include <algorithm>
#include <iostream>
#include <vector>

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "utils/simError.h"

namespace OpenMD {
  ConstraintWriter::ConstraintWriter(SimInfo* info,
                                     const std::string& filename) :
      info_(info) {
    // use a primary - secondary model, only the primary node writes
    // to disk
#ifdef IS_MPI
    if (worldRank == 0) {
#endif
      output_.open(filename.c_str());

      if (!output_) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Could not open %s for Constraint output\n", filename.c_str());
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
    if (worldRank == 0) {
#endif
      output_.close();
#ifdef IS_MPI
    }
#endif
  }

  void ConstraintWriter::writeConstraintForces(
      const std::list<ConstraintPair*>& constraints) {
#ifndef IS_MPI
    std::list<ConstraintPair*>::const_iterator i;
    for (i = constraints.begin(); i != constraints.end(); ++i) {
      if ((*i)->getPrintForce()) {
        output_ << info_->getSnapshotManager()->getCurrentSnapshot()->getTime()
                << "\t" << (*i)->getConsElem1()->getGlobalIndex() << "\t"
                << (*i)->getConsElem2()->getGlobalIndex() << "\t"
                << (*i)->getConstraintForce() << std::endl;
      }
    }
#else

    const int primaryNode = 0;
    int nproc;
    int myNode;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myNode);

    std::vector<int> nConstraints(nproc, 0);
    nConstraints[myNode] = constraints.size();

    // do MPI_ALLREDUCE to exchange the total number of constraints:
    MPI_Allreduce(MPI_IN_PLACE, &nConstraints[0], nproc, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);

    MPI_Status ierr;
    int atom1, atom2, doPrint;
    RealType force;

    if (myNode == primaryNode) {
      std::vector<ConstraintData> constraintData;
      ConstraintData tmpData;
      for (int i = 0; i < nproc; ++i) {
        if (i == primaryNode) {
          std::list<ConstraintPair*>::const_iterator j;
          for (j = constraints.begin(); j != constraints.end(); ++j) {
            tmpData.atom1           = (*j)->getConsElem1()->getGlobalIndex();
            tmpData.atom2           = (*j)->getConsElem2()->getGlobalIndex();
            tmpData.constraintForce = (*j)->getConstraintForce();
            tmpData.printForce      = (*j)->getPrintForce();
            constraintData.push_back(tmpData);
          }

        } else {
          for (int k = 0; k < nConstraints[i]; ++k) {
            MPI_Recv(&atom1, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &ierr);
            MPI_Recv(&atom2, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &ierr);
            MPI_Recv(&force, 1, MPI_REALTYPE, i, 0, MPI_COMM_WORLD, &ierr);
            MPI_Recv(&doPrint, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &ierr);

            tmpData.atom1           = atom1;
            tmpData.atom2           = atom2;
            tmpData.constraintForce = force;
            tmpData.printForce      = (doPrint == 0) ? false : true;
            constraintData.push_back(tmpData);
          }
        }
      }

      std::vector<ConstraintData>::iterator l;
      for (l = constraintData.begin(); l != constraintData.end(); ++l) {
        if (l->printForce) {
          output_
              << info_->getSnapshotManager()->getCurrentSnapshot()->getTime()
              << "\t" << l->atom1 << "\t" << l->atom2 << "\t"
              << l->constraintForce << std::endl;
        }
      }
    } else {
      std::list<ConstraintPair*>::const_iterator j;
      for (j = constraints.begin(); j != constraints.end(); ++j) {
        int atom1                = (*j)->getConsElem1()->getGlobalIndex();
        int atom2                = (*j)->getConsElem2()->getGlobalIndex();
        RealType constraintForce = (*j)->getConstraintForce();
        int printForce           = (int)(*j)->getPrintForce();

        MPI_Send(&atom1, 1, MPI_INT, primaryNode, 0, MPI_COMM_WORLD);
        MPI_Send(&atom2, 1, MPI_INT, primaryNode, 0, MPI_COMM_WORLD);
        MPI_Send(&constraintForce, 1, MPI_REALTYPE, primaryNode, 0,
                 MPI_COMM_WORLD);
        MPI_Send(&printForce, 1, MPI_INT, primaryNode, 0, MPI_COMM_WORLD);
      }
    }
#endif
  }
}  // namespace OpenMD
