/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

#include "io/ZConsWriter.hpp"

#include <algorithm>
#include <iostream>
#include <vector>

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "utils/simError.h"

namespace OpenMD {
  ZConsWriter::ZConsWriter(SimInfo* info, const std::string& filename) :
      info_(info) {
    // use a primary - secondary model, only the primary node writes
    // to disk
#ifdef IS_MPI
    if (worldRank == 0) {
#endif
      output_.open(filename.c_str());

      if (!output_) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Could not open %s for z constrain output_ \n",
                 filename.c_str());
        painCave.isFatal = 1;
        simError();
      }

      output_ << "//time(fs)" << std::endl;
      output_ << "//number of fixed z-constrain molecules" << std::endl;
      output_ << "//global Index of molecule\tzconstrain force\tcurrentZPos"
              << std::endl;
#ifdef IS_MPI
    }
#endif
  }

  ZConsWriter::~ZConsWriter() {
#ifdef IS_MPI
    if (worldRank == 0) {
#endif
      output_.close();
#ifdef IS_MPI
    }
#endif
  }

  void ZConsWriter::writeFZ(const std::list<ZconstraintMol>& fixedZmols) {
#ifndef IS_MPI
    output_ << info_->getSnapshotManager()->getCurrentSnapshot()->getTime()
            << std::endl;
    output_ << fixedZmols.size() << std::endl;

    std::list<ZconstraintMol>::const_iterator i;
    for (i = fixedZmols.begin(); i != fixedZmols.end(); ++i) {
      output_ << i->mol->getGlobalIndex() << "\t" << i->fz << "\t" << i->zpos
              << "\t" << i->param.zTargetPos << std::endl;
    }
#else

    const int primaryNode = 0;
    int nproc;
    int myNode;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myNode);

    std::vector<int> tmpNFixedZmols(nproc, 0);
    std::vector<int> nFixedZmolsInProc(nproc, 0);
    tmpNFixedZmols[myNode] = fixedZmols.size();

    // do MPI_ALLREDUCE to exchange the total number of atoms,
    // rigidbodies and cutoff groups
    MPI_Allreduce(&tmpNFixedZmols[0], &nFixedZmolsInProc[0], nproc, MPI_INT,
                  MPI_SUM, MPI_COMM_WORLD);

    MPI_Status* ierr = NULL;
    int zmolIndex;
    RealType data[3];

    if (myNode == primaryNode) {
      std::vector<ZconsData> zconsData;
      ZconsData tmpData;
      for (int i = 0; i < nproc; ++i) {
        if (i == primaryNode) {
          std::list<ZconstraintMol>::const_iterator j;
          for (j = fixedZmols.begin(); j != fixedZmols.end(); ++j) {
            tmpData.zmolIndex = j->mol->getGlobalIndex();
            tmpData.zforce    = j->fz;
            tmpData.zpos      = j->zpos;
            tmpData.zconsPos  = j->param.zTargetPos;
            zconsData.push_back(tmpData);
          }

        } else {
          for (int k = 0; k < nFixedZmolsInProc[i]; ++k) {
            MPI_Recv(&zmolIndex, 1, MPI_INT, i, 0, MPI_COMM_WORLD, ierr);
            MPI_Recv(data, 3, MPI_REALTYPE, i, 0, MPI_COMM_WORLD, ierr);
            tmpData.zmolIndex = zmolIndex;
            tmpData.zforce    = data[0];
            tmpData.zpos      = data[1];
            tmpData.zconsPos  = data[2];
            zconsData.push_back(tmpData);
          }
        }
      }

      output_ << info_->getSnapshotManager()->getCurrentSnapshot()->getTime()
              << std::endl;
      output_ << zconsData.size() << std::endl;

      std::vector<ZconsData>::iterator l;
      for (l = zconsData.begin(); l != zconsData.end(); ++l) {
        output_ << l->zmolIndex << "\t" << l->zforce << "\t" << l->zpos << "\t"
                << l->zconsPos << std::endl;
      }

    } else {
      std::list<ZconstraintMol>::const_iterator j;
      for (j = fixedZmols.begin(); j != fixedZmols.end(); ++j) {
        zmolIndex = j->mol->getGlobalIndex();
        data[0]   = j->fz;
        data[1]   = j->zpos;
        data[2]   = j->param.zTargetPos;
        MPI_Send(&zmolIndex, 1, MPI_INT, primaryNode, 0, MPI_COMM_WORLD);
        MPI_Send(data, 3, MPI_REALTYPE, primaryNode, 0, MPI_COMM_WORLD);
      }
    }
#endif
  }
}  // namespace OpenMD
