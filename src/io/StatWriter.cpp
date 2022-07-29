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

#include <config.h>

#define _LARGEFILE_SOURCE64
#define _FILE_OFFSET_BITS 64

#ifdef _MSC_VER
#include <amp_math.h>
#endif

#include "brains/Stats.hpp"
#include "io/StatWriter.hpp"
#include "utils/Revision.hpp"
#include "utils/simError.h"

using namespace std;

namespace OpenMD {

  StatWriter::StatWriter(const std::string& filename, Stats* stats) :
      stats_(stats) {
#ifdef IS_MPI
    if (worldRank == 0) {
#endif

      statfile_.open(filename.c_str(), std::ios::out | std::ios::trunc);

      if (!statfile_) {
        sprintf(painCave.errMsg, "Could not open \"%s\" for stat output.\n",
                filename.c_str());
        painCave.isFatal = 1;
        simError();
      }

      // Create the string for embedding the version information in the MetaData
      version.assign("## Last run using OpenMD version: ");
      version.append(OPENMD_VERSION_MAJOR);
      version.append(".");
      version.append(OPENMD_VERSION_MINOR);
      version.append(",");

      std::string rev(revision, strnlen(revision, 40));
      rev.append(40 - rev.length(), ' ');
      version.append(" revision: ");
      // If there's no GIT revision, just call this the RELEASE revision.
      if (!rev.empty()) {
        version.append(rev);
      } else {
        version.append("RELEASE");
      }

      writeTitle();

#ifdef IS_MPI
    }

    sprintf(checkPointMsg, "Sucessfully opened output file for stating.\n");
    errorCheckPoint();
#endif
  }

  StatWriter::~StatWriter() {
#ifdef IS_MPI
    if (worldRank == 0) {
#endif

      statfile_.close();

#ifdef IS_MPI
    }
#endif
  }

  void StatWriter::writeTitle() {
    Stats::StatsBitSet mask = stats_->getStatsMask();

#ifdef IS_MPI
    if (worldRank == 0) {
#endif

      // write revision
      statfile_ << version << std::endl;
      // write title
      statfile_ << "#";
      for (unsigned int i = 0; i < mask.size(); ++i) {
        if (mask[i]) {
          statfile_ << "\t" << stats_->getTitle(i) << "(" << stats_->getUnits(i)
                    << ")";
        }
      }
      statfile_ << std::endl;

#ifdef IS_MPI
    }
#endif
  }

  void StatWriter::writeStat() {
#ifdef IS_MPI
    if (worldRank == 0) {
#endif

      Stats::StatsBitSet mask = stats_->getStatsMask();
      statfile_.precision(stats_->getPrecision());

      for (unsigned int i = 0; i < mask.size(); ++i) {
        if (mask[i]) {
          if (stats_->getDataType(i) == "RealType")
            writeReal(i);
          else if (stats_->getDataType(i) == "Vector3d")
            writeVector(i);
          else if (stats_->getDataType(i) == "potVec")
            writePotVec(i);
          else if (stats_->getDataType(i) == "Mat3x3d")
            writeMatrix(i);
          else if (stats_->getDataType(i) == "Array2d")
            writeArray(i);
          else {
            sprintf(painCave.errMsg,
                    "StatWriter found an unknown data type for: %s ",
                    stats_->getTitle(i).c_str());
            painCave.isFatal = 1;
            simError();
          }
        }
      }

      statfile_ << std::endl;
      statfile_.rdbuf()->pubsync();

#ifdef IS_MPI
    }
    errorCheckPoint();
#endif
  }

  void StatWriter::writeReal(int i) {
    RealType s = stats_->getRealData(i);

    if (!std::isinf(s) && !std::isnan(s)) {
      statfile_ << "\t" << s;
    } else {
      sprintf(painCave.errMsg,
              "StatWriter detected a numerical error writing: %s ",
              stats_->getTitle(i).c_str());
      painCave.isFatal = 1;
      simError();
    }
  }

  void StatWriter::writeVector(int i) {
    Vector3d s = stats_->getVectorData(i);
    if (std::isinf(s[0]) || std::isnan(s[0]) || std::isinf(s[1]) ||
        std::isnan(s[1]) || std::isinf(s[2]) || std::isnan(s[2])) {
      sprintf(painCave.errMsg,
              "StatWriter detected a numerical error writing: %s",
              stats_->getTitle(i).c_str());
      painCave.isFatal = 1;
      simError();
    } else {
      statfile_ << "\t" << s[0] << "\t" << s[1] << "\t" << s[2];
    }
  }

  void StatWriter::writePotVec(int i) {
    potVec s = stats_->getPotVecData(i);

    bool foundError = false;

    for (unsigned int j = 0; j < N_INTERACTION_FAMILIES; j++) {
      if (std::isinf(s[j]) || std::isnan(s[j])) foundError = true;
    }
    if (foundError) {
      sprintf(painCave.errMsg,
              "StatWriter detected a numerical error writing: %s",
              stats_->getTitle(i).c_str());
      painCave.isFatal = 1;
      simError();
    } else {
      for (unsigned int j = 0; j < N_INTERACTION_FAMILIES; j++) {
        statfile_ << "\t" << s[j];
      }
    }
  }

  void StatWriter::writeMatrix(int i) {
    Mat3x3d s = stats_->getMatrixData(i);

    for (unsigned int i1 = 0; i1 < 3; i1++) {
      for (unsigned int j1 = 0; j1 < 3; j1++) {
        if (std::isinf(s(i1, j1)) || std::isnan(s(i1, j1))) {
          sprintf(painCave.errMsg,
                  "StatWriter detected a numerical error writing: %s",
                  stats_->getTitle(i).c_str());
          painCave.isFatal = 1;
          simError();
        } else {
          statfile_ << "\t" << s(i1, j1);
        }
      }
    }
  }

  void StatWriter::writeArray(int i) {
    std::vector<RealType> s = stats_->getArrayData(i);

    for (unsigned int j = 0; j < s.size(); ++j) {
      if (std::isinf(s[j]) || std::isnan(s[j])) {
        sprintf(painCave.errMsg,
                "StatWriter detected a numerical error writing: %s",
                stats_->getTitle(i).c_str());
        painCave.isFatal = 1;
        simError();
      } else {
        statfile_ << "\t" << s[j];
      }
    }
  }

  void StatWriter::writeStatReport() {
#ifdef IS_MPI
    if (worldRank == 0) {
#endif
      reportfile_.open(reportFileName_.c_str(),
                       std::ios::out | std::ios::trunc);

      if (!reportfile_) {
        sprintf(painCave.errMsg, "Could not open \"%s\" for report output.\n",
                reportFileName_.c_str());
        painCave.isFatal = 1;
        simError();
      }

      reportfile_ << stats_->getStatsReport();
      std::cout << stats_->getStatsReport();
      reportfile_.close();

#ifdef IS_MPI
    }
#endif
  }
}  // namespace OpenMD
