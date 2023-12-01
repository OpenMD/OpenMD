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

#include "applications/staticProps/StaticAnalyser.hpp"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <typeindex>
#include <vector>

#include "brains/SimInfo.hpp"
#include "math/Vector3.hpp"
#include "utils/Revision.hpp"
#include "utils/simError.h"

namespace OpenMD {

  StaticAnalyser::StaticAnalyser(SimInfo* info, const std::string& filename,
                                 unsigned int nbins) :
      info_(info),
      dumpFilename_(filename), step_(1), nBins_(nbins) {}

  void StaticAnalyser::writeOutput() {
    std::ofstream ofs(outputFilename_.c_str());

    if (ofs.is_open()) {
      Revision r;
      ofs << "# " << getAnalysisType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      // ofs << "# selection script1: \"" << selectionScript1_ ;
      // ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";

      // write title
      ofs << "#";
      for (const auto& data : data_) {
        if (!data.accumulator.empty()) {
          ofs << "\t" << data.title << "(" << data.units << ")";

          // add some extra tabs for column alignment
          if (data.accumulator[0]->getType() ==
              std::type_index(typeid(Vector3d))) {
            ofs << "\t\t";
          }

          if (data.accumulator[0]->getType() ==
              std::type_index(typeid(std::vector<RealType>))) {
            ofs << "(";
            for (unsigned int type = 0; type < outputTypes_.size(); type++) {
              ofs << outputTypes_[type]->getName() << "\t";
            }
            ofs << ")\t";
          }
        }
      }

      ofs << std::endl;

      ofs.precision(8);

      for (unsigned int bin = 0; bin < nBins_; bin++) {
        for (const auto& data : data_) {
          if (!data.accumulator.empty()) {
            std::string message =
                "StaticAnalyser detected a numerical error writing: " +
                data.title + " for bin " + std::to_string(bin);

            data.accumulator[bin]->writeData(ofs, message, data.dataHandling);
          }
        }
        ofs << std::endl;
      }

      ofs << "#######################################################\n";
      ofs << "# 95% confidence intervals in those quantities follow:\n";
      ofs << "#######################################################\n";

      for (unsigned int bin = 0; bin < nBins_; bin++) {
        ofs << "#";

        for (const auto& data : data_) {
          if (!data.accumulator.empty()) {
            std::string message =
                "StaticAnalyser detected a numerical error writing: " +
                data.title + " std. dev. for bin " + std::to_string(bin);

            data.accumulator[bin]->writeErrorBars(ofs, message);
          }
        }
        ofs << std::endl;
      }

      ofs.flush();
      ofs.close();
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "StaticAnalyser: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }
  }
}  // namespace OpenMD
