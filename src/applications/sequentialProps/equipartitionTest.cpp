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

#include "applications/sequentialProps/equipartitionTest.hpp"

#include <algorithm>
#include <functional>

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Revision.hpp"
#include "utils/simError.h"

namespace OpenMD {

  Equipartition::Equipartition(SimInfo* info, const std::string& filename,
                               const std::string& sele1,
                               const std::string& sele2) :
      SequentialAnalyzer(info, filename, sele1, sele2) {
    setOutputName(getPrefix(filename) + ".temp");
  }

  void Equipartition::doFrame(int) {
    StuntDouble* sd;
    int i;

    if (evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    const RealType kb = 8.31451e-7;
    Vector3d linMom_Temp;
    Vector3d angMom_Temp;
    int count = 0;
    for (sd = seleMan1_.beginSelected(i); sd != NULL;
         sd = seleMan1_.nextSelected(i)) {
      count++;
      Vector3d linMom = sd->getVel() * sd->getMass();
      linMom_Temp[0] += linMom[0] * linMom[0] / (kb * sd->getMass());
      linMom_Temp[1] += linMom[1] * linMom[1] / (kb * sd->getMass());
      linMom_Temp[2] += linMom[2] * linMom[2] / (kb * sd->getMass());
      if (sd->isDirectional()) {
        Vector3d angMom       = sd->getJ();
        Mat3x3d momentInertia = sd->getI();
        angMom_Temp[0] += angMom[0] * angMom[0] / (kb * momentInertia(0, 0));
        angMom_Temp[1] += angMom[1] * angMom[1] / (kb * momentInertia(1, 1));
        angMom_Temp[2] += angMom[2] * angMom[2] / (kb * momentInertia(2, 2));
      }
    }

    linMom_Temp /= count;
    angMom_Temp /= count;
    TempP_.push_back(linMom_Temp);
    TempJ_.push_back(angMom_Temp);
  }

  void Equipartition::writeSequence() {
    std::ofstream ofs(outputFilename_.c_str(), std::ios::binary);

    if (ofs.is_open()) {
      Revision r;

      ofs << "# " << getSequenceType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_;
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";

      ofs << "#time\t Tpx\t Tpy \t Tpz \t Tjx \t Tjy \t Tjz\n";

      for (unsigned int i = 0; i < times_.size(); ++i) {
        ofs << times_[i] << "\t" << TempP_[i].x() << "\t" << TempP_[i].y()
            << "\t" << TempP_[i].z() << "\t" << TempJ_[i].x() << "\t"
            << TempJ_[i].y() << "\t" << TempJ_[i].z() << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Equipartition::writeSequence Error: failed to open %s\n",
               outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    ofs.close();
  }
}  // namespace OpenMD
