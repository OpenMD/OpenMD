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

#include "applications/sequentialProps/CenterOfMass.hpp"

#include <algorithm>
#include <functional>

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Revision.hpp"
#include "utils/simError.h"

namespace OpenMD {

  CenterOfMass::CenterOfMass(SimInfo* info, const std::string& filename,
                             const std::string& sele1,
                             const std::string& sele2) :
      SequentialAnalyzer(info, filename, sele1, sele2) {
    setOutputName(getPrefix(filename) + ".com");
  }

  void CenterOfMass::doFrame(int) {
    StuntDouble* sd;
    int i;

    if (evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    RealType mtot = 0.0;
    Vector3d com(V3Zero);
    RealType mass;

    for (sd = seleMan1_.beginSelected(i); sd != NULL;
         sd = seleMan1_.nextSelected(i)) {
      mass = sd->getMass();
      mtot += mass;
      com += sd->getPos() * mass;
    }

    com /= mtot;

    values_.push_back(com);
  }

  void CenterOfMass::writeSequence() {
    std::ofstream ofs(outputFilename_.c_str(), std::ios::binary);

    if (ofs.is_open()) {
      Revision r;

      ofs << "# " << getSequenceType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";

      ofs << "#time\tvalue\n";

      for (unsigned int i = 0; i < times_.size(); ++i) {
        ofs << times_[i] << "\t" << values_[i].x() << "\t" << values_[i].y()
            << "\t" << values_[i].z() << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "CenterOfMass::writeSequence Error: failed to open %s\n",
               outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    ofs.close();
  }
}  // namespace OpenMD
