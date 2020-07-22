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

#include <algorithm>
#include <functional>
#include "applications/sequentialProps/equipartitionTest.hpp"
#include "utils/simError.h"
#include "utils/Revision.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"

namespace OpenMD {

  Equipartition::Equipartition(SimInfo* info, const std::string& filename,
                             const std::string& sele1, const std::string& sele2)
    : SequentialAnalyzer(info, filename, sele1, sele2) {

    setOutputName(getPrefix(filename) + ".temp");
  }

  void Equipartition::doFrame(int frame) {
    StuntDouble* sd;
    int i;

    if (evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    const RealType kb = 8.31451e-7;    Vector3d linMom_Temp;
    Vector3d angMom_Temp;
    int count = 0;
    for (sd = seleMan1_.beginSelected(i); sd != NULL;
         sd = seleMan1_.nextSelected(i)) {
        count++;
        Vector3d linMom = sd->getVel() * sd->getMass();
        linMom_Temp[0] += linMom[0]*linMom[0]/(kb * sd->getMass());
        linMom_Temp[1] += linMom[1]*linMom[1]/(kb * sd->getMass());
        linMom_Temp[2] += linMom[2]*linMom[2]/(kb * sd->getMass());
        if (sd->isDirectional()) {
          Vector3d angMom = sd->getJ();
          Mat3x3d momentInertia = sd->getI();
          angMom_Temp[0] += angMom[0] * angMom[0]/(kb * momentInertia(0,0));
          angMom_Temp[1] += angMom[1] * angMom[1]/(kb * momentInertia(1,1));
          angMom_Temp[2] += angMom[2] * angMom[2]/(kb * momentInertia(2,2));

        }



  }

    linMom_Temp /= count;
    angMom_Temp /= count;
    TempP_.push_back( linMom_Temp );
    TempJ_.push_back( angMom_Temp);
  }

  void Equipartition::writeSequence() {
    std::ofstream ofs(outputFilename_.c_str(), std::ios::binary);

    if (ofs.is_open()) {

      Revision r;

      ofs << "# " << getSequenceType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_ ;
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";

      ofs << "#time\t Tpx\t Tpy \t Tpz \t Tjx \t Tjy \t Tjz\n";

      for (unsigned int i = 0; i < times_.size(); ++i) {
        ofs << times_[i]
            << "\t" << TempP_[i].x()
            << "\t" << TempP_[i].y()
            << "\t" << TempP_[i].z()
            << "\t" << TempJ_[i].x()
            << "\t" << TempJ_[i].y()
            << "\t" << TempJ_[i].z()
            << "\n";
      }

    } else {
      sprintf(painCave.errMsg,
              "Equipartition::writeSequence Error: failed to open %s\n",
              outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    ofs.close();
  }
}
