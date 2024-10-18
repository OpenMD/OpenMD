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

#include "applications/staticProps/NanoVolume.hpp"
#if defined(HAVE_QHULL)
#include "math/AlphaHull.hpp"
#include "math/ConvexHull.hpp"
#endif
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"

using namespace OpenMD;

NanoVolume::NanoVolume(SimInfo* info, const std::string& filename,
                       const std::string& sele) :
    StaticAnalyser(info, filename, 1),
    selectionScript_(sele), seleMan_(info), evaluator_(info) {
  setOutputName(getPrefix(filename) + ".avol");

  osq.open(getOutputFileName().c_str());

  evaluator_.loadScriptString(sele);
  if (!evaluator_.isDynamic()) {
    seleMan_.setSelectionSet(evaluator_.evaluate());
  }
  frameCounter_ = 0;
}

void NanoVolume::process() {
#if defined(HAVE_QHULL)
  StuntDouble* sd;
  Vector3d vec;
  int i;

  // Do convex hull for now - alpha has issues with perfect structures
  // AlphaHull* thishull = new AlphaHull(2.0);
  ConvexHull* thishull = new ConvexHull();

  DumpReader reader(info_, dumpFilename_);
  int nFrames   = reader.getNFrames();
  frameCounter_ = 0;

  theAtoms_.reserve(info_->getNGlobalAtoms());

  for (int istep = 0; istep < nFrames; istep += step_) {
    reader.readFrame(istep);
    frameCounter_++;
    currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
    RealType time    = currentSnapshot_->getTime();

    // Clear pos vector between each frame.
    theAtoms_.clear();

    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // outer loop is over the selected StuntDoubles:

    for (sd = seleMan_.beginSelected(i); sd != NULL;
         sd = seleMan_.nextSelected(i)) {
      theAtoms_.push_back(sd);
    }

    /* variant below for single atoms, not StuntDoubles:
    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      for (atom = mol->beginAtom(ai); atom != NULL;
           atom = mol->nextAtom(ai)) {
        theAtoms_.push_back(atom);
      }
    }
    */

    // Generate convex hull for this frame.
    thishull->computeHull(theAtoms_);
    RealType volume      = thishull->getVolume();
    RealType surfaceArea = thishull->getArea();

    osq.precision(7);
    if (osq.is_open()) {
      osq << time << "\t" << volume << "\t" << surfaceArea << std::endl;
    }
  }
  osq.close();

#else
  snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
           "NanoVolume: qhull support was not compiled in!\n");
  painCave.isFatal = 1;
  simError();
#endif
}
