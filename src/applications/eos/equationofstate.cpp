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

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "brains/ForceManager.hpp"
#include "brains/Register.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "equationofstateCmd.hpp"
#include "flucq/FluctuatingChargeDamped.hpp"
#include "io/DumpReader.hpp"
#include "io/DumpWriter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/Constants.hpp"
#include "utils/ProgressBar.hpp"
#include "utils/simError.h"

using namespace OpenMD;
using namespace std;

int main(int argc, char* argv[]) {
  gengetopt_args_info args_info;
  string omdFileName;
  string outFileName;

  // parse the command line option
  if (cmdline_parser(argc, argv, &args_info) != 0) { exit(1); }

  // get the omd file name and meta-data file name
  if (args_info.input_given) {
    omdFileName = args_info.input_arg;
  } else {
    strcpy(painCave.errMsg, "No input file name was specified.\n");
    painCave.isFatal = 1;
    simError();
  }

  if (args_info.output_given) {
    outFileName = args_info.output_arg;
  } else {
    strcpy(painCave.errMsg, "No output file name was specified.\n");
    painCave.isFatal = 1;
    simError();
  }

  // convert the input angles to radians for computation
  double initial_affine = args_info.initial_arg;
  double final_affine   = args_info.final_arg;
  int number            = args_info.number_arg;

  RealType affine_step = (final_affine - initial_affine) / (number + 1);

  registerAll();

  SimInfo::MoleculeIterator miter;
  Molecule::IntegrableObjectIterator iiter;
  Molecule::RigidBodyIterator rbIter;
  Molecule* mol;
  StuntDouble* sd;
  StuntDouble* sdNew;
  RigidBody* rb;
  Mat3x3d oldHmat;
  Mat3x3d newHmat;
  Snapshot* oldSnap;
  Snapshot* newSnap;
  Vector3d oldPos;
  Vector3d newPos;
  AtomType* atype;

  // parse omd file and set up the system
  SimCreator oldCreator;
  SimInfo* oldInfo = oldCreator.createSim(omdFileName);
  oldSnap          = oldInfo->getSnapshotManager()->getCurrentSnapshot();
  oldHmat          = oldSnap->getHmat();

  // ProgressBarPtr progressBar {nullptr};
  // progressBar = std::make_unique<ProgressBar>();

  ofstream eos;
  eos.open(outFileName.c_str());

  RealType current_affine = initial_affine;
  int countStep           = 0;
  std::cout << "Calculation for EOS started." << std::endl;
  while (current_affine <= final_affine) {
    // progressBar->setStatus(countStep,number);
    // progressBar->update();
    ++countStep;
    RealType scaling    = std::cbrt(current_affine);
    Mat3x3d scaleMatrix = Mat3x3d(0.0);
    scaleMatrix(0, 0)   = scaling;
    scaleMatrix(1, 1)   = scaling;
    scaleMatrix(2, 2)   = scaling;

    SimInfo* newInfo = oldCreator.createSim(omdFileName);
    newSnap          = newInfo->getSnapshotManager()->getCurrentSnapshot();

    newHmat = scaleMatrix * oldHmat;
    newSnap->setHmat(newHmat);

    int newIndex = 0;
    for (mol = oldInfo->beginMolecule(miter); mol != NULL;
         mol = oldInfo->nextMolecule(miter)) {
      for (sd = mol->beginIntegrableObject(iiter); sd != NULL;
           sd = mol->nextIntegrableObject(iiter)) {
        oldPos = sd->getPos();
        oldSnap->wrapVector(oldPos);
        newPos = scaleMatrix * oldPos;
        sdNew  = newInfo->getIOIndexToIntegrableObject(newIndex);
        sdNew->setPos(newPos);
        sdNew->setVel(sd->getVel());
        if (sd->isAtom()) {
          atype                        = static_cast<Atom*>(sd)->getAtomType();
          FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
          if (fqa.isFluctuatingCharge()) {
            RealType charge = sd->getFlucQPos();
            sdNew->setFlucQPos(charge);
            RealType cv = sd->getFlucQVel();
            sdNew->setFlucQVel(cv);
          }
        }
      }

      newIndex++;
    }

    for (mol = newInfo->beginMolecule(miter); mol != NULL;
         mol = newInfo->nextMolecule(miter)) {
      // change the positions of atoms which belong to the rigidbodies
      for (rb = mol->beginRigidBody(rbIter); rb != NULL;
           rb = mol->nextRigidBody(rbIter)) {
        rb->updateAtoms();
        rb->updateAtomVel();
      }
    }

    ForceManager* fman = new ForceManager(newInfo);
    fman->initialize();

    FluctuatingChargePropagator* flucQ = new FluctuatingChargeDamped(newInfo);
    flucQ->setForceManager(fman);
    flucQ->initialize();

    fman->calcForces();
    Thermo thermo(newInfo);
    RealType totalEnergy(0);
    totalEnergy = thermo.getTotalEnergy();
    eos << current_affine << "\t" << totalEnergy << "\n";

    std::cout << countStep << " data generated." << std::endl;

    current_affine += affine_step;
  }
  eos.close();
}
