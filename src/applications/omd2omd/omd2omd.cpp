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

#include <fstream>
#include <iostream>
#include <string>

#include "brains/ForceManager.hpp"
#include "brains/Register.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "io/DumpReader.hpp"
#include "io/DumpWriter.hpp"
#include "math/Quaternion.hpp"
#include "omd2omdCmd.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

using namespace OpenMD;

using namespace std;

void createMdFile(const std::string& oldMdFileName,
                  const std::string& newMdFileName, std::vector<int> nMol);

int main(int argc, char* argv[]) {
  gengetopt_args_info args_info;
  string dumpFileName;
  string outFileName;

  // parse the command line option
  if (cmdline_parser(argc, argv, &args_info) != 0) { exit(1); }

  // get the dumpfile name and meta-data file name
  if (args_info.input_given) {
    dumpFileName = args_info.input_arg;
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
  double phi   = args_info.rotatePhi_arg * (Constants::PI / 180.0);
  double theta = args_info.rotateTheta_arg * (Constants::PI / 180.0);
  double psi   = args_info.rotatePsi_arg * (Constants::PI / 180.0);

  Mat3x3d rotMatrix = Mat3x3d(0.0);

  rotMatrix.setupRotMat(phi, theta, psi);

  Vector3i repeat = Vector3i(args_info.repeatX_arg, args_info.repeatY_arg,
                             args_info.repeatZ_arg);

  Mat3x3d repeatD = Mat3x3d(0.0);
  repeatD(0, 0)   = repeat.x();
  repeatD(1, 1)   = repeat.y();
  repeatD(2, 2)   = repeat.z();

  Vector3d translate =
      Vector3d(args_info.translateX_arg, args_info.translateY_arg,
               args_info.translateZ_arg);

  // parse md file and set up the system

  SimCreator oldCreator;
  SimInfo* oldInfo   = oldCreator.createSim(dumpFileName, false);
  Globals* simParams = oldInfo->getSimParams();
  std::vector<Component*> components = simParams->getComponents();
  std::vector<int> nMol;
  for (vector<Component*>::iterator i = components.begin();
       i != components.end(); ++i) {
    int nMolOld = (*i)->getNMol();
    int nMolNew = nMolOld * repeat.x() * repeat.y() * repeat.z();
    nMol.push_back(nMolNew);
  }

  createMdFile(dumpFileName, outFileName, nMol);

  SimCreator newCreator;
  SimInfo* newInfo = newCreator.createSim(outFileName, false);

  DumpReader* dumpReader = new DumpReader(oldInfo, dumpFileName);
  int nframes            = dumpReader->getNFrames();

  DumpWriter* writer = new DumpWriter(newInfo, outFileName);
  if (writer == NULL) {
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "error in creating DumpWriter");
    painCave.isFatal = 1;
    simError();
  }

  SimInfo::MoleculeIterator miter;
  Molecule::IntegrableObjectIterator iiter;
  Molecule::RigidBodyIterator rbIter;
  Molecule* mol;
  StuntDouble* sd;
  StuntDouble* sdNew;
  RigidBody* rb;
  Mat3x3d oldHmat;
  Mat3x3d rotHmat;
  Mat3x3d newHmat;
  Snapshot* oldSnap;
  Snapshot* newSnap;
  Vector3d oldPos;
  Vector3d newPos;
  Vector3d relPos;
  Vector3d COM;
  Vector3d molCOM;
  Thermo thermo(oldInfo);
  AtomType* atype;

  for (int i = 0; i < nframes; i++) {
    cerr << "frame = " << i << "\n";
    dumpReader->readFrame(i);
    oldSnap = oldInfo->getSnapshotManager()->getCurrentSnapshot();
    newSnap = newInfo->getSnapshotManager()->getCurrentSnapshot();

    newSnap->setID(oldSnap->getID());
    newSnap->setTime(oldSnap->getTime());

    oldHmat = oldSnap->getHmat();
    rotHmat = rotMatrix * oldHmat;
    newHmat = repeatD * rotHmat;
    newSnap->setHmat(newHmat);

    newSnap->setThermostat(oldSnap->getThermostat());
    newSnap->setBarostat(oldSnap->getBarostat());

    COM = thermo.getCom();

    int newIndex = 0;
    for (mol = oldInfo->beginMolecule(miter); mol != NULL;
         mol = oldInfo->nextMolecule(miter)) {
      if (args_info.repairMolecules_arg == 1) { molCOM = mol->getCom(); }

      for (int ii = 0; ii < repeat.x(); ii++) {
        for (int jj = 0; jj < repeat.y(); jj++) {
          for (int kk = 0; kk < repeat.z(); kk++) {
            Vector3d trans = Vector3d(ii, jj, kk);
            for (sd = mol->beginIntegrableObject(iiter); sd != NULL;
                 sd = mol->nextIntegrableObject(iiter)) {
              if (args_info.repairMolecules_arg == 1) {
                relPos = sd->getPos() - molCOM;
                oldPos = molCOM - COM + translate;
                oldSnap->wrapVector(relPos);
                oldSnap->wrapVector(oldPos);
                oldPos += relPos;
              } else {
                oldPos = sd->getPos() - COM + translate;
                oldSnap->wrapVector(oldPos);
              }

              newPos = rotMatrix * oldPos + trans * oldHmat;
              sdNew  = newInfo->getIOIndexToIntegrableObject(newIndex);
              sdNew->setPos(newPos);
              sdNew->setVel(rotMatrix * sd->getVel());

              if (sd->isDirectional()) {
                Mat3x3d bodyRotMat = sd->getA();
                bodyRotMat         = bodyRotMat * rotMatrix.inverse();
                sdNew->setA(bodyRotMat);

                sdNew->setJ(rotMatrix * sd->getJ());
              }

              if (sd->isAtom()) {
                atype = static_cast<Atom*>(sd)->getAtomType();
                FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
                if (fqa.isFluctuatingCharge()) {
                  RealType charge = sd->getFlucQPos();
                  sdNew->setFlucQPos(charge);
                  RealType cv = sd->getFlucQVel();
                  sdNew->setFlucQVel(cv);
                }
              }

              newIndex++;
            }
          }
        }
      }
    }

    // update atoms of rigidbody
    for (mol = newInfo->beginMolecule(miter); mol != NULL;
         mol = newInfo->nextMolecule(miter)) {
      // change the positions of atoms which belong to the rigidbodies
      for (rb = mol->beginRigidBody(rbIter); rb != NULL;
           rb = mol->nextRigidBody(rbIter)) {
        rb->updateAtoms();
        rb->updateAtomVel();
      }
    }

    writer->writeDump();
  }
  // deleting the writer will put the closing at the end of the dump file.
  delete writer;
  delete oldInfo;
}

void createMdFile(const std::string& oldMdFileName,
                  const std::string& newMdFileName, std::vector<int> nMol) {
  ifstream oldMdFile;
  ofstream newMdFile;
  const int MAXLEN = 65535;
  char buffer[MAXLEN];

  // create new .omd file based on old .omd file

  oldMdFile.open(oldMdFileName.c_str());
  newMdFile.open(newMdFileName.c_str());

  oldMdFile.getline(buffer, MAXLEN);

  std::size_t i = 0;
  while (!oldMdFile.eof()) {
    // correct molecule number
    if (strstr(buffer, "nMol") != NULL) {
      if (i < nMol.size()) {
        snprintf(buffer, MAXLEN, "\tnMol = %i;", nMol.at(i));
        newMdFile << buffer << std::endl;
        i++;
      }
    } else
      newMdFile << buffer << std::endl;

    oldMdFile.getline(buffer, MAXLEN);
  }

  oldMdFile.close();
  newMdFile.close();

  if (i != nMol.size()) {
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "Couldn't replace the correct number of nMol\n"
             "\tstatements in component blocks.");
    painCave.isFatal = 1;
    simError();
  }
}
