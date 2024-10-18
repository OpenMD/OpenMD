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

#include <fstream>
#include <iostream>
#include <string>

#include "Dump2XYZCmd.hpp"
#include "brains/ForceManager.hpp"
#include "brains/Register.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "io/DumpReader.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "utils/simError.h"
#include "visitors/AtomNameVisitor.hpp"
#include "visitors/AtomVisitor.hpp"
#include "visitors/CompositeVisitor.hpp"
#include "visitors/LipidTransVisitor.hpp"
#include "visitors/OtherVisitor.hpp"
#include "visitors/ReplacementVisitor.hpp"
#include "visitors/RigidBodyVisitor.hpp"
#include "visitors/ZconsVisitor.hpp"

using namespace OpenMD;

using namespace std;
int main(int argc, char* argv[]) {
  gengetopt_args_info args_info;
  string dumpFileName;
  string xyzFileName;

  bool printVel(false);
  bool printFrc(false);
  bool printVec(false);
  bool printChrg(false);
  bool printField(false);
  bool printGlobalID(false);

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
    xyzFileName = args_info.output_arg;
  } else {
    xyzFileName = dumpFileName;
    xyzFileName = xyzFileName.substr(0, xyzFileName.rfind(".")) + ".xyz";
  }

  // parse md file and set up the system
  SimCreator creator;
  SimInfo* info          = creator.createSim(dumpFileName, false);
  ForceManager* forceMan = new ForceManager(info);

  // create visitor list
  CompositeVisitor* compositeVisitor = new CompositeVisitor();

  // create RigidBody Visitor
  if (args_info.rigidbody_flag) {
    RBCOMVisitor* rbCOMVisitor = new RBCOMVisitor(info);
    compositeVisitor->addVisitor(rbCOMVisitor, 900);
  }

  // create SSD atom visitor
  SSDAtomVisitor* ssdVisitor = new SSDAtomVisitor(info);
  compositeVisitor->addVisitor(ssdVisitor, 800);

  // create GBtail atom visitor
  GBtailVisitor* gbtVisitor = new GBtailVisitor(info);
  compositeVisitor->addVisitor(gbtVisitor, 790);

  // create GBhead atom visitor
  GBheadVisitor* gbhVisitor = new GBheadVisitor(info);
  compositeVisitor->addVisitor(gbhVisitor, 789);

  // create default atom visitor
  DefaultAtomVisitor* defaultAtomVisitor = new DefaultAtomVisitor(info);
  compositeVisitor->addVisitor(defaultAtomVisitor, 700);

  // if we gave the -w option, we want to skip the waters:
  if (!args_info.water_given) {
    // create waterType visitor
    if (args_info.watertype_flag) {
      WaterTypeVisitor* waterTypeVisitor = new WaterTypeVisitor;
      compositeVisitor->addVisitor(waterTypeVisitor, 600);
    }
  }

  if (args_info.basetype_flag) {
    AtomNameVisitor* atomNameVisitor = new AtomNameVisitor(info);
    compositeVisitor->addVisitor(atomNameVisitor, 550);
    // When debugging visitors, you may find this helpful:
    //    cout << compositeVisitor->toString();
  }

  // create ZconsVisitor
  if (args_info.zconstraint_flag) {
    ZConsVisitor* zconsVisitor = new ZConsVisitor(info);

    if (zconsVisitor->haveZconsMol()) {
      compositeVisitor->addVisitor(zconsVisitor, 500);
    } else {
      delete zconsVisitor;
    }
  }

  // create wrapping visitor

  // if(args_info.periodicBox_flag){
  //  WrappingVisitor* wrappingVisitor = new WrappingVisitor(info);
  //  compositeVisitor->addVisitor(wrappingVisitor, 400);
  //}

  // create replicate visitor
  if (args_info.repeatX_given > 0 || args_info.repeatY_given > 0 ||
      args_info.repeatZ_given > 0) {
    Vector3i replicateOpt(args_info.repeatX_arg, args_info.repeatY_arg,
                          args_info.repeatZ_arg);
    ReplicateVisitor* replicateVisitor =
        new ReplicateVisitor(info, replicateOpt);
    compositeVisitor->addVisitor(replicateVisitor, 300);
  }

  // create rotation visitor
  if (args_info.refsele_given && args_info.originsele_given) {
    compositeVisitor->addVisitor(
        new LipidTransVisitor(info, args_info.originsele_arg,
                              args_info.refsele_arg),
        250);
  } else if (args_info.refsele_given || args_info.originsele_given) {
    strcpy(
        painCave.errMsg,
        "The --refsele and --originsele arguments should appear together.\n");
    painCave.isFatal = 1;
    simError();
  }

  // create xyzVisitor
  XYZVisitor* xyzVisitor;

  if (args_info.selection_given) {
    xyzVisitor = new XYZVisitor(info, args_info.selection_arg);
  } else {
    xyzVisitor = new XYZVisitor(info);
  }

  if (args_info.velocities_flag) {
    printVel = true;
    xyzVisitor->doVelocities(printVel);
  }
  if (args_info.forces_flag) {
    printFrc = true;
    xyzVisitor->doForces(printFrc);
  }
  if (args_info.vectors_flag) {
    printVec = true;
    xyzVisitor->doVectors(printVec);
  }
  if (args_info.charges_flag) {
    printChrg = true;
    xyzVisitor->doCharges(printChrg);
  }
  if (args_info.efield_flag) {
    printField = true;
    xyzVisitor->doElectricFields(printField);
  }
  if (args_info.globalID_flag) {
    printGlobalID = true;
    xyzVisitor->doGlobalIDs(printGlobalID);
  }

  compositeVisitor->addVisitor(xyzVisitor, 200);

  // create prepareVisitor
  PrepareVisitor* prepareVisitor = new PrepareVisitor();

  // open dump file
  DumpReader* dumpReader = new DumpReader(info, dumpFileName);
  int nframes            = dumpReader->getNFrames();

  ofstream xyzStream(xyzFileName.c_str());

  SimInfo::MoleculeIterator miter;
  Molecule::IntegrableObjectIterator iiter;
  Molecule::RigidBodyIterator rbIter;
  Molecule* mol;
  StuntDouble* sd;
  RigidBody* rb;
  Vector3d molCom;
  Vector3d newMolCom;
  Vector3d displacement;
  Mat3x3d hmat;
  Snapshot* currentSnapshot;

  for (int i = 0; i < nframes; i += args_info.frame_arg) {
    dumpReader->readFrame(i);

    if (printFrc) forceMan->calcForces();

    // wrapping the molecule
    if (args_info.periodicBox_flag) {
      currentSnapshot = info->getSnapshotManager()->getCurrentSnapshot();
      for (mol = info->beginMolecule(miter); mol != NULL;
           mol = info->nextMolecule(miter)) {
        molCom    = mol->getCom();
        newMolCom = molCom;
        currentSnapshot->wrapVector(newMolCom);
        displacement = newMolCom - molCom;

        for (sd = mol->beginIntegrableObject(iiter); sd != NULL;
             sd = mol->nextIntegrableObject(iiter)) {
          sd->setPos(sd->getPos() + displacement);
        }
      }
    }

    // update atoms of rigidbody
    for (mol = info->beginMolecule(miter); mol != NULL;
         mol = info->nextMolecule(miter)) {
      // change the positions of atoms which belong to the rigidbodies
      for (rb = mol->beginRigidBody(rbIter); rb != NULL;
           rb = mol->nextRigidBody(rbIter)) {
        rb->updateAtoms();
        if (printVel) rb->updateAtomVel();
      }
    }

    // prepare visit
    for (mol = info->beginMolecule(miter); mol != NULL;
         mol = info->nextMolecule(miter)) {
      for (sd = mol->beginIntegrableObject(iiter); sd != NULL;
           sd = mol->nextIntegrableObject(iiter)) {
        sd->accept(prepareVisitor);
      }
    }

    // update visitor
    compositeVisitor->update();

    // visit stuntdouble
    for (mol = info->beginMolecule(miter); mol != NULL;
         mol = info->nextMolecule(miter)) {
      for (sd = mol->beginIntegrableObject(iiter); sd != NULL;
           sd = mol->nextIntegrableObject(iiter)) {
        sd->accept(compositeVisitor);
      }
    }

    xyzVisitor->writeFrame(xyzStream);
    xyzVisitor->clear();

  }  // end for (int i = 0; i < nframes; i += args_info.frame_arg)

  xyzStream.close();

  delete forceMan;
  delete compositeVisitor;
  delete prepareVisitor;
  delete dumpReader;

  delete info;
}
