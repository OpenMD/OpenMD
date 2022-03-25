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

#include "HydroCmd.hpp"
#include "applications/hydrodynamics/AnalyticalModel.hpp"
#include "applications/hydrodynamics/BeadModel.hpp"
#include "applications/hydrodynamics/HydrodynamicsModel.hpp"
#include "applications/hydrodynamics/HydrodynamicsModelCreator.hpp"
#include "applications/hydrodynamics/HydrodynamicsModelFactory.hpp"
#include "applications/hydrodynamics/RoughShell.hpp"
#include "applications/hydrodynamics/ShapeBuilder.hpp"
#include "brains/Register.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "utils/StringUtils.hpp"
#include "utils/simError.h"

using namespace OpenMD;

struct SDShape {
  StuntDouble* sd;
  Shape* shape;
};
void registerHydrodynamicsModels();
void writeHydroProps(std::ostream& os);

int main(int argc, char* argv[]) {
  registerHydrodynamicsModels();

  gengetopt_args_info args_info;
  std::string dumpFileName;
  std::string mdFileName;
  std::string prefix;

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
    prefix = args_info.output_arg;
  } else {
    prefix = getPrefix(dumpFileName);
  }
  std::string outputFilename = prefix + ".hydro";

  // parse md file and set up the system
  SimCreator creator;
  SimInfo* info = creator.createSim(dumpFileName, true);

  SimInfo::MoleculeIterator mi;
  Molecule* mol;
  Molecule::IntegrableObjectIterator ii;
  StuntDouble* sd;
  Mat3x3d identMat;
  identMat(0, 0) = 1.0;
  identMat(1, 1) = 1.0;
  identMat(2, 2) = 1.0;

  Globals* simParams = info->getSimParams();
  RealType temperature(0.0);
  RealType viscosity(0.0);

  if (simParams->haveViscosity()) {
    viscosity = simParams->getViscosity();
  } else {
    sprintf(painCave.errMsg, "viscosity must be set\n");
    painCave.isFatal = 1;
    simError();
  }

  if (simParams->haveTargetTemp()) {
    temperature = simParams->getTargetTemp();
  } else {
    sprintf(painCave.errMsg, "target temperature must be set\n");
    painCave.isFatal = 1;
    simError();
  }

  std::map<std::string, SDShape> uniqueStuntDoubles;

  for (mol = info->beginMolecule(mi); mol != NULL;
       mol = info->nextMolecule(mi)) {
    for (sd = mol->beginIntegrableObject(ii); sd != NULL;
         sd = mol->nextIntegrableObject(ii)) {
      if (uniqueStuntDoubles.find(sd->getType()) == uniqueStuntDoubles.end()) {
        sd->setPos(V3Zero);
        // sd->setA(identMat);
        if (sd->isRigidBody()) {
          sd->setA(identMat);
          RigidBody* rb = static_cast<RigidBody*>(sd);
          rb->updateAtoms();
        }

        SDShape tmp;
        tmp.shape = ShapeBuilder::createShape(sd);
        tmp.sd    = sd;
        uniqueStuntDoubles.insert(
            std::map<std::string, SDShape>::value_type(sd->getType(), tmp));
      }
    }
  }

  std::ofstream outputHydro;
  outputHydro.open(outputFilename.c_str());

  std::map<std::string, SDShape>::iterator si;
  for (si = uniqueStuntDoubles.begin(); si != uniqueStuntDoubles.end(); ++si) {
    HydrodynamicsModel* model = NULL;
    Shape* shape              = si->second.shape;
    StuntDouble* sd           = si->second.sd;

    // if (shape->hasAnalyticalSolution()) {
    //  model = new AnalyticalModel(sd, info);
    //} else {
    if (args_info.model_given) {
      std::string modelName;
      switch (args_info.model_arg) {
      case model_arg_RoughShell:
        modelName = "RoughShell";
        break;
      case model_arg_BeadModel:
      default:
        modelName = "BeadModel";
        break;
      }

      model =
          HydrodynamicsModelFactory::getInstance()->createHydrodynamicsModel(
              modelName, sd, info);

      if (args_info.model_arg == model_arg_RoughShell) {
        RealType bs = args_info.beadSize_arg;
        dynamic_cast<RoughShell*>(model)->setSigma(bs);
      }
    }
    //}

    if (model != NULL) {
      model->init();

      std::ofstream ofs;
      std::stringstream outputBeads;
      outputBeads << prefix << "_" << model->getStuntDoubleName() << ".xyz";
      ofs.open(outputBeads.str().c_str());
      model->writeBeads(ofs);
      ofs.close();

      // if beads option is turned on, skip the calculation
      if (!args_info.beads_flag) {
        model->calcHydroProps(shape, viscosity, temperature);
        model->writeHydroProps(outputHydro);
      }

      delete model;

    } else {
      sprintf(painCave.errMsg, "Could not create HydrodynamicsModel\n");
      painCave.isFatal = 1;
      simError();
    }
  }

  outputHydro.close();

  delete info;
}

void registerHydrodynamicsModels() {
  HydrodynamicsModelFactory::getInstance()->registerHydrodynamicsModel(
      new HydrodynamicsModelBuilder<RoughShell>("RoughShell"));
  HydrodynamicsModelFactory::getInstance()->registerHydrodynamicsModel(
      new HydrodynamicsModelBuilder<BeadModel>("BeadModel"));
  HydrodynamicsModelFactory::getInstance()->registerHydrodynamicsModel(
      new HydrodynamicsModelBuilder<AnalyticalModel>("AnalyticalModel"));
}
