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

#include <fstream>
#include <iostream>
#include <string>

#include "HydroCmd.hpp"
#include "brains/Register.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "hydrodynamics/AnalyticalModel.hpp"
#include "hydrodynamics/AtomicBeadModel.hpp"
#include "hydrodynamics/BoundaryElementModel.hpp"
#include "hydrodynamics/CompositeShape.hpp"
#include "hydrodynamics/HydroIO.hpp"
#include "hydrodynamics/HydrodynamicsModel.hpp"
#include "hydrodynamics/HydrodynamicsModelCreator.hpp"
#include "hydrodynamics/HydrodynamicsModelFactory.hpp"
#include "hydrodynamics/Mesh.hpp"
#include "hydrodynamics/RoughShell.hpp"
#include "hydrodynamics/ShapeBuilder.hpp"
#include "hydrodynamics/Sphere.hpp"
#include "io/MSMSFormat.hpp"
#include "io/XYZFormat.hpp"
#include "io/stl_reader.h"
#include "utils/ElementsTable.hpp"
#include "utils/MemoryUtils.hpp"
#include "utils/StringTokenizer.hpp"
#include "utils/StringUtils.hpp"
#include "utils/Trim.hpp"
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
  // parse the command line option
  if (cmdline_parser(argc, argv, &args_info) != 0) { exit(1); }

  std::string inFileName;
  std::string modelName;

  bool hasInput = false;
  bool hasModel = false;

  Shape* shape;

  // get the dumpfile name and meta-data file name
  if (args_info.input_given) {
    inFileName = args_info.input_arg;
    hasInput   = true;
  } else {
    if (args_info.stl_given) {
      inFileName = args_info.stl_arg;
      modelName  = "BoundaryElementModel";
      hasInput   = true;
      hasModel   = true;
    } else if (args_info.msms_given) {
      inFileName = args_info.msms_arg;
      modelName  = "BoundaryElementModel";
      hasInput   = true;
      hasModel   = true;
    } else {
      if (args_info.xyz_given) {
        inFileName = args_info.xyz_arg;
        modelName  = "AtomicBeadModel";
        hasInput   = true;
        hasModel   = true;
      }
    }
  }

  if (!hasInput) {
    strcpy(painCave.errMsg, "No input file name was specified.\n");
    painCave.isFatal = 1;
    simError();
  }

  std::string prefix;
  if (args_info.output_given) {
    prefix = args_info.output_arg;
  } else {
    prefix = getPrefix(inFileName);
  }

  if (args_info.model_given) {
    switch (args_info.model_arg) {
    case model_arg_RoughShell:
      modelName = "RoughShell";
      break;
    case model_arg_BoundaryElement:
      modelName = "BoundaryElementModel";
      break;
    case model_arg_AtomicBead:
    default:
      modelName = "AtomicBeadModel";
      break;
    }
  }

  std::string outputFilename = prefix + ".hydro";

  RealType temperature = args_info.temperature_arg;
  RealType viscosity   = args_info.viscosity_arg;

  std::ofstream outputHydro;
  outputHydro.open(outputFilename.c_str());
  HydroIO* hio = new HydroIO();
  hio->openWriter(outputHydro);

  if (args_info.stl_given) {
    std::vector<RealType> coords, normals;
    std::vector<unsigned int> tris, solids;
    shape = new Mesh();

    try {
      stl_reader::ReadStlFile(inFileName.c_str(), coords, normals, tris,
                              solids);
      const size_t numTris = tris.size() / 3;
      for (size_t itri = 0; itri < numTris; ++itri) {
        RealType* c0 = &coords[3 * tris[3 * itri + 0]];
        Vector3d vert0(c0[0], c0[1], c0[2]);
        RealType* c1 = &coords[3 * tris[3 * itri + 1]];
        Vector3d vert1(c1[0], c1[1], c1[2]);
        RealType* c2 = &coords[3 * tris[3 * itri + 2]];
        Vector3d vert2(c2[0], c2[1], c2[2]);
        dynamic_cast<Mesh*>(shape)->add(vert0, vert1, vert2);
      }
    } catch (std::exception& e) { std::cout << e.what() << std::endl; }
    HydrodynamicsModel* model =
        HydrodynamicsModelFactory::getInstance()->createHydrodynamicsModel(
            modelName);

    if (model != NULL) {
      model->init();
      model->setShape(shape);

      HydroProp* hp = model->calcHydroProps(viscosity);
      hio->writeHydroProp(hp, viscosity, temperature, outputHydro);
      hio->interpretHydroProp(hp, viscosity, temperature);
      delete model;

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Could not create HydrodynamicsModel\n");
      painCave.isFatal = 1;
      simError();
    }

    hio->closeWriter(outputHydro);

    outputHydro.close();
  } else if (args_info.msms_given) {
    MSMSFormat* msms = new MSMSFormat(inFileName.c_str());
    shape            = msms->ReadShape();

    HydrodynamicsModel* model =
        HydrodynamicsModelFactory::getInstance()->createHydrodynamicsModel(
            modelName);

    if (model != NULL) {
      model->init();
      model->setShape(shape);

      HydroProp* hp = model->calcHydroProps(viscosity);
      hio->writeHydroProp(hp, viscosity, temperature, outputHydro);
      hio->interpretHydroProp(hp, viscosity, temperature);
      delete model;

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Could not create HydrodynamicsModel\n");
      painCave.isFatal = 1;
      simError();
    }

    hio->closeWriter(outputHydro);

    outputHydro.close();

  } else if (args_info.xyz_given) {
    ifstream in(inFileName);
    if (!in) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Could not open XYZ file\n");
      painCave.isFatal = 1;
      simError();
    }

    XYZFormat* xyz = new XYZFormat();
    xyz->ReadMolecule(in);

    shape = new CompositeShape();
    shape->setName(xyz->title_);

    Shape* currShape = NULL;

    size_t natoms = xyz->mol_.size();
    for (size_t iatom = 0; iatom < natoms; ++iatom) {
      Vector3d pos      = xyz->mol_[iatom]->pos;
      int anum          = xyz->mol_[iatom]->atomicNum;
      std::string atype = xyz->mol_[iatom]->type;
      currShape         = new Sphere(pos, etab.GetVdwRad(anum));
      currShape->setName(atype);
      if (currShape != NULL) {
        dynamic_cast<CompositeShape*>(shape)->addShape(currShape);
      }
    }

    std::cerr << "model = " << modelName << "\n";
    HydrodynamicsModel* model =
        HydrodynamicsModelFactory::getInstance()->createHydrodynamicsModel(
            modelName);

    if (args_info.model_arg == model_arg_RoughShell) {
      RealType bs = args_info.beadSize_arg;
      dynamic_cast<RoughShell*>(model)->setSigma(bs);
    }

    if (model != NULL) {
      model->init();
      model->setShape(shape);

      HydroProp* hp = model->calcHydroProps(viscosity);
      hio->writeHydroProp(hp, viscosity, temperature, outputHydro);
      hio->interpretHydroProp(hp, viscosity, temperature);
      delete model;

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Could not create HydrodynamicsModel\n");
      painCave.isFatal = 1;
      simError();
    }

    hio->closeWriter(outputHydro);

    outputHydro.close();
  } else {
    // parse md file and set up the system
    SimCreator creator;
    SimInfo* info = creator.createSim(inFileName, true);

    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::IntegrableObjectIterator ii;
    StuntDouble* sd;
    Mat3x3d identMat;
    identMat(0, 0) = 1.0;
    identMat(1, 1) = 1.0;
    identMat(2, 2) = 1.0;

    Globals* simParams = info->getSimParams();

    if (simParams->haveViscosity()) {
      viscosity = simParams->getViscosity();
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "viscosity must be set\n");
      painCave.isFatal = 1;
      simError();
    }

    if (simParams->haveTargetTemp()) {
      temperature = simParams->getTargetTemp();
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "target temperature must be set\n");
      painCave.isFatal = 1;
      simError();
    }

    std::map<std::string, SDShape> uniqueStuntDoubles;

    for (mol = info->beginMolecule(mi); mol != NULL;
         mol = info->nextMolecule(mi)) {
      for (sd = mol->beginIntegrableObject(ii); sd != NULL;
           sd = mol->nextIntegrableObject(ii)) {
        if (uniqueStuntDoubles.find(sd->getType()) ==
            uniqueStuntDoubles.end()) {
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

    std::map<std::string, SDShape>::iterator si;
    for (si = uniqueStuntDoubles.begin(); si != uniqueStuntDoubles.end();
         ++si) {
      HydrodynamicsModel* model = NULL;
      shape                     = si->second.shape;
      // StuntDouble* sd = si->second.sd;

      // if (shape->hasAnalyticalSolution()) {
      //  model = new AnalyticalModel(sd, info);
      //} else {
      if (args_info.model_given) {
        std::string modelName;
        switch (args_info.model_arg) {
        case model_arg_RoughShell:
          modelName = "RoughShell";
          break;
        case model_arg_AtomicBead:
        default:
          modelName = "AtomicBeadModel";
          break;
        }

        model =
            HydrodynamicsModelFactory::getInstance()->createHydrodynamicsModel(
                modelName);

        if (args_info.model_arg == model_arg_RoughShell) {
          RealType bs = args_info.beadSize_arg;
          dynamic_cast<RoughShell*>(model)->setSigma(bs);
        }
      }
      //}

      if (model != NULL) {
        model->init();
        model->setShape(shape);

        std::ofstream ofs;
        std::stringstream outputBeads;
        outputBeads << prefix << "_" << shape->getName() << ".xyz";
        ofs.open(outputBeads.str().c_str());
        model->writeElements(ofs);
        ofs.close();

        // if beads option is turned on, skip the calculation
        if (!args_info.beads_flag) {
          HydroProp* hp = model->calcHydroProps(viscosity);
          hio->writeHydroProp(hp, viscosity, temperature, outputHydro);
          hio->interpretHydroProp(hp, viscosity, temperature);
        }

        delete model;

      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Could not create HydrodynamicsModel\n");
        painCave.isFatal = 1;
        simError();
      }
    }
    hio->closeWriter(outputHydro);

    outputHydro.close();
    delete info;
  }
}

void registerHydrodynamicsModels() {
  HydrodynamicsModelFactory::getInstance()->registerHydrodynamicsModel(
      new HydrodynamicsModelBuilder<RoughShell>("RoughShell"));
  HydrodynamicsModelFactory::getInstance()->registerHydrodynamicsModel(
      new HydrodynamicsModelBuilder<AtomicBeadModel>("AtomicBeadModel"));
  HydrodynamicsModelFactory::getInstance()->registerHydrodynamicsModel(
      new HydrodynamicsModelBuilder<AnalyticalModel>("AnalyticalModel"));
  HydrodynamicsModelFactory::getInstance()->registerHydrodynamicsModel(
      new HydrodynamicsModelBuilder<BoundaryElementModel>(
          "BoundaryElementModel"));
}
