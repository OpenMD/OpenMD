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
#include "stl_reader.h"
#include "utils/ElementsTable.hpp"
#include "utils/MemoryUtils.hpp"
#include "utils/StringTokenizer.hpp"
#include "utils/StringUtils.hpp"
#include "utils/Trim.hpp"
#include "utils/simError.h"

using namespace OpenMD;

void registerHydrodynamicsModels();
void writeHydroProps(std::ostream& os);

int main(int argc, char* argv[]) {
  registerHydrodynamicsModels();

  gengetopt_args_info args_info;
  // parse the command line option
  if (cmdline_parser(argc, argv, &args_info) != 0) { exit(1); }

  std::string inFileName;
  std::string modelName;
  std::string modelSpecified;
  std::string modelRequired;

  bool hasInput = false;
  bool hasModel = false;

  Shape* shape;

  // Check to make sure the model has been given
  if (args_info.model_given) {
    switch (args_info.model_arg) {
    case model_arg_RoughShell:
      modelSpecified = "RoughShell";
      break;
    case model_arg_BoundaryElement:
      modelSpecified = "BoundaryElementModel";
      break;
    case model_arg_AtomicBead:
    default:
      modelSpecified = "AtomicBeadModel";
      break;
    }
    hasModel = true;
  }

  // figure out what kind of input file we have

  if (args_info.input_given) {
    inFileName    = args_info.input_arg;
    modelRequired = modelSpecified;
    hasInput      = true;
  } else if (args_info.stl_given) {
    inFileName    = args_info.stl_arg;
    modelRequired = "BoundaryElementModel";
    hasInput      = true;
  } else if (args_info.msms_given) {
    inFileName    = args_info.msms_arg;
    modelRequired = "BoundaryElementModel";
    hasInput      = true;
  } else if (args_info.xyz_given) {
    inFileName    = args_info.xyz_arg;
    modelRequired = "AtomicBeadModel";
    hasInput      = true;
  }

  if (!hasInput) {
    strcpy(painCave.errMsg, "No input file name was specified.\n");
    painCave.isFatal = 1;
    simError();
  }

  if (hasModel) {
    if (modelSpecified.compare(modelRequired) != 0) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Specified model (%s) does not match model required for input "
               "type (%s).\n",
               modelSpecified.c_str(), modelRequired.c_str());
      painCave.isFatal = 1;
      simError();
    } else {
      modelName = modelSpecified;
    }
  } else {
    modelName = modelRequired;
    hasModel  = true;
  }

  std::string prefix;
  if (args_info.output_given) {
    prefix = args_info.output_arg;
  } else {
    prefix = getPrefix(inFileName);
  }

  std::string outputFilename = prefix + ".hydro";

  RealType temperature = args_info.temperature_arg;
  RealType viscosity   = args_info.viscosity_arg;
  RealType beadSize    = args_info.beadSize_arg;

  // read the files and create the shapes:

  std::map<std::string, Shape*> uniqueShapes;

  if (args_info.stl_given) {
    try {
      stl_reader::StlMesh<RealType, unsigned int> mesh(inFileName.c_str());
      for (size_t isolid = 0; isolid < mesh.num_solids(); ++isolid) {
        shape = new Mesh();
        std::cout << "solid " << isolid << std::endl;
        for (size_t itri = mesh.solid_tris_begin(isolid);
             itri < mesh.solid_tris_end(isolid); ++itri) {
          const RealType* c0 = mesh.tri_corner_coords(itri, 0);
          const RealType* c1 = mesh.tri_corner_coords(itri, 1);
          const RealType* c2 = mesh.tri_corner_coords(itri, 2);
          Vector3d vert0(c0[0], c0[1], c0[2]);
          Vector3d vert1(c1[0], c1[1], c1[2]);
          Vector3d vert2(c2[0], c2[1], c2[2]);
          dynamic_cast<Mesh*>(shape)->add(vert0, vert1, vert2);
        }
        std::string solidName;
        if (mesh.num_solids() > 1) {
          solidName = prefix + "_" + std::to_string(isolid);
        } else {
          solidName = prefix;
        }

        shape->setName(solidName);
        uniqueShapes.insert(
            std::map<std::string, Shape*>::value_type(solidName, shape));
      }
    } catch (std::exception& e) { std::cout << e.what() << std::endl; }

  } else if (args_info.msms_given) {
    MSMSFormat* msms = new MSMSFormat(inFileName.c_str());
    shape            = msms->ReadShape();
    shape->setName(prefix);
    uniqueShapes.insert(
        std::map<std::string, Shape*>::value_type(prefix, shape));
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
    uniqueShapes.insert(
        std::map<std::string, Shape*>::value_type(xyz->title_, shape));
  } else {
    // parse md file and set up the system
    SimCreator creator;
    SimInfo* info = creator.createSim(inFileName, true);

    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::IntegrableObjectIterator ii;
    StuntDouble* sd;

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

    for (mol = info->beginMolecule(mi); mol != NULL;
         mol = info->nextMolecule(mi)) {
      for (sd = mol->beginIntegrableObject(ii); sd != NULL;
           sd = mol->nextIntegrableObject(ii)) {
        if (uniqueShapes.find(sd->getType()) == uniqueShapes.end()) {
          sd->setPos(V3Zero);
          if (sd->isRigidBody()) {
            sd->setA(SquareMatrix3<RealType>::identity());
            RigidBody* rb = static_cast<RigidBody*>(sd);
            rb->updateAtoms();
          }

          Shape* tmp = ShapeBuilder::createShape(sd);
          uniqueShapes.insert(
              std::map<std::string, Shape*>::value_type(sd->getType(), tmp));
        }
      }
    }
    delete info;
  }

  HydrodynamicsModel* model =
      HydrodynamicsModelFactory::getInstance()->createHydrodynamicsModel(
          modelName);

  if (model == NULL) {
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "Could not create HydrodynamicsModel\n");
    painCave.isFatal = 1;
    simError();
  } else {
    model->init();

    if (modelName.compare("RoughShell") == 0) {
      dynamic_cast<RoughShell*>(model)->setSigma(beadSize);
    }

    std::ofstream outputHydro;
    outputHydro.open(outputFilename.c_str());
    HydroIO* hio = new HydroIO();
    hio->openWriter(outputHydro);

    std::map<std::string, Shape*>::iterator si;
    for (si = uniqueShapes.begin(); si != uniqueShapes.end(); ++si) {
      shape = si->second;
      model->setShape(shape);

      // write out the elements making up the shape:
      std::ofstream ofs;
      std::stringstream elementFile;
      elementFile << prefix << "_" << shape->getName();
      if (modelName.compare("BoundaryElementModel") == 0) {
        elementFile << ".stl";
      } else {
        elementFile << ".xyz";
      }
      ofs.open(elementFile.str().c_str());
      model->writeElements(ofs);
      ofs.close();

      // if elements option is turned on, skip the calculation
      if (!args_info.elements_flag) {
        HydroProp* hp = model->calcHydroProps(viscosity);
        hio->writeHydroProp(hp, viscosity, temperature, outputHydro);
        hio->interpretHydroProp(hp, viscosity, temperature);
      }
    }

    hio->closeWriter(outputHydro);
    outputHydro.close();

    delete model;
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
