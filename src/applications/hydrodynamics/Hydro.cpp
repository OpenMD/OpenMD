/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
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

#include <fstream>
#include <iostream>
#include <string>

#include "HydroCmd.hpp"
#include "hydrodynamics/AnalyticalModel.hpp"
#include "hydrodynamics/AtomicBeadModel.hpp"
#include "hydrodynamics/HydrodynamicsModel.hpp"
#include "hydrodynamics/HydrodynamicsModelCreator.hpp"
#include "hydrodynamics/HydrodynamicsModelFactory.hpp"
#include "hydrodynamics/RoughShell.hpp"
#include "hydrodynamics/ShapeBuilder.hpp"
#include "hydrodynamics/HydroIO.hpp"
#include "hydrodynamics/Mesh.hpp"
#include "brains/Register.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "utils/MemoryUtils.hpp"
#include "utils/StringUtils.hpp"
#include "utils/simError.h"

#include "io/stl_reader.h"

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
  bool hasViscosity = false;
  bool hasTemperature = false;
  
  // get the dumpfile name and meta-data file name
  if (args_info.input_given) {
    inFileName = args_info.input_arg;
    hasInput = true;
  }  else {
    if (args_info.stl_given) {
      inFileName = args_info.stl_arg;
      modelName = "BoundaryElement";
      hasInput = true;
      hasModel = true;	  
    } else {
      if (args_info.xyz_given) {
	inFileName = args_info.xyz_arg;
	modelName = "AtomicBead";
	hasInput = true;
	hasModel = true;
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
  
  std::string outputFilename = prefix + ".hydro";

  RealType temperature(0.0);
  RealType viscosity(0.0);

  if (args_info.viscosity_given) {
    viscosity = args_info.viscosity_arg;
    hasViscosity = true;
  }
  if (args_info.temperature_given) {
    temperature = args_info.temperature_arg;
    hasTemperature = true;
  }

  std::ofstream outputHydro;
  outputHydro.open(outputFilename.c_str());  
  HydroIO* hio = new HydroIO();
  hio->openWriter(outputHydro);

  if (args_info.stl_given) {
    std::vector<RealType> coords, normals;
    std::vector<unsigned int> tris, solids;
    Mesh* mesh = new Mesh();

    try {
      stl_reader::ReadStlFile(inFileName.c_str(), coords, normals, tris,
			      solids);
      const size_t numTris = tris.size() / 3;
      for(size_t itri = 0; itri < numTris; ++itri) {
	std::cout << "coordinates of triangle " << itri << ": ";	
	for(size_t icorner = 0; icorner < 3; ++icorner) {
          RealType* c0 = &coords[3 * tris [3 * itri + 0]];
	  Vector3d vert0 (c0[0], c0[1], c0[2]);
          RealType* c1 = &coords[3 * tris [3 * itri + 1]];
	  Vector3d vert1 (c1[0], c1[1], c1[2]);
          RealType* c2 = &coords[3 * tris [3 * itri + 2]];
	  Vector3d vert2 (c2[0], c2[1], c2[2]);
	  mesh->add(vert0, vert1, vert2);	  
	}
      }
    }
    catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
    HydrodynamicsModel* model =
      HydrodynamicsModelFactory::getInstance()->createHydrodynamicsModel(modelName);
    
    if (model != NULL) {
      model->init();
      model->setShape(mesh);

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
  }
  
  
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
      case model_arg_AtomicBeadModel:
      default:
        modelName = "AtomicBeadModel";
        break;
      }

      model =
          HydrodynamicsModelFactory::getInstance()->createHydrodynamicsModel(modelName);
      
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

  // Utils::deletePointers(shapes);
  delete info;
}

void registerHydrodynamicsModels() {
  HydrodynamicsModelFactory::getInstance()->registerHydrodynamicsModel(
      new HydrodynamicsModelBuilder<RoughShell>("RoughShell"));
  HydrodynamicsModelFactory::getInstance()->registerHydrodynamicsModel(
      new HydrodynamicsModelBuilder<AtomicBeadModel>("AtomicBeadModel"));
  HydrodynamicsModelFactory::getInstance()->registerHydrodynamicsModel(
      new HydrodynamicsModelBuilder<AnalyticalModel>("AnalyticalModel"));
}

class XYZFormat {
public:
  XYZFormat() {}
  bool ReadMolecule(OBBase* pOb, OBConversion* pConv) override;
};

bool XYZFormat::ReadMolecule(istream& ifs) {
  char buffer[BUFF_SIZE]; 
  unsigned int natoms;

  if (!ifs)
      return false; // we're attempting to read past the end of the file

  if (!ifs.getline(buffer,BUFF_SIZE)) {
    strcpy(painCave.errMsg, "Problems reading an XYZ file: Cannot read the first line.\n");
    painCave.isFatal = 1;
    simError();
  }
  
  if (sscanf(buffer, "%d", &natoms) == 0 || !natoms) {
    strcpy(painCave.errMsg, "Problems reading an XYZ file: The first line must contain the number of atoms.\n");
    painCave.isFatal = 1;
    simError();
  }

  mol.ReserveAtoms(natoms);

  // The next line contains a title string for the molecule. Use this
  // as the title for the molecule if the line is not
  // empty. Otherwise, use the title given by the calling function.
  if (!ifs.getline(buffer,BUFF_SIZE)) {
    strcpy(painCave.errMsg, "Problems reading an XYZ file: Could not read the second line (title/comments).\n");
    painCave.isFatal = 1;
    simError();
  }
  string readTitle(buffer);
  string::size_type location = readTitle.find("Energy");
  if (location != string::npos)
    readTitle.erase(location);
  Trim(readTitle);
  
  location = readTitle.find_first_not_of(" \t\n\r");
  if (location != string::npos)
    mol.SetTitle(readTitle);
  else
    mol.SetTitle(title);
  
  mol.BeginModify();

  // The next lines contain four items each, separated by white
  // spaces: the atom type, and the coordinates of the atom
  vector<string> vs;
  for (unsigned int i = 1; i <= natoms; i ++) {
    if (!ifs.getline(buffer,BUFF_SIZE)) {
      errorMsg << "Problems reading an XYZ file: "
	       << "Could not read line #" << i+2 << ", file error." << endl
	       << " According to line one, there should be " << natoms
	       << " atoms, and therefore " << natoms+2 << " lines in the file.";
      
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
      return(false);
    }
    tokenize(vs,buffer);
    if (vs.size() < 4) // ignore extra columns which some applications add
      {
	errorMsg << "Problems reading an XYZ file: "
		 << "Could not read line #" << i+2 << "." << endl
		 << "OpenBabel found the line '" << buffer << "'" << endl
		 << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
		 << "However, OpenBabel found " << vs.size() << " items.";
	
	obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
	return(false);
      }

    // Atom Type: get the atomic number from the element table, using
    // the first entry in the currently read line. If the entry makes
    // sense, set the atomic number and leave the atomic type open
    // (the type is then later faulted in when atom->GetType() is
    // called). If the entry does not make sense to use, set the atom
    // type manually, assuming that the author of the xyz-file had
    // something "special" in mind.
    OBAtom *atom  = mol.NewAtom();
    
    int atomicNum = OBElements::GetAtomicNum(vs[0].c_str());
    //set atomic number, or '0' if the atom type is not recognized
    if (atomicNum == 0) {
      // Sometimes people call this an XYZ file, but it's actually Unichem
      // i.e., the first column is the atomic number, not a symbol
      // so we'll first check if we can convert this to an element number
      atomicNum = atoi(vs[0].c_str());
    }
    
    atom->SetAtomicNum(atomicNum);
    if (atomicNum == 0) // still strange, try using an atom type
      atom->SetType(vs[0]);
    
    // Read the atom coordinates
    char *endptr;
    double x = strtod((char*)vs[1].c_str(),&endptr);
    if (endptr == (char*)vs[1].c_str()) {
      
      errorMsg << "Problems reading an XYZ file: "
	       << "Could not read line #" << i+2 << "." << endl
	       << "OpenBabel found the line '" << buffer << "'" << endl
	       << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
	       << "OpenBabel could not interpret item #1 as a number.";
      
      strcpy(painCave.errMsg, errorMsg.c_str());
      painCave.isFatal = 1;
      simError();
    }
    double y = strtod((char*)vs[2].c_str(),&endptr);
    if (endptr == (char*)vs[2].c_str()) {

      errorMsg << "Problems reading an XYZ file: "
	       << "Could not read line #" << i+2 << "." << endl
	       << "OpenBabel found the line '" << buffer << "'" << endl
	       << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
	       << "OpenBabel could not interpret item #2 as a number.";
      
      strcpy(painCave.errMsg, errorMsg.c_str());
      painCave.isFatal = 1;
      simError();
    }
    double z = strtod((char*)vs[3].c_str(),&endptr);
    if (endptr == (char*)vs[3].c_str()) {
      
      errorMsg << "Problems reading an XYZ file: "
	       << "Could not read line #" << i+2 << "." << endl
	       << "OpenBabel found the line '" << buffer << "'" << endl
	       << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
	       << "OpenBabel could not interpret item #3 as a number.";
      
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
      return(false);
    }
    atom->SetVector(x,y,z); //set coordinates
    
    // OK, sometimes there's sym x y z charge -- accepted by Jmol
    if (vs.size() > 5) {
      string::size_type decimal = vs[4].find('.');
      if (decimal !=string::npos) { // period found
	double charge = strtod((char*)vs[4].c_str(),&endptr);
	if (endptr != (char*)vs[4].c_str())
	  atom->SetPartialCharge(charge);
      }
    } // attempt to parse charges
  }
  
  // clean out any remaining blank lines
  std::streampos ipos;
  do {
    ipos = ifs.tellg();
    ifs.getline(buffer,BUFF_SIZE);
  }
  while(strlen(buffer) == 0 && !ifs.eof() );
  ifs.seekg(ipos);
  
  if (!pConv->IsOption("b",OBConversion::INOPTIONS))
    mol.ConnectTheDots();
  if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
    mol.PerceiveBondOrders();
  
  mol.EndModify();
  
  return(true);
}
