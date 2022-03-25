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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>

#include "brains/Register.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "io/DumpWriter.hpp"
#include "lattice/Lattice.hpp"
#include "lattice/LatticeFactory.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "randomBuilderCmd.hpp"
#include "utils/MoLocator.hpp"
#include "utils/StringUtils.hpp"

using namespace std;
using namespace OpenMD;

void createMdFile(const std::string& oldMdFileName,
                  const std::string& newMdFileName, std::vector<int> nMol);

int main(int argc, char* argv[]) {
  registerLattice();

  gengetopt_args_info args_info;
  std::string latticeType;
  std::string inputFileName;
  std::string outputFileName;
  Lattice* simpleLat;
  RealType latticeConstant;
  std::vector<RealType> lc;
  const RealType rhoConvertConst = 1.661;
  RealType density;
  int nx, ny, nz;
  Mat3x3d hmat;
  MoLocator* locator;
  std::vector<Vector3d> latticePos;
  std::vector<Vector3d> latticeOrt;
  int nMolPerCell;
  DumpWriter* writer;

  // parse command line arguments
  if (cmdline_parser(argc, argv, &args_info) != 0) exit(1);

  density = args_info.density_arg;

  // get lattice type
  latticeType = "FCC";
  if (args_info.lattice_given) { latticeType = args_info.lattice_arg; }

  simpleLat = LatticeFactory::getInstance().createLattice(latticeType);

  if (simpleLat == NULL) {
    sprintf(painCave.errMsg, "Lattice Factory can not create %s lattice\n",
            latticeType.c_str());
    painCave.isFatal = 1;
    simError();
  }
  nMolPerCell = simpleLat->getNumSitesPerCell();

  // get the number of unit cells in each direction:

  nx = args_info.nx_arg;

  if (nx <= 0) {
    sprintf(painCave.errMsg, "The number of unit cells in the x direction "
                             "must be greater than 0.");
    painCave.isFatal = 1;
    simError();
  }

  ny = args_info.ny_arg;

  if (ny <= 0) {
    sprintf(painCave.errMsg, "The number of unit cells in the y direction "
                             "must be greater than 0.");
    painCave.isFatal = 1;
    simError();
  }

  nz = args_info.nz_arg;

  if (nz <= 0) {
    sprintf(painCave.errMsg, "The number of unit cells in the z direction "
                             "must be greater than 0.");
    painCave.isFatal = 1;
    simError();
  }

  int nSites = nMolPerCell * nx * ny * nz;

  // get input file name
  if (args_info.inputs_num)
    inputFileName = args_info.inputs[0];
  else {
    sprintf(painCave.errMsg, "No input .omd file name was specified "
                             "on the command line");
    painCave.isFatal = 1;
    simError();
  }

  // parse md file and set up the system

  SimCreator oldCreator;
  SimInfo* oldInfo   = oldCreator.createSim(inputFileName, false);
  Globals* simParams = oldInfo->getSimParams();

  // Calculate lattice constant (in Angstroms)

  std::vector<Component*> components = simParams->getComponents();
  std::vector<RealType> molFractions;
  std::vector<RealType> molecularMasses;
  std::vector<int> nMol;
  std::size_t nComponents = components.size();

  if (nComponents == 1) {
    molFractions.push_back(1.0);
  } else {
    if (args_info.molFraction_given == nComponents) {
      for (std::size_t i = 0; i < nComponents; i++) {
        molFractions.push_back(args_info.molFraction_arg[i]);
      }
    } else if (args_info.molFraction_given == nComponents - 1) {
      RealType remainingFraction = 1.0;
      for (std::size_t i = 0; i < nComponents - 1; i++) {
        molFractions.push_back(args_info.molFraction_arg[i]);
        remainingFraction -= molFractions[i];
      }
      molFractions.push_back(remainingFraction);
    } else {
      sprintf(painCave.errMsg,
              "randomBuilder can't figure out molFractions "
              "for all of the components in the <MetaData> block.");
      painCave.isFatal = 1;
      simError();
    }
  }

  // do some sanity checking:

  RealType totalFraction = 0.0;

  for (std::size_t i = 0; i < nComponents; i++) {
    if (molFractions.at(i) < 0.0) {
      sprintf(painCave.errMsg, "One of the requested molFractions was"
                               " less than zero!");
      painCave.isFatal = 1;
      simError();
    }
    if (molFractions.at(i) > 1.0) {
      sprintf(painCave.errMsg, "One of the requested molFractions was"
                               " greater than one!");
      painCave.isFatal = 1;
      simError();
    }
    totalFraction += molFractions.at(i);
  }
  if (abs(totalFraction - 1.0) > 1e-6) {
    sprintf(painCave.errMsg,
            "The sum of molFractions was not close enough to 1.0");
    painCave.isFatal = 1;
    simError();
  }

  int remaining = nSites;
  for (std::size_t i = 0; i < nComponents - 1; i++) {
    nMol.push_back(int((RealType)nSites * molFractions.at(i)));
    remaining -= nMol.at(i);
  }
  nMol.push_back(remaining);

  // recompute actual mol fractions and perform final sanity check:

  int totalMolecules = 0;
  RealType totalMass = 0.0;
  for (std::size_t i = 0; i < nComponents; i++) {
    molFractions[i] = (RealType)(nMol.at(i)) / (RealType)nSites;
    totalMolecules += nMol.at(i);
    molecularMasses.push_back(MoLocator::getMolMass(
        oldInfo->getMoleculeStamp(i), oldInfo->getForceField()));
    totalMass += (RealType)(nMol.at(i)) * molecularMasses.at(i);
  }
  RealType avgMass = totalMass / (RealType)totalMolecules;

  if (totalMolecules != nSites) {
    sprintf(painCave.errMsg, "Computed total number of molecules is not equal "
                             "to the number of lattice sites!");
    painCave.isFatal = 1;
    simError();
  }

  latticeConstant = pow(rhoConvertConst * nMolPerCell * avgMass / density,
                        (RealType)(1.0 / 3.0));

  // Set the lattice constant

  lc.push_back(latticeConstant);
  simpleLat->setLatticeConstant(lc);

  // Calculate the lattice sites and fill the lattice vector.

  // Get the standard orientations of the cell sites

  latticeOrt = simpleLat->getLatticePointsOrt();

  vector<Vector3d> sites;
  vector<Vector3d> orientations;

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        // Get the position of the cell sites

        simpleLat->getLatticePointsPos(latticePos, i, j, k);

        for (int l = 0; l < nMolPerCell; l++) {
          sites.push_back(latticePos[l]);
          orientations.push_back(latticeOrt[l]);
        }
      }
    }
  }

  outputFileName = args_info.output_arg;

  // create a new .omd file on the fly which corrects the number of molecules

  createMdFile(inputFileName, outputFileName, nMol);

  delete oldInfo;

  // We need to read in the new SimInfo object, then Parse the
  // md file and set up the system

  SimCreator newCreator;
  SimInfo* newInfo = newCreator.createSim(outputFileName, false);

  // fill Hmat

  hmat(0, 0) = nx * latticeConstant;
  hmat(0, 1) = 0.0;
  hmat(0, 2) = 0.0;

  hmat(1, 0) = 0.0;
  hmat(1, 1) = ny * latticeConstant;
  hmat(1, 2) = 0.0;

  hmat(2, 0) = 0.0;
  hmat(2, 1) = 0.0;
  hmat(2, 2) = nz * latticeConstant;

  // Set Hmat

  newInfo->getSnapshotManager()->getCurrentSnapshot()->setHmat(hmat);

  // place the molecules

  // Randomize a vector of ints:

  vector<int> ids;
  for (std::size_t i = 0; i < sites.size(); i++)
    ids.push_back(i);

  /* Set up the random number generator engine */
  std::random_device rd;   // Non-deterministic, uniformly-distributed integer
                           // random number generator
  std::mt19937 gen(rd());  // 32-bit Mersenne Twister random number engine

  std::shuffle(ids.begin(), ids.end(), gen);

  Molecule* mol;
  int l = 0;
  for (std::size_t i = 0; i < nComponents; i++) {
    locator =
        new MoLocator(newInfo->getMoleculeStamp(i), newInfo->getForceField());
    for (int n = 0; n < nMol.at(i); n++) {
      mol = newInfo->getMoleculeByGlobalIndex(l);
      locator->placeMol(sites[ids[l]], orientations[ids[l]], mol);
      l++;
    }
  }

  // Create DumpWriter and write out the coordinates

  writer = new DumpWriter(newInfo, outputFileName);

  if (writer == NULL) {
    sprintf(painCave.errMsg, "error in creating DumpWriter");
    painCave.isFatal = 1;
    simError();
  }

  writer->writeDump();

  // deleting the writer will put the closing at the end of the dump file.

  delete writer;

  sprintf(painCave.errMsg,
          "A new OpenMD file called \"%s\" has been "
          "generated.\n",
          outputFileName.c_str());
  painCave.isFatal  = 0;
  painCave.severity = OPENMD_INFO;
  simError();
  return 0;
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
        sprintf(buffer, "\tnMol = %i;", nMol.at(i));
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
    sprintf(painCave.errMsg,
            "Couldn't replace the correct number of nMol\n"
            "\tstatements in component blocks.  Make sure that all\n"
            "\tcomponents in the template file have nMol=1");
    painCave.isFatal = 1;
    simError();
  }
}
