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

#include <config.h>

#include <algorithm>
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
#include "nanorod_pentBuilderCmd.hpp"
#include "shapedLatticePentRod.hpp"
#include "shapedLatticeRod.hpp"
#include "utils/Constants.hpp"
#include "utils/MoLocator.hpp"
#include "utils/StringUtils.hpp"

using namespace std;
using namespace OpenMD;
void createMdFile(const std::string& oldMdFileName,
                  const std::string& newMdFileName, std::vector<int> numMol);

int main(int argc, char* argv[]) {
  registerLattice();

  gengetopt_args_info args_info;
  std::string latticeType;
  std::string inputFileName;
  std::string outputFileName;
  MoLocator* locator;
  int nComponents;
  double latticeConstant;
  RealType rodRadius;
  RealType rodLength;
  Mat3x3d hmat;
  DumpWriter* writer;

  // Parse Command Line Arguments
  if (cmdline_parser(argc, argv, &args_info) != 0) exit(1);

  /* get lattice type */
  latticeType = "FCC";

  /* get input file name */
  if (args_info.inputs_num)
    inputFileName = args_info.inputs[0];
  else {
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "No input .omd file name was specified "
             "on the command line");
    painCave.isFatal = 1;
    cmdline_parser_print_help();
    simError();
  }

  /* parse md file and set up the system */
  SimCreator oldCreator;
  SimInfo* oldInfo = oldCreator.createSim(inputFileName, false);

  latticeConstant    = args_info.latticeConstant_arg;
  rodRadius          = args_info.radius_arg;
  rodLength          = args_info.length_arg;
  Globals* simParams = oldInfo->getSimParams();

  /* Create nanorod */
  shapedLatticePentRod nanoRod(latticeConstant, latticeType, rodRadius,
                               rodLength);

  /* Set up the random number generator engine */
  std::random_device rd;   // Non-deterministic, uniformly-distributed integer
                           // random number generator
  std::mt19937 gen(rd());  // 32-bit Mersenne Twister random number engine

  /* Build a lattice and get lattice points for this lattice constant */

  // Rotation angles for lattice
  RealType phi, theta, psi;

  /*
    RealType cphi, sphi, ctheta, stheta, cpsi, spsi;

    cphi = cos(phi);
    sphi = sin(phi);
    ctheta = cos(theta);
    stheta = sin(theta);
    cpsi = cos(psi);
    spsi = sin(psi);
  */

  // Rotates 45 degrees about z-axis
  RotMat3x3d rotation45(45.0 * Constants::PI / 180.0, 0.0, 0.0);

  /*rotation45[0][0] = sqrt(2)/2;
  rotation45[0][1] = -sqrt(2)/2;
  rotation45[0][2] = 0;
  rotation45[1][0] = sqrt(2)/2;
  rotation45[1][1] = sqrt(2)/2;
  rotation45[1][2] = 0;
  rotation45[2][0] = 0;
  rotation45[2][1] = 0;
  rotation45[2][2] = 1;*/

  phi   = 0.0;
  theta = 72.0 * Constants::PI / 180.0;
  psi   = 0.0;

  // Rotates 72 degrees about y-axis
  RotMat3x3d rotation72(phi, theta, psi);

  /*rotation72[0][0] = sqrt(5)/4 - 0.25;
  rotation72[0][1] = 0;
  rotation72[0][2] = sqrt(2*(sqrt(5) + 5))/4;
  rotation72[1][0] = 0;
  rotation72[1][1] = 1;
  rotation72[1][2] = 0;
  rotation72[2][0] = -sqrt(2*(sqrt(5) + 5))/4;
  rotation72[2][1] = 0;
  rotation72[2][2] = sqrt(5)/4 - 0.25;*/

  vector<Vector3d> getsites        = nanoRod.getSites();
  vector<Vector3d> getorientations = nanoRod.getOrientations();
  vector<Vector3d> sites;
  vector<Vector3d> orientations;

  for (unsigned int index = 0; index < getsites.size(); index++) {
    Vector3d mySite   = getsites[index];
    Vector3d myOrient = getorientations[index];
    Vector3d mySite2  = rotation45 * mySite;
    Vector3d o2       = rotation45 * myOrient;
    sites.push_back(mySite2);
    orientations.push_back(o2);

    mySite2 = rotation72 * mySite2;
    o2      = rotation72 * o2;
    sites.push_back(mySite2);
    orientations.push_back(o2);

    mySite2 = rotation72 * mySite2;
    o2      = rotation72 * o2;
    sites.push_back(mySite2);
    orientations.push_back(o2);

    mySite2 = rotation72 * mySite2;
    o2      = rotation72 * o2;
    sites.push_back(mySite2);
    orientations.push_back(o2);

    mySite2 = rotation72 * mySite2;
    o2      = rotation72 * o2;
    sites.push_back(mySite2);
    orientations.push_back(o2);
  }

  int nCenter = int((rodLength + 1.154700538 * rodRadius) / 2.88);

  for (unsigned int index = 0; index <= 0.5 * nCenter; index++) {
    Vector3d myLoc_top(2.88 * index, 0.0, 0.0);
    sites.push_back(myLoc_top);
    orientations.push_back(Vector3d(0.0));
  }

  for (unsigned int index = 1; index <= 0.5 * nCenter; index++) {
    Vector3d myLoc_bottom(-2.88 * index, 0.0, 0.0);
    sites.push_back(myLoc_bottom);
    orientations.push_back(Vector3d(0.0));
  }

  std::vector<std::size_t> vacancyTargets;
  vector<bool> isVacancy;

  Vector3d myLoc;
  RealType myR;

  for (unsigned int i = 0; i < sites.size(); i++)
    isVacancy.push_back(false);

  // cerr << "checking vacancyPercent" << "\n";
  if (args_info.vacancyPercent_given) {
    // cerr << "vacancyPercent given" << "\n";
    if (args_info.vacancyPercent_arg < 0.0 ||
        args_info.vacancyPercent_arg > 100.0) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "vacancyPercent was set to a non-sensical value.");
      painCave.isFatal = 1;
      simError();
    } else {
      RealType vF = args_info.vacancyPercent_arg / 100.0;
      //  cerr << "vacancyPercent = " << vF << "\n";
      RealType vIR;
      RealType vOR;
      if (args_info.vacancyInnerRadius_given) {
        vIR = args_info.vacancyInnerRadius_arg;
      } else {
        vIR = 0.0;
      }
      if (args_info.vacancyOuterRadius_given) {
        vOR = args_info.vacancyOuterRadius_arg;
      } else {
        vOR = rodRadius;
      }
      if (vIR >= 0.0 && vOR <= rodRadius && vOR >= vIR) {
        for (unsigned int i = 0; i < sites.size(); i++) {
          myLoc = sites[i];
          myR   = myLoc.length();
          if (myR >= vIR && myR <= vOR) { vacancyTargets.push_back(i); }
        }
        std::shuffle(vacancyTargets.begin(), vacancyTargets.end(), gen);

        int nTargets = vacancyTargets.size();
        vacancyTargets.resize((int)(vF * nTargets));

        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Removing %d atoms from randomly-selected\n"
                 "\tsites between %lf and %lf.",
                 (int)vacancyTargets.size(), vIR, vOR);
        painCave.severity = OPENMD_INFO;
        painCave.isFatal  = 0;
        simError();

        isVacancy.clear();
        for (std::size_t i = 0; i < sites.size(); i++) {
          bool vac = false;
          for (std::size_t j = 0; j < vacancyTargets.size(); j++) {
            if (i == vacancyTargets[j]) vac = true;
          }
          isVacancy.push_back(vac);
        }

      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Something is strange about the vacancy\n"
                 "\tinner or outer radii.  Check their values.");
        painCave.isFatal = 1;
        simError();
      }
    }
  }

  /* Get number of lattice sites */
  int nSites = sites.size() - vacancyTargets.size();

  // cerr << "sites.size() = " << sites.size() << "\n";
  // cerr << "nSites = " << nSites << "\n";
  // cerr << "vacancyTargets = " << vacancyTargets.size() << "\n";

  std::vector<Component*> components = simParams->getComponents();
  std::vector<RealType> molFractions;
  std::vector<RealType> shellRadii;
  std::vector<int> nMol;
  std::map<int, int> componentFromSite;
  nComponents = components.size();
  // cerr << "nComponents = " << nComponents << "\n";

  if (args_info.molFraction_given && args_info.shellRadius_given) {
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "Specify either molFraction or shellRadius "
             "arguments, but not both!");
    painCave.isFatal = 1;
    simError();
  }

  if (nComponents == 1) {
    molFractions.push_back(1.0);
    shellRadii.push_back(rodRadius);
  } else if (args_info.molFraction_given) {
    if ((int)args_info.molFraction_given == nComponents) {
      for (int i = 0; i < nComponents; i++) {
        molFractions.push_back(args_info.molFraction_arg[i]);
      }
    } else if ((int)args_info.molFraction_given == nComponents - 1) {
      RealType remainingFraction = 1.0;
      for (int i = 0; i < nComponents - 1; i++) {
        molFractions.push_back(args_info.molFraction_arg[i]);
        remainingFraction -= molFractions[i];
      }
      molFractions.push_back(remainingFraction);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "nanorodBuilder can't figure out molFractions "
               "for all of the components in the <MetaData> block.");
      painCave.isFatal = 1;
      simError();
    }
  } else if ((int)args_info.shellRadius_given) {
    if ((int)args_info.shellRadius_given == nComponents) {
      for (int i = 0; i < nComponents; i++) {
        shellRadii.push_back(args_info.shellRadius_arg[i]);
      }
    } else if ((int)args_info.shellRadius_given == nComponents - 1) {
      for (int i = 0; i < nComponents - 1; i++) {
        shellRadii.push_back(args_info.shellRadius_arg[i]);
      }
      shellRadii.push_back(rodRadius);
    } else {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "nanorodBuilder can't figure out the\n"
          "\tshell radii for all of the components in the <MetaData> block.");
      painCave.isFatal = 1;
      simError();
    }
  } else {
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "You have a multi-component <MetaData> block,\n"
             "\tbut have not specified either molFraction or shellRadius "
             "arguments.");
    painCave.isFatal = 1;
    simError();
  }

  if (args_info.molFraction_given) {
    RealType totalFraction = 0.0;

    /* Do some simple sanity checking*/

    for (int i = 0; i < nComponents; i++) {
      if (molFractions.at(i) < 0.0) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "One of the requested molFractions was"
                 " less than zero!");
        painCave.isFatal = 1;
        simError();
      }
      if (molFractions.at(i) > 1.0) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "One of the requested molFractions was"
                 " greater than one!");
        painCave.isFatal = 1;
        simError();
      }
      totalFraction += molFractions.at(i);
    }
    if (abs(totalFraction - 1.0) > 1e-6) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "The sum of molFractions was not close enough to 1.0");
      painCave.isFatal = 1;
      simError();
    }

    int remaining = nSites;
    for (int i = 0; i < nComponents - 1; i++) {
      nMol.push_back(int((RealType)nSites * molFractions.at(i)));
      remaining -= nMol.at(i);
    }
    nMol.push_back(remaining);

    // recompute actual mol fractions and perform final sanity check:

    int totalMolecules = 0;
    for (int i = 0; i < nComponents; i++) {
      molFractions[i] = (RealType)(nMol.at(i)) / (RealType)nSites;
      totalMolecules += nMol.at(i);
    }
    if (totalMolecules != nSites) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Computed total number of molecules is not equal "
               "to the number of lattice sites!");
      painCave.isFatal = 1;
      simError();
    }
  } else {
    for (unsigned int i = 0; i < shellRadii.size(); i++) {
      if (shellRadii.at(i) > rodRadius + 1e-6) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "One of the shellRadius values exceeds the rod Radius.");
        painCave.isFatal = 1;
        simError();
      }
      if (shellRadii.at(i) <= 0.0) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "One of the shellRadius values is smaller than zero!");
        painCave.isFatal = 1;
        simError();
      }
    }
  }

  vector<int> ids;
  if ((int)args_info.molFraction_given) {
    //  cerr << "molFraction given 2" << "\n";
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "Creating a randomized spherically-capped nanorod.");
    painCave.isFatal  = 0;
    painCave.severity = OPENMD_INFO;
    simError();
    /* Random rod is the default case*/

    for (unsigned int i = 0; i < sites.size(); i++)
      if (!isVacancy[i]) ids.push_back(i);

    std::shuffle(ids.begin(), ids.end(), gen);

  } else {
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "Creating an fcc nanorod.");
    painCave.isFatal  = 0;
    painCave.severity = OPENMD_INFO;
    simError();

    // RealType smallestSoFar;
    int myComponent = -1;
    nMol.clear();
    nMol.resize(nComponents);

    // cerr << "shellRadii[0] " << shellRadii[0] << "\n";
    //  cerr << "rodRadius " << rodRadius << "\n";

    for (unsigned int i = 0; i < sites.size(); i++) {
      myLoc = sites[i];
      myR   = myLoc.length();
      // smallestSoFar = rodRadius;
      // cerr << "vac = " << isVacancy[i]<< "\n";

      if (!isVacancy[i]) {
        // for (int j = 0; j < nComponents; j++) {
        //   if (myR <= shellRadii[j]) {
        //     if (shellRadii[j] <= smallestSoFar) {
        //       smallestSoFar = shellRadii[j];
        //       myComponent = j;
        //     }
        //   }
        // }
        myComponent          = 0;
        componentFromSite[i] = myComponent;
        nMol[myComponent]++;
        //	cerr << "nMol for myComp(" << myComponent<<") = " <<
        // nMol[myComponent] <<
        //"\n";
      }
    }
  }
  //     cerr << "nMol = " << nMol.at(0) << "\n";

  outputFileName = args_info.output_arg;

  // creat new .omd file on fly which corrects the number of molecule

  createMdFile(inputFileName, outputFileName, nMol);

  delete oldInfo;

  SimCreator newCreator;
  SimInfo* NewInfo = newCreator.createSim(outputFileName, false);

  // Place molecules
  Molecule* mol;
  SimInfo::MoleculeIterator mi;
  mol = NewInfo->beginMolecule(mi);

  int l = 0;

  for (int i = 0; i < nComponents; i++) {
    locator =
        new MoLocator(NewInfo->getMoleculeStamp(i), NewInfo->getForceField());

    //   cerr << "nMol = " << nMol.at(i) << "\n";
    if (!args_info.molFraction_given) {
      for (unsigned int n = 0; n < sites.size(); n++) {
        if (!isVacancy[n]) {
          if (componentFromSite[n] == i) {
            mol = NewInfo->getMoleculeByGlobalIndex(l);
            locator->placeMol(sites[n], orientations[n], mol);
            l++;
          }
        }
      }
    } else {
      for (int n = 0; n < nMol.at(i); n++) {
        mol = NewInfo->getMoleculeByGlobalIndex(l);
        locator->placeMol(sites[ids[l]], orientations[ids[l]], mol);
        l++;
      }
    }
  }

  // fill Hmat
  hmat(0, 0) = 10.0 * rodRadius;
  hmat(0, 1) = 0.0;
  hmat(0, 2) = 0.0;

  hmat(1, 0) = 0.0;
  hmat(1, 1) = 10.0 * rodRadius;
  hmat(1, 2) = 0.0;

  hmat(2, 0) = 0.0;
  hmat(2, 1) = 0.0;
  hmat(2, 2) = 5.0 * rodLength + 2.0 * rodRadius;

  // set Hmat
  NewInfo->getSnapshotManager()->getCurrentSnapshot()->setHmat(hmat);

  // create dumpwriter and write out the coordinates
  writer = new DumpWriter(NewInfo, outputFileName);

  if (writer == NULL) {
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "Error in creating dumpwriter object ");
    painCave.isFatal = 1;
    simError();
  }

  writer->writeDump();

  // deleting the writer will put the closing at the end of the dump file

  delete writer;

  // cleanup a by calling sim error.....
  snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
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

  unsigned int i = 0;
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
             "\tstatements in component blocks.  Make sure that all\n"
             "\tcomponents in the template file have nMol=1");
    painCave.isFatal = 1;
    simError();
  }
}
