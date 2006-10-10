/* Copyright (c) 2006 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 *
 *  randomBuilder.cpp
 *
 *  Created by Charles F. Vardeman II on 10 Apr 2006.
 *  @author  Charles F. Vardeman II
 *  @version $Id: randomBuilder.cpp,v 1.4 2006-10-10 02:44:13 gezelter Exp $
 *
 */


#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <fstream>

#include "applications/randomBuilder/randomBuilderCmd.h"
#include "lattice/LatticeFactory.hpp"
#include "utils/MoLocator.hpp"
#include "lattice/Lattice.hpp"
#include "brains/Register.hpp"
#include "brains/SimInfo.hpp"
#include "brains/SimCreator.hpp"
#include "io/DumpWriter.hpp"
#include "math/Vector3.hpp"
#include "math/SquareMatrix3.hpp"
#include "utils/StringUtils.hpp"

using namespace std;
using namespace oopse;

void createMdFile(const std::string&oldMdFileName, 
                  const std::string&newMdFileName,
                  int components, int* numMol);

int main(int argc, char *argv []) {

  // register force fields
  registerForceFields();
  registerLattice();
    
  gengetopt_args_info args_info;
  std::string latticeType;
  std::string inputFileName;
  std::string outputFileName;
  Lattice *simpleLat;
  int* numMol;
  RealType latticeConstant;
  std::vector<RealType> lc;
  RealType mass;
  const RealType rhoConvertConst = 1.661;
  RealType density;
  int nx, ny, nz;
  Mat3x3d hmat;
  MoLocator *locator;
  std::vector<Vector3d> latticePos;
  std::vector<Vector3d> latticeOrt;
  int numMolPerCell;
  int curMolIndex;
  DumpWriter *writer;

  // parse command line arguments
  if (cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  density = args_info.density_arg;

  //get lattice type
  latticeType = UpperCase(args_info.latticetype_arg);

  simpleLat = LatticeFactory::getInstance()->createLattice(latticeType);
    
  if (simpleLat == NULL) {
    sprintf(painCave.errMsg, "Lattice Factory can not create %s lattice\n",
	    latticeType.c_str());
    painCave.isFatal = 1;
    simError();
  }

  //get the number of unit cells in each direction:

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

  //get input file name
  if (args_info.inputs_num)
    inputFileName = args_info.inputs[0];
  else {
    sprintf(painCave.errMsg, "No input .md file name was specified "
            "on the command line");
    painCave.isFatal = 1;
    simError();
  }

  //parse md file and set up the system

  SimCreator oldCreator;
  SimInfo* oldInfo = oldCreator.createSim(inputFileName, false);
  Globals* simParams = oldInfo->getSimParams();

  int nComponents =simParams->getNComponents();
  if (oldInfo->getNMoleculeStamp() > 2) {
    sprintf(painCave.errMsg, "randomBuilder can't yet build a system with "
            "more than two components.");
    painCave.isFatal = 1;
    simError();
  }

  //get mass of molecule. 

  mass = getMolMass(oldInfo->getMoleculeStamp(0), oldInfo->getForceField());

  // Create the lattice

  simpleLat = LatticeFactory::getInstance()->createLattice(latticeType);

  if (simpleLat == NULL) {
    sprintf(painCave.errMsg, "Error in creating lattice.");
    painCave.isFatal = 1;
    simError();
  }

  numMolPerCell = simpleLat->getNumSitesPerCell();

  // Calculate lattice constant (in Angstroms)

  latticeConstant = pow(rhoConvertConst * numMolPerCell * mass / density,
			(RealType)(1.0 / 3.0));

  // Set the lattice constant

  lc.push_back(latticeConstant);
  simpleLat->setLatticeConstant(lc);

  // Calculate the total number of molecules

  int totMol = nx * ny * nz * numMolPerCell;

  // Calculate the lattice sites and fill the lattice vector.

  // Get the standard orientations of the cell sites

  latticeOrt = simpleLat->getLatticePointsOrt();

  vector<Vector3d> sites;
  vector<Vector3d> orientations;
  
  for(int i = 0; i < nx; i++) {
    for(int j = 0; j < ny; j++) {
      for(int k = 0; k < nz; k++) {

	// Get the position of the cell sites

	simpleLat->getLatticePointsPos(latticePos, i, j, k);

	for(int l = 0; l < numMolPerCell; l++) {
	  sites.push_back(latticePos[l]);
          orientations.push_back(latticeOrt[l]);
	}
      }
    }
  }

  int numSites = sites.size();

  numMol = new int[nComponents];
  if (nComponents != args_info.molFraction_given && nComponents != 1){
    sprintf(painCave.errMsg, "There needs to be the same number of "
            "molFraction arguments as there are components in the "
            "<MetaData> block.");
    painCave.isFatal = 1;
    simError();
  }
  int totComponents = 0;
  for (int i = 0;i<nComponents-1;i++){
    numMol[i] = int((RealType)numSites * args_info.molFraction_arg[i]);
    std::cout<<numMol[i]<<std::endl;
    totComponents += numMol[i];
  }
  numMol[nComponents-1] = numSites - totComponents;
  
  outputFileName = args_info.output_arg;

  // create a new .md file on the fly which corrects the number of molecules

  createMdFile(inputFileName, outputFileName, nComponents, numMol);

  if (oldInfo != NULL)
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

  curMolIndex = 0;

  // Randomize a vector of ints:

  vector<int> ids;
  for (int i = 0; i < sites.size(); i++) ids.push_back(i);
  std::random_shuffle(ids.begin(), ids.end());

  Molecule* mol;
  int l = 0;
  for (int i = 0; i < nComponents; i++){
    locator = new MoLocator(newInfo->getMoleculeStamp(i), 
                            newInfo->getForceField());
    for (int n = 0; n < numMol[i]; n++) {
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

  sprintf(painCave.errMsg, "A new OOPSE MD file called \"%s\" has been "
          "generated.", outputFileName.c_str());
  painCave.isFatal = 0;
  simError();
  return 0;
}

void createMdFile(const std::string&oldMdFileName, 
                  const std::string&newMdFileName, 
                  int components, int* numMol) {
  ifstream oldMdFile;
  ofstream newMdFile;
  const int MAXLEN = 65535;
  char buffer[MAXLEN];
  
  //create new .md file based on old .md file

  oldMdFile.open(oldMdFileName.c_str());
  newMdFile.open(newMdFileName.c_str());
  
  oldMdFile.getline(buffer, MAXLEN);
 
  int i = 0;
  while (!oldMdFile.eof()) {
    
    //correct molecule number
    if (strstr(buffer, "nMol") != NULL) {
      if(i<components){
	sprintf(buffer, "\tnMol = %i;", numMol[i]);
	newMdFile << buffer << std::endl;
	i++;
      }
    } else
      newMdFile << buffer << std::endl;
    
    oldMdFile.getline(buffer, MAXLEN);
  }
  
  oldMdFile.close();
  newMdFile.close();
}

