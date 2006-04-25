/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 */

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <algorithm>

#include "config.h"
#include "shapedLatticeSpherical.hpp"
#include "nanoparticleBuilderCmd.h"
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
                  int components,int* numMol);

int main(int argc, char *argv []) {
  
  //register force fields
  registerForceFields();
  registerLattice();
  
  gengetopt_args_info args_info;
  std::string latticeType;
  std::string inputFileName;
  std::string outPrefix;
  std::string outMdFileName;
  std::string outInitFileName;

  
  
  Lattice *simpleLat;
  MoLocator* locator;
  int* numMol;
  int nComponents;
  double latticeConstant;
  std::vector<double> lc;
  double mass;											    
  const double rhoConvertConst = 1.661;
  double density;
  double particleRadius;
  
  

  Mat3x3d hmat;
  std::vector<Vector3d> latticePos;
  std::vector<Vector3d> latticeOrt;
  int numMolPerCell;
  int nShells; /* Number of shells in nanoparticle*/
  int numSites;
  
  DumpWriter *writer;
  
  // Parse Command Line Arguments
  if (cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);
  
	
	
  /* get lattice type */
  latticeType = UpperCase(args_info.latticetype_arg);
    
  /* get input file name */
  if (args_info.inputs_num)
    inputFileName = args_info.inputs[0];
  else {
    std::cerr << "You must specify a input file name.\n" << std::endl;
    cmdline_parser_print_help();
    exit(1);
  }
  
  /* parse md file and set up the system */
  SimCreator oldCreator;
  SimInfo* oldInfo = oldCreator.createSim(inputFileName, false);
  
  
  /*calculate lattice constant (in Angstrom)
    latticeConstant = pow(rhoConvertConst * numMolPerCell * mass / density,
    1.0 / 3.0);*/
  
  latticeConstant = args_info.latticeCnst_arg;
  particleRadius = args_info.radius_arg;
  Globals* simParams = oldInfo->getSimParams();
  
  /* Find out how many different components in this simualtion */
  nComponents =simParams->getNComponents();
  
  /*determine the output file names*/
  if (args_info.output_given){
    outInitFileName = args_info.output_arg;
  }else{
    outInitFileName = getPrefix(inputFileName.c_str()) + ".in";
  }
   
  std::cout <<"Before build shaped lattice. "<<std::endl;
  
  /* create Molocators */
  locator = new MoLocator(oldInfo->getMoleculeStamp(0), oldInfo->getForceField());
  
  /* Create nanoparticle */
  shapedLatticeSpherical nanoParticle(latticeConstant,latticeType,particleRadius);
  
  std::cout <<"Before build getPoints. "<<std::endl;
  /* Build a lattice and get lattice points for this lattice constant */
  vector<Vector3d> nanoParticleSites = nanoParticle.getPoints();
  
  /* Get number of lattice sites */
  numSites = nanoParticleSites.size();
 
  //std::cout <<"numSites are %d "<<numSites<<std::endl;
  // std::cout <<"nComponents are %d "<<nComponents<<std::endl;
  numMol = new int[nComponents];
 
  
  /* Random particle is the default case*/
  if (!args_info.ShellRadius_given){
    std::cout << "Creating a random nanoparticle" << std::endl;
    /* Check to see if we have enough components */
    if (nComponents != args_info.molFraction_given && nComponents != 1){
      std::cerr << "Number of components does not equal molFraction occurances." << std::endl;
      exit(1);
    }
    /* Build the mole fractions and number of molecules of each type */   
    int totComponents = 0;
    for (int i = 0;i<nComponents-1;i++){ /* Figure out Percent for each component */
      numMol[i] = int((double)numSites * args_info.molFraction_arg[i]);
      std::cout<<numMol[i]<<std::endl;
      totComponents += numMol[i];
    }
    numMol[nComponents-1] = numSites - totComponents;
   
    /* do the iPod thing, Shuffle da vector */
    std::random_shuffle(nanoParticleSites.begin(), nanoParticleSites.end());
  } else{ /*Handle core-shell with multiple components.*/
    std::cout << "Creating a core-shell nanoparticle." << std::endl;
    if (nComponents != args_info.ShellRadius_given + 1){
      std::cerr << "Number of components does not equal ShellRadius occurances." << std::endl;
      exit(1);
    }   
    
  }

  
  //get the orientation of the cell sites
  //for the same type of molecule in same lattice, it will not change
  latticeOrt = nanoParticle.getPointsOrt();
  std::cout<<"Orientational vector Size: "<< std::endl;
  std::cout<<latticeOrt.size()<< std::endl;
  
  
  
  // needed for writing out new md file.
  
  outPrefix = getPrefix(inputFileName.c_str()) + "_" + latticeType;
  outMdFileName = outPrefix + ".md";
  
  //creat new .md file on fly which corrects the number of molecule     
  createMdFile(inputFileName, outMdFileName, nComponents,numMol);
  
  if (oldInfo != NULL)
    delete oldInfo;
  
  
  // We need to read in new siminfo object.	
  //parse md file and set up the system
  //SimCreator NewCreator;
  SimCreator newCreator;
  SimInfo* NewInfo = newCreator.createSim(outMdFileName, false);
  
  
  // Place molecules
  Molecule* mol;
  SimInfo::MoleculeIterator mi;
  mol = NewInfo->beginMolecule(mi);
  int l = 0;
  for (mol = NewInfo->beginMolecule(mi); mol != NULL; mol = NewInfo->nextMolecule(mi)) {
    locator->placeMol(nanoParticleSites[l], latticeOrt[l], mol);
    l++;
  }
 


  
  //fill Hmat
  hmat(0, 0)=  latticeConstant;
  hmat(0, 1) = 0.0;
  hmat(0, 2) = 0.0;
  
  hmat(1, 0) = 0.0;
  hmat(1, 1) =  latticeConstant;
  hmat(1, 2) = 0.0;
  
  hmat(2, 0) = 0.0;
  hmat(2, 1) = 0.0;
  hmat(2, 2) =  latticeConstant;
  
  //set Hmat
  NewInfo->getSnapshotManager()->getCurrentSnapshot()->setHmat(hmat);
  
  
  //create dumpwriter and write out the coordinates
  NewInfo->setFinalConfigFileName(outInitFileName);
  writer = new DumpWriter(NewInfo);
  
  if (writer == NULL) {
    std::cerr << "error in creating DumpWriter" << std::endl;
    exit(1);
  }
  
  writer->writeDump();
  std::cout << "new initial configuration file: " << outInitFileName
            << " is generated." << std::endl;
  
  //delete objects
  
  //delete oldInfo and oldSimSetup
  
  if (NewInfo != NULL)
    delete NewInfo;
  
  if (writer != NULL)
    delete writer;	
  cmdline_parser_free(&args_info);
  return 0;
}

void createMdFile(const std::string&oldMdFileName, const std::string&newMdFileName,
                  int components,int* numMol) {
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

