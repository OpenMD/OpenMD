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

#include "config.h"

#include "nanoparticleBuilderCmd.h"
#include "sphericalNanoparticle.hpp"
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
                  int numMol);

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
  int numMol;
  double latticeConstant;
  std::vector<double> lc;
  double mass;
  const double rhoConvertConst = 1.661;
  double density;
  
  

  Mat3x3d hmat;
  MoLocator *locator;
  sphericalNanoparticle *nanoparticle;
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
  
  nShells = 0;
  if (args_info.coreShellRadius_given){
    nShells = args_info.coreShellRadius_given;
  }
  
  nComponents = oldInfo->getNMoleculeStamp();
  
  /* Check to see if we have enough components to build that many shells. */
  if (nShells){
    if (oldInfo->getNMoleculeStamp() != nShells) {
      std::cerr << "Not enough components present in MD file to build specified number of shells"
      << std::endl;
      exit(1);
    }
  }
  
  
  //creat lattice
  simpleLat = LatticeFactory::getInstance()->createLattice(latticeType);
  
  if (simpleLat == NULL) {
    std::cerr << "Error in creating lattice" << std::endl;
    exit(1);
  }
  
  numMolPerCell = simpleLat->getNumSitesPerCell();
  
  /*calculate lattice constant (in Angstrom)
  latticeConstant = pow(rhoConvertConst * numMolPerCell * mass / density,
                        1.0 / 3.0);*/
  
  latticeConstant = args_info.latticeCnst_arg;
  particleRadius = args_info.radius_arg;
  particleDiameter = 2.0 * particleRadius;
  
  /* set lattice constant */
  lc.push_back(latticeConstant);
  simpleLat->setLatticeConstant(lc);
  
  
  /*determine the output file names*/
  if (args_info.output_given){
    outInitFileName = args_info.output_arg;
  }else{
    outInitFileName = getPrefix(inputFileName.c_str()) + ".in";
  }
  
         
  
  
  
  /* create Molocators */
  locator = new MoLocator(oldInfo->getMoleculeStamp(0), oldInfo->getForceField());
  
  /* create a new spherical nanoparticle */
  nanoparticle = new sphericalNanoparticle(particleRadius,latticeConstant);
  /* Build a nanoparticle to see how many sites are there */
  numSites = new int[nComponents]
  int numSites = nanoparticle.getNMol(simpleLat);
  

  /* Random particle is the default case*/
  if (!args_info.ShellRadius_given){
    std::cout << "Creating a random nanoparticle" << std::endl;
    /* Check to see if we have enough components */
    if (nComponents != args_info.molFraction_given + 1){
      std::cerr << "Number of components does not equal molFraction occurances." << std::endl;
      exit 1;
    }
    int totComponents = 0;
    for (int i = 0;i<nComponents-2;i++){ /* Figure out Percent for each component */
      numMol[i] = int((double)numSites * args_info.molFraction_arg[i]);
      totComponents += numMol[i];
    }
    numMol[nComponents-1] = numSites - totComponents;

  } else{ /*Handle core-shell with multiple components.*/
    std::cout << "Creating a core-shell nanoparticle." << std::endl;
    if (nComponents != args_info.ShellRadius_given + 1){
      std::cerr << "Number of components does not equal ShellRadius occurances." << std::endl;
      exit 1;
    }
    
    
    
  }

   //get the orientation of the cell sites
  //for the same type of molecule in same lattice, it will not change
   latticeOrt = simpleLat->getLatticePointsOrt();
  
  
  
  // needed for writing out new md file.
  
    outPrefix = getPrefix(inputFileName.c_str()) + "_" + latticeType;
    outMdFileName = outPrefix + ".md";
  
    //creat new .md file on fly which corrects the number of molecule     
    createMdFile(inputFileName, outMdFileName, numcomponents,numMol);
  
  if (oldInfo != NULL)
    delete oldInfo;
  
  
  // We need to read in new siminfo object.	
  //parse md file and set up the system
  //SimCreator NewCreator;
  
  SimInfo* NewInfo = oldCreator.createSim(outMdFileName, false);
  
  // This was so much fun the first time, lets do it again.
  
  Molecule* mol;
  SimInfo::MoleculeIterator mi;
  mol = NewInfo->beginMolecule(mi);


  for(int i = -nx; i < nx; i++) {
     for(int j = -ny; j < ny; j++) {
        for(int k = -nz; k < nz; k++) {
           
           //get the position of the cell sites
           simpleLat->getLatticePointsPos(latticePos, i, j, k);
           
           for(int l = 0; l < numMolPerCell; l++) {
#ifdef HAVE_CGAL              
              if (myGeometry->isInsidePolyhedron(latticePos[l][0],latticePos[l][1],latticePos[l][2])){
#endif                              
                 if (mol != NULL) {
                    locator->placeMol(latticePos[l], latticeOrt[l], mol);
                 } else {
                    std::cerr<<"Error in placing molecule " << std::endl;                    
                 }
                 mol = NewInfo->nextMolecule(mi);
#ifdef HAVE_CGAL                
              }
#endif              
           }
        }
     }
  }
  

  
  //fill Hmat
  hmat(0, 0)= nx * latticeConstant;
  hmat(0, 1) = 0.0;
  hmat(0, 2) = 0.0;
  
  hmat(1, 0) = 0.0;
  hmat(1, 1) = ny * latticeConstant;
  hmat(1, 2) = 0.0;
  
  hmat(2, 0) = 0.0;
  hmat(2, 1) = 0.0;
  hmat(2, 2) = nz * latticeConstant;
  
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
  delete simpleLat;	
  cmdline_parser_free(&args_info);
  return 0;
}

void createMdFile(const std::string&oldMdFileName, const std::string&newMdFileName,
                  int components,int &nummol) {
  ifstream oldMdFile;
  ofstream newMdFile;
  const int MAXLEN = 65535;
  char buffer[MAXLEN];
  
  //create new .md file based on old .md file
  oldMdFile.open(oldMdFileName.c_str());
  newMdFile.open(newMdFileName.c_str());
  
  oldMdFile.getline(buffer, MAXLEN);
  
  while (!oldMdFile.eof()) {
    
    //correct molecule number
    if (strstr(buffer, "nMol") != NULL) {
      sprintf(buffer, "\tnMol = %i;", numMol);				
      newMdFile << buffer << std::endl;
    } else
      newMdFile << buffer << std::endl;
    
    oldMdFile.getline(buffer, MAXLEN);
  }
  
  oldMdFile.close();
  newMdFile.close();
}

