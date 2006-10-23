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

#include "nanorodBuilderCmd.h"
#ifdef HAVE_CGAL
#include "GeometryBuilder.hpp"
#endif
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
  std::string outputFileName;
  std::string outGeomFileName;
  
  
  Lattice *simpleLat;
  int numMol;
  RealType latticeConstant;
  std::vector<RealType> lc;
  RealType mass;
  const RealType rhoConvertConst = 1.661;
  RealType density;
  RealType rodLength;
  RealType rodDiameter;
  
  
  int nx, ny, nz;
  Mat3x3d hmat;
  MoLocator *locator;
  std::vector<Vector3d> latticePos;
  std::vector<Vector3d> latticeOrt;
  int numMolPerCell;
  int nComponents;
  DumpWriter *writer;
  
  // parse command line arguments
  if (cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);
  
  
  // Check for lib CGAL, if we don't have it, we should exit....
	
#ifndef HAVE_CGAL
  std::cerr << "nanoRodBuilder requires libCGAL to function, please rebuild OOPSE with libCGAL."
	    << std::endl;
  exit(1);
#endif
	
	
	
  //get lattice type
  latticeType = UpperCase(args_info.latticetype_arg);
    
  /* get input file name */
  if (args_info.inputs_num)
    inputFileName = args_info.inputs[0];
  else {
    sprintf(painCave.errMsg, "No input .md file name was specified "
            "on the command line");
    painCave.isFatal = 1;
    cmdline_parser_print_help();
    simError();
  }

  //parse md file and set up the system
  SimCreator oldCreator;
  SimInfo* oldInfo = oldCreator.createSim(inputFileName, false);
  Globals* simParams = oldInfo->getSimParams();

  nComponents = simParams->getNComponents();
  if (nComponents> 1) {
    sprintf(painCave.errMsg, "Nanorods can only contain a single component ");
    painCave.isFatal = 1;
    simError();
  }
  
  //get mass of molecule. 
  
  mass = getMolMass(oldInfo->getMoleculeStamp(0), oldInfo->getForceField());
  
  //creat lattice
  simpleLat = LatticeFactory::getInstance()->createLattice(latticeType);
  
  if (simpleLat == NULL) {
    sprintf(painCave.errMsg, "Error in creating lattice. ");
    painCave.isFatal = 1;
    simError();
  }
  
  numMolPerCell = simpleLat->getNumSitesPerCell();
  
  //calculate lattice constant (in Angstrom)
  //latticeConstant = pow(rhoConvertConst * numMolPerCell * mass / density,
  //                      1.0 / 3.0);
  
  latticeConstant = args_info.latticeCnst_arg;
  rodLength = args_info.length_arg;
  rodDiameter = args_info.width_arg;
  
  //set lattice constant
  lc.push_back(latticeConstant);
  simpleLat->setLatticeConstant(lc);
  
  
  //determine the output file names  
  if (args_info.output_given){
    outputFileName = args_info.output_arg;
  }else{
    sprintf(painCave.errMsg, "No output file name was specified "
            "on the command line");
    painCave.isFatal = 1;
    cmdline_parser_print_help();
    simError();
  }
  
  
         
  
  
  
  //creat Molocator
  locator = new MoLocator(oldInfo->getMoleculeStamp(0), oldInfo->getForceField());
  
  /*
    Assume we are carving nanorod out of a cublic block of material and that
    the shape the material will fit within that block.... 
    The model in geometry builder assumes the long axis is in the y direction and the x-z plane is the 
    diameter of the particle.		 
  */
  // Number of Unit Cells in Length first
  ny = (int)(rodLength/latticeConstant);
  // Number of unit cells in Width
  nx = (int)(rodDiameter/latticeConstant);
  nz = (int)(rodDiameter/latticeConstant);
  
  
  
  // Create geometry for nanocrystal
#ifdef HAVE_CGAL
  GeometryBuilder *myGeometry;
  // GeometryBuilder myGeometry(rodLength,rodDiameter);
  if (args_info.twinnedCrystal_flag){
     myGeometry = new GeometryBuilder(rodLength,rodDiameter);
  }
  else{
     myGeometry = new GeometryBuilder(rodLength,rodDiameter);
  }
  
  if (args_info.genGeomview_given){
     if (args_info.genGeomview_flag){
        outGeomFileName = getPrefix(inputFileName.c_str()) + ".off";
        myGeometry->dumpGeometry(outGeomFileName);
     }
  }
  
#endif
  
  /*
    We have to build the system first to figure out how many molecules
    there are then create a md file and then actually build the
    system.
  */
  
  //place the molecules
  

  
  //get the orientation of the cell sites
  //for the same type of molecule in same lattice, it will not change
  latticeOrt = simpleLat->getLatticePointsOrt();
  
  
  
  numMol = 0;
  for(int i = -nx; i < nx; i++) {     
    for(int j = -ny; j < ny; j++) {       
      for(int k = -nz; k < nz; k++) {
        //if (oldInfo->getNGlobalMolecules() != numMol) {
	
	
	
        //get the position of the cell sites
        simpleLat->getLatticePointsPos(latticePos, i, j, k);
	
        for(int l = 0; l < numMolPerCell; l++) {

#ifdef HAVE_CGAL
          if (myGeometry->isInsidePolyhedron(latticePos[l][0],latticePos[l][1],latticePos[l][2])){
            numMol++;
          }
#else
           numMol++;
#endif
        }
      }
    }
  }

  
  // needed for writing out new md file.
  


  
  //creat new .md file on fly which corrects the number of molecule     
  createMdFile(inputFileName, outputFileName, numMol);
  
  if (oldInfo != NULL)
    delete oldInfo;
  
  
  // We need to read in new siminfo object.	
  //parse md file and set up the system
  SimCreator newCreator;
  SimInfo* newInfo = newCreator.createSim(outputFileName, false);
  
  // This was so much fun the first time, lets do it again.
  
  Molecule* mol;
  SimInfo::MoleculeIterator mi;
  mol = newInfo->beginMolecule(mi);


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
		   sprintf(painCave.errMsg, "Error in placing molecule onto lattice ");
		   painCave.isFatal = 1;
		   simError();
                 }
                 mol = newInfo->nextMolecule(mi);
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
  newInfo->getSnapshotManager()->getCurrentSnapshot()->setHmat(hmat);
  
  
  //create dumpwriter and write out the coordinates
  newInfo->setFinalConfigFileName(outputFileName);
  writer = new DumpWriter(newInfo);
  
  if (writer == NULL) {
    sprintf(painCave.errMsg, "Error in creating DumpWrite object. ");
    painCave.isFatal = 1;
    simError();
  }
  
  writer->writeEor();
  std::cout << "new initial configuration file: " << outputFileName
            << " is generated." << std::endl;
  
  //delete objects
  
  //delete oldInfo and oldSimSetup
  
  if (writer != NULL)
    delete writer;
  delete simpleLat;	
  cmdline_parser_free(&args_info);
  return 0;
}

void createMdFile(const std::string&oldMdFileName, const std::string&newMdFileName,
                  int numMol) {
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

