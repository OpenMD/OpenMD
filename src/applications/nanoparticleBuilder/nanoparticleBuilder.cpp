/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
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
using namespace OpenMD;
void createMdFile(const std::string&oldMdFileName, 
                  const std::string&newMdFileName,
                  std::vector<int> numMol);

int main(int argc, char *argv []) {
  
  //register force fields
  registerForceFields();
  registerLattice();
  
  gengetopt_args_info args_info;
  std::string latticeType;
  std::string inputFileName;
  std::string outputFileName;

  MoLocator* locator;
  int nComponents;
  double latticeConstant;
  std::vector<double> lc;

  RealType particleRadius;

  Mat3x3d hmat;
  std::vector<Vector3d> latticePos;
  std::vector<Vector3d> latticeOrt;
 
  DumpWriter *writer;
  
  // Parse Command Line Arguments
  if (cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);
         
  /* get lattice type */
  latticeType = "FCC";

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
  
  /* parse md file and set up the system */
  SimCreator oldCreator;
  SimInfo* oldInfo = oldCreator.createSim(inputFileName, false);
  
  latticeConstant = args_info.latticeConstant_arg;
  particleRadius = args_info.radius_arg;
  Globals* simParams = oldInfo->getSimParams();
  
  /* Create nanoparticle */
  shapedLatticeSpherical nanoParticle(latticeConstant, latticeType,
                                      particleRadius);
  
  /* Build a lattice and get lattice points for this lattice constant */
  vector<Vector3d> sites = nanoParticle.getSites();
  vector<Vector3d> orientations = nanoParticle.getOrientations();
  std::vector<int> vacancyTargets;
  vector<bool> isVacancy;
  
  Vector3d myLoc;
  RealType myR;
 
  for (int i = 0; i < sites.size(); i++) 
    isVacancy.push_back(false);

  if (args_info.vacancyPercent_given) {
    if (args_info.vacancyPercent_arg < 0.0 || args_info.vacancyPercent_arg > 100.0) {
      sprintf(painCave.errMsg, "vacancyPercent was set to a non-sensical value.");
      painCave.isFatal = 1;
      simError();
    } else {
      RealType vF = args_info.vacancyPercent_arg / 100.0;
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
        vOR = particleRadius;
      }
      if (vIR >= 0.0 && vOR <= particleRadius && vOR >= vIR) {
        
        for (int i = 0; i < sites.size(); i++) {
          myLoc = sites[i];
          myR = myLoc.length();
          if (myR >= vIR && myR <= vOR) {
            vacancyTargets.push_back(i);
          }          
        }
        std::random_shuffle(vacancyTargets.begin(), vacancyTargets.end());
        
        int nTargets = vacancyTargets.size();
        vacancyTargets.resize((int)(vF * nTargets));
        
                  
        sprintf(painCave.errMsg, "Removing %d atoms from randomly-selected\n"
                "\tsites between %lf and %lf.", (int) vacancyTargets.size(), 
                vIR, vOR); 
        painCave.isFatal = 0;
        simError();

        isVacancy.clear();
        for (int i = 0; i < sites.size(); i++) {
          bool vac = false;
          for (int j = 0; j < vacancyTargets.size(); j++) {
            if (i == vacancyTargets[j]) vac = true;
          }
          isVacancy.push_back(vac);
        }
               
      } else {
        sprintf(painCave.errMsg, "Something is strange about the vacancy\n"
                "\tinner or outer radii.  Check their values.");
        painCave.isFatal = 1;
        simError();
      }
    }
  }

  /* Get number of lattice sites */
  int nSites = sites.size() - vacancyTargets.size();

  std::vector<Component*> components = simParams->getComponents();
  std::vector<RealType> molFractions;
  std::vector<RealType> shellRadii;
  std::vector<RealType> molecularMasses;
  std::vector<int> nMol;
  std::map<int, int> componentFromSite;
  nComponents = components.size();

  if (args_info.molFraction_given && args_info.shellRadius_given) {
    sprintf(painCave.errMsg, "Specify either molFraction or shellRadius "
            "arguments, but not both!");
    painCave.isFatal = 1;
    simError();
  }
  
  if (nComponents == 1) {
    molFractions.push_back(1.0);    
    shellRadii.push_back(particleRadius);
  } else if (args_info.molFraction_given) {
    if ((int)args_info.molFraction_given == nComponents) {
      for (int i = 0; i < nComponents; i++) {
        molFractions.push_back(args_info.molFraction_arg[i]);
      }
    } else if ((int)args_info.molFraction_given == nComponents-1) {
      RealType remainingFraction = 1.0;
      for (int i = 0; i < nComponents-1; i++) {
        molFractions.push_back(args_info.molFraction_arg[i]);
        remainingFraction -= molFractions[i];
      }
      molFractions.push_back(remainingFraction);
    } else {    
      sprintf(painCave.errMsg, "nanoparticleBuilder can't figure out molFractions "
              "for all of the components in the <MetaData> block.");
      painCave.isFatal = 1;
      simError();
    }
  } else if ((int)args_info.shellRadius_given) {
    if ((int)args_info.shellRadius_given == nComponents) {
      for (int i = 0; i < nComponents; i++) {
        shellRadii.push_back(args_info.shellRadius_arg[i]);
      }
    } else if ((int)args_info.shellRadius_given == nComponents-1) {
      for (int i = 0; i < nComponents-1; i++) {
        shellRadii.push_back(args_info.shellRadius_arg[i]);
      }
      shellRadii.push_back(particleRadius);
    } else {    
      sprintf(painCave.errMsg, "nanoparticleBuilder can't figure out the\n"
              "\tshell radii for all of the components in the <MetaData> block.");
      painCave.isFatal = 1;
      simError();
    }
  } else {
    sprintf(painCave.errMsg, "You have a multi-component <MetaData> block,\n"
            "\tbut have not specified either molFraction or shellRadius arguments.");
    painCave.isFatal = 1;
    simError();
  }
    
  if (args_info.molFraction_given) {
    RealType totalFraction = 0.0;
    
    /* Do some simple sanity checking*/
    
    for (int i = 0; i < nComponents; i++) {
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
      sprintf(painCave.errMsg, "The sum of molFractions was not close enough to 1.0");
      painCave.isFatal = 1;
      simError();
    }
    
    int remaining = nSites;
    for (int i=0; i < nComponents-1; i++) {    
      nMol.push_back(int((RealType)nSites * molFractions.at(i)));
      remaining -= nMol.at(i);
    }
    nMol.push_back(remaining);
    
    // recompute actual mol fractions and perform final sanity check:
    
    int totalMolecules = 0;
    for (int i=0; i < nComponents; i++) {
      molFractions[i] = (RealType)(nMol.at(i))/(RealType)nSites;
      totalMolecules += nMol.at(i);
    }
    
    if (totalMolecules != nSites) {
      sprintf(painCave.errMsg, "Computed total number of molecules is not equal "
              "to the number of lattice sites!");
      painCave.isFatal = 1;
      simError();
    }
  } else {

    for (int i = 0; i < shellRadii.size(); i++) {
      if (shellRadii.at(i) > particleRadius + 1e-6 ) {
        sprintf(painCave.errMsg, "One of the shellRadius values exceeds the particle Radius.");
        painCave.isFatal = 1;
        simError();
      } 
      if (shellRadii.at(i) <= 0.0 ) {
        sprintf(painCave.errMsg, "One of the shellRadius values is smaller than zero!");
        painCave.isFatal = 1;
        simError();
      }
    }
  }

  vector<int> ids;           
  if ((int)args_info.molFraction_given){
    sprintf(painCave.errMsg, "Creating a randomized spherical nanoparticle.");
    painCave.isFatal = 0;
    simError();
    /* Random particle is the default case*/

    for (int i = 0; i < sites.size(); i++) 
      if (!isVacancy[i]) ids.push_back(i);
    
    std::random_shuffle(ids.begin(), ids.end());
    
  } else{ 
    sprintf(painCave.errMsg, "Creating a core-shell spherical nanoparticle.");
    painCave.isFatal = 0;
    simError();

    RealType smallestSoFar;
    int myComponent = -1;
    nMol.clear();
    nMol.resize(nComponents);

    for (int i = 0; i < sites.size(); i++) {
      myLoc = sites[i];
      myR = myLoc.length();
      smallestSoFar = particleRadius;      
      if (!isVacancy[i]) {
        for (int j = 0; j < nComponents; j++) {
          if (myR <= shellRadii[j]) {
            if (shellRadii[j] <= smallestSoFar) {
              smallestSoFar = shellRadii[j];
              myComponent = j;
            }
          }
        }
        componentFromSite[i] = myComponent;
        nMol[myComponent]++;
      }
    }       
  }
  
  outputFileName = args_info.output_arg;
   
  //creat new .md file on fly which corrects the number of molecule     
  createMdFile(inputFileName, outputFileName, nMol);
  
  if (oldInfo != NULL)
    delete oldInfo;
  
  SimCreator newCreator;
  SimInfo* NewInfo = newCreator.createSim(outputFileName, false);
    
  // Place molecules
  Molecule* mol;
  SimInfo::MoleculeIterator mi;
  mol = NewInfo->beginMolecule(mi);

  int l = 0;
  int whichSite = 0;

  for (int i = 0; i < nComponents; i++){
    locator = new MoLocator(NewInfo->getMoleculeStamp(i), 
                            NewInfo->getForceField());
    
    if (!args_info.molFraction_given) {
      for (int n = 0; n < sites.size(); n++) {
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
  
  //fill Hmat
  hmat(0, 0)=  10.0*particleRadius;
  hmat(0, 1) = 0.0;
  hmat(0, 2) = 0.0;
  
  hmat(1, 0) = 0.0;
  hmat(1, 1) =  10.0*particleRadius;
  hmat(1, 2) = 0.0;
  
  hmat(2, 0) = 0.0;
  hmat(2, 1) = 0.0;
  hmat(2, 2) =  10.0*particleRadius;
  
  //set Hmat
  NewInfo->getSnapshotManager()->getCurrentSnapshot()->setHmat(hmat);
  
  
  //create dumpwriter and write out the coordinates
  writer = new DumpWriter(NewInfo, outputFileName);
  
  if (writer == NULL) {
    sprintf(painCave.errMsg, "Error in creating dumpwriter object ");
    painCave.isFatal = 1;
    simError();
  }
  
  writer->writeDump();

  // deleting the writer will put the closing at the end of the dump file

  delete writer;

  // cleanup a by calling sim error.....
  sprintf(painCave.errMsg, "A new OpenMD file called \"%s\" has been "
          "generated.\n", outputFileName.c_str());
  painCave.isFatal = 0;
  simError();
  return 0;
}

void createMdFile(const std::string&oldMdFileName, 
                  const std::string&newMdFileName,
                  std::vector<int> nMol) {
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
      if(i<nMol.size()){
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
    sprintf(painCave.errMsg, "Couldn't replace the correct number of nMol\n"
            "\tstatements in component blocks.  Make sure that all\n"
            "\tcomponents in the template file have nMol=1");
    painCave.isFatal = 1;
    simError();
  }
    
}

