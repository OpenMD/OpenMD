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

#include "applications/simpleBuilder/simpleBuilderCmd.h"
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
void createMdFile(const std::string&oldMdFileName, const std::string&newMdFileName,
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
    int nx,
    ny,
    nz;
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

    //get the number of unit cell
    nx = args_info.nx_arg;

    if (nx <= 0) {
        std::cerr << "The number of unit cell in h direction must be greater than 0" << std::endl;
        exit(1);
    }

    ny = args_info.ny_arg;

    if (ny <= 0) {
        std::cerr << "The number of unit cell in l direction must be greater than 0" << std::endl;
        exit(1);
    }

    nz = args_info.nz_arg;

    if (nz <= 0) {
        std::cerr << "The number of unit cell in k direction must be greater than 0" << std::endl;
        exit(1);
    }

    //get input file name
    if (args_info.inputs_num)
        inputFileName = args_info.inputs[0];
    else {
        std::cerr << "You must specify a input file name.\n" << std::endl;
        cmdline_parser_print_help();
        exit(1);
    }

    //parse md file and set up the system
    SimCreator oldCreator;
    SimInfo* oldInfo = oldCreator.createSim(inputFileName, false);

    if (oldInfo->getNMoleculeStamp()>= 2) {
        std::cerr << "can not build the system with more than two components"
            << std::endl;
        exit(1);
    }

    //get mass of molecule. 

    mass = getMolMass(oldInfo->getMoleculeStamp(0), oldInfo->getForceField());

    //creat lattice
    simpleLat = LatticeFactory::getInstance()->createLattice(latticeType);

    if (simpleLat == NULL) {
        std::cerr << "Error in creating lattice" << std::endl;
        exit(1);
    }

    numMolPerCell = simpleLat->getNumSitesPerCell();

    //calculate lattice constant (in Angstrom)
    latticeConstant = pow(rhoConvertConst * numMolPerCell * mass / density,
                          1.0 / 3.0);

    //set lattice constant
    lc.push_back(latticeConstant);
    simpleLat->setLatticeConstant(lc);

    //calculate the total number of molecules
    numMol = nx * ny * nz * numMolPerCell;

    if (oldInfo->getNGlobalMolecules() != numMol) {
        outPrefix = getPrefix(inputFileName.c_str()) + "_" + latticeType;
        outMdFileName = outPrefix + ".md";

        //creat new .md file on fly which corrects the number of molecule     
        createMdFile(inputFileName, outMdFileName, numMol);
        std::cerr
            << "SimpleBuilder Error: the number of molecule and the density are not matched"
            << std::endl;
        std::cerr << "A new .md file: " << outMdFileName
            << " is generated, use it to rerun the simpleBuilder" << std::endl;
	exit(1);
    }

    //determine the output file names  
    if (args_info.output_given)
        outInitFileName = args_info.output_arg;
    else
        outInitFileName = getPrefix(inputFileName.c_str()) + ".in";
    
    //creat Molocator
    locator = new MoLocator(oldInfo->getMoleculeStamp(0), oldInfo->getForceField());

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
    oldInfo->getSnapshotManager()->getCurrentSnapshot()->setHmat(hmat);

    //place the molecules

    curMolIndex = 0;

    //get the orientation of the cell sites
    //for the same type of molecule in same lattice, it will not change
    latticeOrt = simpleLat->getLatticePointsOrt();

    Molecule* mol;
    SimInfo::MoleculeIterator mi;
    mol = oldInfo->beginMolecule(mi);
    for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < nz; k++) {

                //get the position of the cell sites
                simpleLat->getLatticePointsPos(latticePos, i, j, k);

                for(int l = 0; l < numMolPerCell; l++) {
                    if (mol != NULL) {
                        locator->placeMol(latticePos[l], latticeOrt[l], mol);
                    } else {
                        std::cerr << std::endl;                    
                    }
                    mol = oldInfo->nextMolecule(mi);
                }
            }
        }
    }

    //create dumpwriter and write out the coordinates
    oldInfo->setFinalConfigFileName(outInitFileName);
    writer = new DumpWriter(oldInfo);

    if (writer == NULL) {
        std::cerr << "error in creating DumpWriter" << std::endl;
        exit(1);
    }

    writer->writeDump();
    std::cout << "new initial configuration file: " << outInitFileName
        << " is generated." << std::endl;

    //delete objects

    //delete oldInfo and oldSimSetup
    if (oldInfo != NULL)
        delete oldInfo;

    if (writer != NULL)
        delete writer;
    
    delete simpleLat;

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
            sprintf(buffer, "\t\tnMol = %d;", numMol);
            newMdFile << buffer << std::endl;
        } else
            newMdFile << buffer << std::endl;

        oldMdFile.getline(buffer, MAXLEN);
    }

    oldMdFile.close();
    newMdFile.close();
}

