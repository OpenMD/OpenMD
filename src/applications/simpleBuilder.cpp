#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <fstream>

#include "io/Globals.hpp"
#include "brains/SimInfo.hpp"
#include "brains/SimSetup.hpp"
#include "applications/simpleBuilderCmd.h"
#include "utils/StringUtils.hpp"
#include "applications/LatticeFactory.hpp"
#include "applications/Vector3d.hpp"
#include "applications/MoLocator.hpp"
#include "applications/Lattice.hpp"

using namespace std;

void createMdFile(const string& oldMdFileName, const string& newMdFileName, int numMol);
double getMolMass(MoleculeStamp* molStamp, ForceFields* myFF);

int main( int argc, char* argv[]){

  gengetopt_args_info args_info;
  string latticeType;
  string inputFileName;
  string outPrefix;
  string outMdFileName;
  string outInitFileName;
  SimInfo* oldInfo;
  SimSetup* oldSimSetup;
  BaseLattice* simpleLat;
  int numMol;
  double latticeConstant;
  vector<double> lc;
  double mass;
  const double rhoConvertConst = 1.661;
  double density;
  int nx, ny, nz;
  double Hmat[3][3];
  MoLocator *locator;
  vector<Vector3d> latticePos;
  vector<Vector3d> latticeOrt;
  int numMolPerCell;
  int curMolIndex;
  DumpWriter* writer;
  
  // parse command line arguments
  if (cmdline_parser (argc, argv, &args_info) != 0)
    exit(1) ;
  
  density = args_info.density_arg;

  //get lattice type
  latticeType = UpperCase(args_info.latticetype_arg);
  if(!LatticeFactory::getInstance()->hasLatticeCreator(latticeType)){
    cerr << latticeType << " is an invalid lattice type" << endl;
    cerr << LatticeFactory::getInstance()->toString() << endl;
    exit(1);
  }

  //get the number of unit cell
  nx = args_info.nx_arg;
  if(nx <= 0){
    cerr << "The number of unit cell in h direction must be greater than 0" << endl;
    exit(1);
  }

  ny = args_info.ny_arg;
  if(ny <= 0){
    cerr << "The number of unit cell in l direction must be greater than 0" << endl;
    exit(1);
  }

  nz = args_info.nz_arg;
  if(nz <= 0){
    cerr << "The number of unit cell in k direction must be greater than 0" << endl;
    exit(1);
  }
	
  //get input file name
  if (args_info.inputs_num) 
    inputFileName = args_info.inputs[0];
  else {		 
    cerr <<"You must specify a input file name.\n" << endl;
    cmdline_parser_print_help();
    exit(1);
  }

  
  //parse md file and set up the system
  oldInfo = new SimInfo;
  if(oldInfo == NULL){
     cerr << "error in creating SimInfo" << endl;
     exit(1);
  }

  oldSimSetup = new SimSetup();  
  if(oldSimSetup == NULL){
     cerr << "error in creating SimSetup" << endl;
     exit(1);
  }

  oldSimSetup->suspendInit();
  oldSimSetup->setSimInfo(oldInfo );
  oldSimSetup->parseFile(&inputFileName[0] );
  oldSimSetup->createSim(); 
  
  if(oldInfo->nComponents >=2){
      cerr << "can not build the system with more than two components" << endl;
      exit(1);
  }
  
  //get mass of molecule. 
  //Due to the design of OOPSE, given atom type, we have to query forcefiled to get the mass
  mass = getMolMass(oldInfo->compStamps[0], oldSimSetup->getForceField());
  
  //creat lattice
	simpleLat = LatticeFactory::getInstance()->createLattice(latticeType);
	if(simpleLat == NULL){
		cerr << "Error in creating lattice" << endl;
 		exit(1);
	}

  numMolPerCell = simpleLat->getNumSitesPerCell();
  
  //calculate lattice constant (in Angstrom)
  latticeConstant = pow(rhoConvertConst * numMolPerCell * mass /density, 1.0/3.0);
  
  //set lattice constant
  lc.push_back(latticeConstant);
  simpleLat->setLatticeConstant(lc);
  
  //calculate the total number of molecules
  numMol = nx * ny * nz * numMolPerCell;

  if (oldInfo->n_mol != numMol){

    outPrefix = getPrefix(inputFileName.c_str()) + "_" + latticeType;
    outMdFileName = outPrefix + ".md";

    //creat new .md file on fly which corrects the number of molecule     
    createMdFile(inputFileName, outMdFileName, numMol);
    cerr << "SimpleBuilder Error: the number of molecule and the density are not matched" <<endl;
    cerr << "A new .md file: " << outMdFileName << " is generated, use it to rerun the simpleBuilder" << endl;
    exit(1);
  }
  
  //determine the output file names  
  if (args_info.output_given)
    outInitFileName = args_info.output_arg;
  else
    outInitFileName = getPrefix(inputFileName.c_str())  + ".in";
  
  
  //allocat memory for storing pos, vel and etc
  oldInfo->getConfiguration()->createArrays(oldInfo->n_atoms);
  for (int i = 0; i < oldInfo->n_atoms; i++)
    oldInfo->atoms[i]->setCoords();  

  //creat Molocator
  locator = new MoLocator(oldInfo->compStamps[0], oldSimSetup->getForceField());

  //fill Hmat
  Hmat[0][0] = nx * latticeConstant;
  Hmat[0][1] = 0.0;
  Hmat[0][2] = 0.0;

  Hmat[1][0] = 0.0;
  Hmat[1][1] = ny * latticeConstant;
  Hmat[1][2] = 0.0;

  Hmat[2][0] = 0.0;
  Hmat[2][1] = 0.0;
  Hmat[2][2] = nz * latticeConstant ;

  //set Hmat
  oldInfo->setBoxM(Hmat);
  
  //place the molecules

  curMolIndex = 0;

  //get the orientation of the cell sites
  //for the same type of molecule in same lattice, it will not change
  latticeOrt = simpleLat->getLatticePointsOrt();

  for(int i =0; i < nx; i++){
    for(int j=0; j < ny; j++){
       for(int k = 0; k < nz; k++){

          //get the position of the cell sites
          simpleLat->getLatticePointsPos(latticePos, i, j, k);

          for(int l = 0; l < numMolPerCell; l++)
            locator->placeMol(latticePos[l], latticeOrt[l], &(oldInfo->molecules[curMolIndex++]));
       }
    }
  }

  //create dumpwriter and write out the coordinates
  oldInfo->finalName = outInitFileName;
  writer = new DumpWriter( oldInfo );
  if(writer == NULL){
    cerr << "error in creating DumpWriter" << endl;
    exit(1);    
  }
  writer->writeFinal(0);
  cout << "new initial configuration file: " << outInitFileName <<" is generated." << endl;
  //delete objects

  //delete oldInfo and oldSimSetup
  if(oldInfo != NULL)
     delete oldInfo;
  
  if(oldSimSetup != NULL)
     delete oldSimSetup; 
  
  if (writer != NULL)
    delete writer;
  return 0;
}

void createMdFile(const string& oldMdFileName, const string& newMdFileName, int numMol){
  ifstream oldMdFile;
  ofstream newMdFile;
  const int MAXLEN = 65535;
  char buffer[MAXLEN];

  //create new .md file based on old .md file
  oldMdFile.open(oldMdFileName.c_str());
  newMdFile.open(newMdFileName.c_str());

  oldMdFile.getline(buffer, MAXLEN);
  while(!oldMdFile.eof()){

    //correct molecule number
    if(strstr(buffer, "nMol") !=NULL){      
      sprintf(buffer, "\t\tnMol = %d;", numMol);
      newMdFile << buffer << endl;
    }
    else
      newMdFile << buffer << endl;

    oldMdFile.getline(buffer, MAXLEN);
  }

  oldMdFile.close();
  newMdFile.close();

}

double getMolMass(MoleculeStamp* molStamp, ForceFields* myFF){
  int nAtoms;
  AtomStamp* currAtomStamp;
  double totMass;
  
  totMass = 0;
  nAtoms = molStamp->getNAtoms();

  for(size_t i=0; i<nAtoms; i++){
    currAtomStamp = molStamp->getAtom(i);
    totMass += myFF->getAtomTypeMass(currAtomStamp->getType());
  }

  return totMass;
}
