/*
 *  Shapes_FF.cpp
 *  oopse
 *
 *  Created by Chris Fennell on 10/20/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <map>
#include <iostream>

using namespace std;
using namespace oopse;

#ifdef IS_MPI
#include <mpi.h>
#endif //is_mpi

#include "UseTheForce/ForceFields.hpp"
#include "primitives/SRI.hpp"
#include "utils/simError.h"
#include "io/basic_ifstrstream.hpp"
#include "math/RealSphericalHarmonic.hpp"
#include "math/SquareMatrix3.hpp"

#include "UseTheForce/DarkSide/atype_interface.h"
#include "UseTheForce/DarkSide/shapes_interface.h"

#ifdef IS_MPI
#include "UseTheForce/mpiForceField.h"
#endif // is_mpi

Shapes_FF::Shapes_FF() {
  Shapes_FF("");
}

Shapes_FF::Shapes_FF(char* the_variant){
  ffPath_env = "FORCE_PARAM_PATH";

}

Shapes_FF::~Shapes_FF(){

#ifdef IS_MPI
  if( worldRank == 0 ){
#endif // is_mpi
    
    fclose( frcFile );
    
#ifdef IS_MPI
  }
#endif // is_mpi
}


void Shapes_FF::calcRcut( void ){
  
#ifdef IS_MPI
  double tempShapesRcut = shapesRcut;
  MPI_Allreduce( &tempShapesRcut, &shapesRcut, 1, MPI_DOUBLE, MPI_MAX,
		 MPI_COMM_WORLD);
#endif  //is_mpi
  entry_plug->setDefaultRcut(shapesRcut);
}


void Shapes_FF::initForceField(){
  initFortran(0);
}


void Shapes_FF::readForceFile( void ){

  char readLine[1024];
  char fileName[200];
  char temp[200];
  char* ffPath;
  char *shapeFileName;
  char *nameToken;
  char *delim = " ,;\t\n";
  int nTokens, i;
  int nContact = 0;
  int nRange = 0;
  int nStrength = 0;
  int myATID;
  string nameString;
  ShapeType* st;
  map<string, AtomType*> atomTypeMap;
  map<string, AtomType*>::iterator iter;

  // vectors for shape transfer to fortran
  vector<RealSphericalHarmonic> tempSHVector;
  vector<int> contactL;
  vector<int> contactM;
  vector<int> contactFunc;
  vector<double> contactCoeff;
  vector<int> rangeL;
  vector<int> rangeM;
  vector<int> rangeFunc;
  vector<double> rangeCoeff;
  vector<int> strengthL;
  vector<int> strengthM;
  vector<int> strengthFunc;
  vector<double> strengthCoeff;

  // build a basic file reader
  ifstrsteam frcReader;
   
  // generate the force file name   
  strcpy( fileName, "Shapes.frc" );
  
  // attempt to open the file in the current directory first.
  frcReader.open( fileName );
  
  if( frcReader == NULL ){
    
    // next see if the force path enviorment variable is set
    
    ffPath = getenv( ffPath_env );
    if( ffPath == NULL ) {
      STR_DEFINE(ffPath, FRC_PATH );
    }
        
    strcpy( temp, ffPath );
    strcat( temp, "/" );
    strcat( temp, fileName );
    strcpy( fileName, temp );
    
    frcReader.open( fileName );
    
    if( frcFile == NULL ){
      
      sprintf( painCave.errMsg,
	       "Error opening the force field parameter file:\n"
	       "\t%s\n"
	       "\tHave you tried setting the FORCE_PARAM_PATH environment "
	       "variable?\n",
	       fileName );
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  }
  
  // read in the shape types.
  
  findBegin( "ShapeTypes" );
  
  while( !frcReader.eof() ){
    frcReader.getline( readLine, sizeof(readLine) );

    // toss comment lines
    if( readLine[0] != '!' && readLine[0] != '#' ){
      
      if (isEndLine(readLine)) break;
      
      nTokens = count_tokens(readLine, " ,;\t");
      if (nTokens != 0) {
	if (nTokens < 2) {
	  sprintf( painCave.errMsg,
		   "Not enough data on a ShapeTypes line in file: %s\n"
		   fileName );
	  painCave.severity = OOPSE_ERROR;
	  painCave.isFatal = 1;
	  simError();
	}
	
	nameToken = strtok( readLine, delim );
	shapeFileName = strtok( NULL, delim );

	// strings are not char arrays!
	nameString = nameToken;

	// does this AtomType name already exist in the map?
	iter = atomTypeMap.find(nameString);	
	if (iter == atomTypeMap.end()) {
	  // no, it doesn't, so we may proceed:
	  
	  st = new ShapeType();

	  st->setName(nameString);
	  myATID = atomTypeMap.size();
	  st->setIdent(myATID);
	  parseShapeFile(shapeFileName, st);
	  st->complete();
	  atomTypeMap.insert(make_pair(nameString, st));
	  
	} else {
	  // atomType map already contained this string (i.e. it was
	  // declared in a previous block, and we just need to add
	  // the shape-specific information for this AtomType:

	  st = (ShapeType*)(iter->second);
	  parseShapeFile(shapeFileName, st);
	}
      }
    }
  }

#ifdef IS_MPI

  // looks like all the processors have their ShapeType vectors ready...
  sprintf( checkPointMsg,
	   "Shapes_FF shape objects read successfully." );
  MPIcheckPoint();

#endif // is_mpi

  // pack up and send the necessary info to fortran
  for (iter = atomTypeMap.begin(); iter != atomTypeMap.end(); ++iter){

    at = (AtomType*)(iter.second);

    if (at->isDirectional()) {

      dat = (DirectionalAtomType*)at;

      if (dat->isShape()) {

	st = (ShapeAtomType*)at;
	
	contactL.clear();
	contactM.clear();
	contactFunc.clear();
	contactCoeff.clear();
	
	tempSHVector = st->getContactFuncs();
	
	nContact = tempSHVector.size();
	for (i=0; i<nContact; i++){
	  contactL.push_back(tempSHVector[i].getL());
	  contactM.push_back(tempSHVector[i].getM());
	  contactFunc.push_back(tempSHVector[i].getFunctionType());
	  contactCoeff.push_back(tempSHVector[i].getCoefficient());
	}
	
	rangeL.clear();
	rangeM.clear();
	rangeFunc.clear();
	rangeCoeff.clear();
	
	tempSHVector = st->getRangeFuncs();
	
	nRange = tempSHVector.size();
	for (i=0; i<nRange; i++){
	  rangeL.push_back(tempSHVector[i].getL());
	  rangeM.push_back(tempSHVector[i].getM());
	  rangeFunc.push_back(tempSHVector[i].getFunctionType());
	  rangeCoeff.push_back(tempSHVector[i].getCoefficient());
	}
	
	strengthL.clear();
	strengthM.clear();
	strengthFunc.clear();
	strengthCoeff.clear();
	
	tempSHVector = st->getStrengthFuncs();
	
	nStrength = tempSHVector.size();
	for (i=0; i<nStrength; i++){
	  strengthL.push_back(tempSHVector[i].getL());
	  strengthM.push_back(tempSHVector[i].getM());
	  strengthFunc.push_back(tempSHVector[i].getFunctionType());
	  strengthCoeff.push_back(tempSHVector[i].getCoefficient());
	}
	
	isError = 0;
	myATID = at->getIdent();
	makeShape( &nContact, &contactL, &contactM, &contactFunc, 
		   &contactCoeff,
		   &nRange, &rangeL, &rangeM, &rangeFunc, &rangeCoeff,
		   &nStrength, &strengthL, &strengthM,
		   &strengthFunc, &strengthCoeff,
		   &myATID, 
		   &isError);
	if( isError ){
	  sprintf( painCave.errMsg,
		   "Error initializing the \"%s\" shape in fortran\n",
		   marker->first );
	  painCave.isFatal = 1;
	  simError();
	}
      }
    }
  }
  
#ifdef IS_MPI
  sprintf( checkPointMsg,
	   "Shapes_FF atom structures successfully sent to fortran\n" );
  MPIcheckPoint();
#endif // is_mpi

}

void SHAPES_FF::initializeAtoms( int nAtoms, Atom** the_atoms ){
  
  int i,j,k;

  // initialize the atoms
  DirectionalAtom* dAtom;
  AtomType* at;
  DirectionalAtomType* dat;
  double mySigma;
  double ji[3];
  double inertialMat[3][3];
  Mat3x3d momInt;
  string myTypeString;

  for( i=0; i<nAtoms; i++ ){
    
    myTypeString = the_atoms[i]->getType();
    iter = atomTypeMap.find(myTypeString);

    if (iter == atomTypeMap.end()) {
      sprintf( painCave.errMsg, 
	       "AtomType error, %s not found in force file.\n",
	       the_atoms[i]->getType() );
      painCave.isFatal = 1;
      simError();
    } else {

      at = (AtomType*)(iter->second);

      the_atoms[i]->setMass( at->getMass() );
      the_atoms[i]->setIdent( at->getIdent() );

      if( at->isLennardJones() ) {
	mySigma = at->properties->getPropertyByName("sigma");
	if( bigSigma < mySigma ) bigSigma = mySigma;
      }

      the_atoms[i]->setHasCharge(at->isCharge());

      if( at->isDirectional() ){

	dat = (DirectionalAtomType*)at;
	dAtom = (DirectionalAtom *) the_atoms[i];

	momInt = dat->getI();

	// zero out the moments of inertia matrix
	for( j=0; j<3; j++ )
	  for( k=0; k<3; k++ )
	    inertialMat[j][k] = momInt(j,k);
	dAtom->setI( inertialMat );       	

	ji[0] = 0.0;
	ji[1] = 0.0;
	ji[2] = 0.0;
	dAtom->setJ( ji );

	if (dat->isDipole()) {
	  dAtom->setHasDipole( dat->isDipole() );
	  entry_plug->n_dipoles++;
	}	  	
      }
    }
  }
}

void Shapes_FF::initializeBonds( int nBonds, Bond** BondArray,
				 bond_pair* the_bonds ){
  
  if( nBonds ){
    sprintf( painCave.errMsg,
	     "Shapes_FF does not support bonds.\n" );
    painCave.isFatal = 1;
    simError();
  }
}

void Shapes_FF::initializeBends( int nBends, Bend** bendArray,
				 bend_set* the_bends ){
  
  if( nBends ){
    sprintf( painCave.errMsg,
	     "Shapes_FF does not support bends.\n" );
    painCave.isFatal = 1;
    simError();
  }
}

void Shapes_FF::initializeTorsions( int nTorsions, Torsion** torsionArray,
				    torsion_set* the_torsions ){
  
  if( nTorsions ){
    sprintf( painCave.errMsg,
	     "Shapes_FF does not support torsions.\n" );
    painCave.isFatal = 1;
    simError();
  }
}

int Shapes_FF::parseShapeFile(char *shapeFile, ShapeAtomType* st){
  const int MAXLEN = 1024;
  char inLine[MAXLEN];   
  char *token;
  char *delim = " ,;\t\n";
  Mat3x3d momInert;
  RealSphericalHarmonic tempHarmonic;
  vector<RealSphericalHarmonic> functionVector;

  ifstrstream shapeFile;
  
  // first grab the values in the ShapeInfo section
  findBegin(shapeFile, "ShapeInfo");
  
  shapeFile.getline(inLine, MAXLEN);
  while( !shapeFile.eof() ) {
    // toss comment lines
    if( inLine[0] != '!' && inLine[0] != '#' ){
      // end marks section completion
      if (isEndLine(inLine)) break;
      
      nTokens = count_tokens(inLine, delim);
      if (nTokens != 0) {
	if (nTokens < 5) {
	  sprintf( painCave.errMsg,
		   "Not enough data on a ShapeInfo line in file: %s\n"
		   fileName );
	  painCave.severity = OOPSE_ERROR;
	  painCave.isFatal = 1;
	  simError();

	  token = strtok(inLine, delim);
	  token = strtok(NULL, delim);
	  st->setMass(atof(token));
	  token = strtok(NULL, delim);
	  momInert(0,0) = atof(token);
	  token = strtok(NULL, delim);
	  momInert(1,1) = atof(token);
	  token = strtok(NULL, delim);
	  momInert(2,2) = atof(token);
	  st->setI(momInert);
	}
      }
    }
  }

  // now grab the contact functions
  findBegin(shapeFile, "ContactFunctions");
  functionVector.clear();
  
  shapeFile.getline(inLine, MAXLEN);
  while( !shapeFile.eof() ) {
    // toss comment lines
    if( inLine[0] != '!' && inLine[0] != '#' ){
      // end marks section completion
      if (isEndLine(inLine)) break;
      
      nTokens = count_tokens(inLine, delim);
      if (nTokens != 0) {
	if (nTokens < 4) {
	  sprintf( painCave.errMsg,
		   "Not enough data on a ContactFunctions line in file: %s\n"
		   fileName );
	  painCave.severity = OOPSE_ERROR;
	  painCave.isFatal = 1;
	  simError();
	  
	  // read in a spherical harmonic function
	  token = strtok(inLine, delim);
	  tempHarmonic.setL(atoi(token));
	  token = strtok(NULL, delim);
	  tempHarmonic.setM(atoi(token));
	  token = strtok(NULL, delim);
	  if (!strcasecmp("sin",token))
	    tempHarmonic.makeSinFunction();
	  else
	    tempHarmonic.makeCosFunction();		
	  token = strtok(NULL, delim);
	  tempHarmonic.setCoefficient(atof(token));
	  
	  functionVector.push_back(tempHarmonic);
	}
      }
    }
  }

  // pass contact functions to ShapeType
  st->setContactFuncs(functionVector);

  // now grab the range functions
  findBegin(shapeFile, "RangeFunctions");
  functionVector.clear();
  
  shapeFile.getline(inLine, MAXLEN);
  while( !shapeFile.eof() ) {
    // toss comment lines
    if( inLine[0] != '!' && inLine[0] != '#' ){
      // end marks section completion
      if (isEndLine(inLine)) break;
      
      nTokens = count_tokens(inLine, delim);
      if (nTokens != 0) {
	if (nTokens < 4) {
	  sprintf( painCave.errMsg,
		   "Not enough data on a RangeFunctions line in file: %s\n"
		   fileName );
	  painCave.severity = OOPSE_ERROR;
	  painCave.isFatal = 1;
	  simError();
	  
	  // read in a spherical harmonic function
	  token = strtok(inLine, delim);
	  tempHarmonic.setL(atoi(token));
	  token = strtok(NULL, delim);
	  tempHarmonic.setM(atoi(token));
	  token = strtok(NULL, delim);
	  if (!strcasecmp("sin",token))
	    tempHarmonic.makeSinFunction();
	  else
	    tempHarmonic.makeCosFunction();		
	  token = strtok(NULL, delim);
	  tempHarmonic.setCoefficient(atof(token));
	  
	  functionVector.push_back(tempHarmonic);
	}
      }
    }
  }

  // pass range functions to ShapeType
  st->setRangeFuncs(functionVector);

  // finally grab the strength functions
  findBegin(shapeFile, "StrengthFunctions");
  functionVector.clear();
  
  shapeFile.getline(inLine, MAXLEN);
  while( !shapeFile.eof() ) {
    // toss comment lines
    if( inLine[0] != '!' && inLine[0] != '#' ){
      // end marks section completion
      if (isEndLine(inLine)) break;
      
      nTokens = count_tokens(inLine, delim);
      if (nTokens != 0) {
	if (nTokens < 4) {
	  sprintf( painCave.errMsg,
		   "Not enough data on a StrengthFunctions line in file: %s\n"
		   fileName );
	  painCave.severity = OOPSE_ERROR;
	  painCave.isFatal = 1;
	  simError();
	  
	  // read in a spherical harmonic function
	  token = strtok(inLine, delim);
	  tempHarmonic.setL(atoi(token));
	  token = strtok(NULL, delim);
	  tempHarmonic.setM(atoi(token));
	  token = strtok(NULL, delim);
	  if (!strcasecmp("sin",token))
	    tempHarmonic.makeSinFunction();
	  else
	    tempHarmonic.makeCosFunction();		
	  token = strtok(NULL, delim);
	  tempHarmonic.setCoefficient(atof(token));
	  
	  functionVector.push_back(tempHarmonic);
	}
      }
    }
  }

  // pass strength functions to ShapeType
  st->setStrengthFuncs(functionVector);

  // we're done reading from this file
  shapeFile.close();
  return 0;
}

