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
#include <cmath>
#include <iostream>

#ifdef IS_MPI
#include <mpi.h>
#endif //is_mpi

#include "UseTheForce/ForceFields.hpp"
#include "primitives/SRI.hpp"
#include "utils/simError.h"
#include "utils/StringUtils.hpp"
#include "io/basic_ifstrstream.hpp"
#include "math/RealSphericalHarmonic.hpp"
#include "math/SquareMatrix3.hpp"
#include "types/ShapeAtomType.hpp"
#include "UseTheForce/DarkSide/atype_interface.h"
#include "UseTheForce/DarkSide/shapes_interface.h"

#ifdef IS_MPI
#include "UseTheForce/mpiForceField.h"
#endif // is_mpi

using namespace std;
using namespace oopse;

Shapes_FF::~Shapes_FF(){

#ifdef IS_MPI
  if( worldRank == 0 ){
#endif // is_mpi
    
    forceFile.close();
    
#ifdef IS_MPI
  }
#endif // is_mpi
}


void Shapes_FF::calcRcut( void ){
  
#ifdef IS_MPI
  double tempShapesRcut = bigContact;
  MPI_Allreduce( &tempShapesRcut, &shapesRcut, 1, MPI_DOUBLE, MPI_MAX,
		 MPI_COMM_WORLD);
#else
  shapesRcut = bigContact;
#endif  //is_mpi
  entry_plug->setDefaultRcut(shapesRcut);
}


void Shapes_FF::initForceField(){
  initFortran(0);
}


void Shapes_FF::readParams( void ){

  char readLine[1024];

  string fileName;
  string shapeFileName;
  string tempString;

  char *nameToken;
  char *delim = " ,;\t\n";
  int nTokens, i;
  int nContact = 0;
  int nRange = 0;
  int nStrength = 0;
  int myATID;
  int isError;
  string nameString;
  AtomType* at;
  DirectionalAtomType* dat;
  ShapeAtomType* st;

  map<string, AtomType*>::iterator iter;

  // vectors for shape transfer to fortran
  vector<RealSphericalHarmonic*> tempSHVector;
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
   
  // generate the force file name   
  fileName = "Shapes.frc";
  
  // attempt to open the file in the current directory first.
  forceFile.open( fileName.c_str() );
  
  if( forceFile == NULL ){
    
    tempString = ffPath;
    tempString += "/";
    tempString += fileName;
    fileName = tempString;
    
    forceFile.open( fileName.c_str() );
    
    if( forceFile == NULL ){
      
      sprintf( painCave.errMsg,
	       "Error opening the force field parameter file:\n"
	       "\t%s\n"
	       "\tHave you tried setting the FORCE_PARAM_PATH environment "
	       "variable?\n",
	       fileName.c_str() );
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  }
  
  // read in the shape types.
  
  findBegin( forceFile, "ShapeTypes" );
  
  while( !forceFile.eof() ){
    forceFile.getline( readLine, sizeof(readLine) );

    // toss comment lines
    if( readLine[0] != '!' && readLine[0] != '#' ){
      
      if (isEndLine(readLine)) break;
      
      nTokens = countTokens(readLine, " ,;\t");
      if (nTokens != 0) {
	if (nTokens < 2) {
	  sprintf( painCave.errMsg,
		   "Not enough data on a ShapeTypes line in file: %s\n",
		   fileName.c_str() );
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
	  
	  st = new ShapeAtomType();

	  st->setName(nameString);
	  myATID = atomTypeMap.size() + 1;
	  st->setIdent(myATID);
	  parseShapeFile(shapeFileName, st);
	  st->complete();
	  atomTypeMap.insert(make_pair(nameString, st));
	  
	} else {
	  // atomType map already contained this string (i.e. it was
	  // declared in a previous block, and we just need to add
	  // the shape-specific information for this AtomType:

	  st = (ShapeAtomType*)(iter->second);
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

    at = (AtomType*)(iter->second);

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
	  contactL.push_back(tempSHVector[i]->getL());
	  contactM.push_back(tempSHVector[i]->getM());
	  contactFunc.push_back(tempSHVector[i]->getFunctionType());
	  contactCoeff.push_back(tempSHVector[i]->getCoefficient());
	}
	
	rangeL.clear();
	rangeM.clear();
	rangeFunc.clear();
	rangeCoeff.clear();
	
	tempSHVector = st->getRangeFuncs();
	
	nRange = tempSHVector.size();
	for (i=0; i<nRange; i++){
	  rangeL.push_back(tempSHVector[i]->getL());
	  rangeM.push_back(tempSHVector[i]->getM());
	  rangeFunc.push_back(tempSHVector[i]->getFunctionType());
	  rangeCoeff.push_back(tempSHVector[i]->getCoefficient());
	}
	
	strengthL.clear();
	strengthM.clear();
	strengthFunc.clear();
	strengthCoeff.clear();
	
	tempSHVector = st->getStrengthFuncs();
	
	nStrength = tempSHVector.size();
	for (i=0; i<nStrength; i++){
	  strengthL.push_back(tempSHVector[i]->getL());
	  strengthM.push_back(tempSHVector[i]->getM());
	  strengthFunc.push_back(tempSHVector[i]->getFunctionType());
	  strengthCoeff.push_back(tempSHVector[i]->getCoefficient());
	}
	
	isError = 0;
	myATID = at->getIdent();

	makeShape( &nContact, &contactL[0], &contactM[0], &contactFunc[0], 
		   &contactCoeff[0],
		   &nRange, &rangeL[0], &rangeM[0], &rangeFunc[0], 
                   &rangeCoeff[0],
		   &nStrength, &strengthL[0], &strengthM[0],
		   &strengthFunc[0], &strengthCoeff[0],
		   &myATID, 
		   &isError);

	if( isError ){
	  sprintf( painCave.errMsg,
		   "Error initializing the \"%s\" shape in fortran\n",
		   (iter->first).c_str() );
	  painCave.isFatal = 1;
	  simError();
	}
      }
    }
  }
  
  isError = 0;
  completeShapeFF(&isError);
  if( isError ){
    sprintf( painCave.errMsg,
	     "Error completing Shape FF in fortran\n");
    painCave.isFatal = 1;
    simError();
  }
  
#ifdef IS_MPI
  sprintf( checkPointMsg,
	   "Shapes_FF atom structures successfully sent to fortran\n" );
  MPIcheckPoint();
#endif // is_mpi

}

void Shapes_FF::cleanMe( void ){

}

void Shapes_FF::initializeAtoms( int nAtoms, Atom** the_atoms ){
  
  int i,j,k;
  map<string, AtomType*>::iterator iter;

  // initialize the atoms
  DirectionalAtom* dAtom;
  AtomType* at;
  DirectionalAtomType* dat;
  ShapeAtomType* sat;
  double longCutoff;
  double ji[3];
  double inertialMat[3][3];
  Mat3x3d momInt;
  string myTypeString;

  bigContact = 0.0;

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
      
      if ( at->isShape() ) {
        
        sat = (ShapeAtomType*)at;
        longCutoff = findCutoffDistance(sat);
        if (longCutoff > bigContact) bigContact = longCutoff;  
	cout << bigContact << " is the cutoff value\n";
	
	entry_plug->useShapes = 1;
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

void Shapes_FF::parseShapeFile(string shapeFileName, ShapeAtomType* st){
  const int MAXLEN = 1024;
  char inLine[MAXLEN];   
  char *token;
  char *delim = " ,;\t\n";
  int nTokens;
  Mat3x3d momInert;
  RealSphericalHarmonic* rsh;
  vector<RealSphericalHarmonic*> functionVector;
  ifstrstream shapeFile;
  string tempString;

  shapeFile.open( shapeFileName.c_str() );
  
  if( shapeFile == NULL ){
    
    tempString = ffPath;
    tempString += "/";
    tempString += shapeFileName;
    shapeFileName = tempString;
        
    shapeFile.open( shapeFileName.c_str() );
    
    if( shapeFile == NULL ){
      
      sprintf( painCave.errMsg,
	       "Error opening the shape file:\n"
	       "\t%s\n"
	       "\tHave you tried setting the FORCE_PARAM_PATH environment "
	       "variable?\n",
	       shapeFileName.c_str() );
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  }
  
  // read in the shape types. 
  
  // first grab the values in the ShapeInfo section
  findBegin( shapeFile, "ShapeInfo");
  
  shapeFile.getline(inLine, MAXLEN);
  while( !shapeFile.eof() ) {
    // toss comment lines
    if( inLine[0] != '!' && inLine[0] != '#' ){
      // end marks section completion
      if (isEndLine(inLine)) break;
      
      nTokens = countTokens(inLine, delim);
      if (nTokens != 0) {
	if (nTokens < 5) {
	  sprintf( painCave.errMsg,
		   "Not enough data on a ShapeInfo line in file: %s\n",
		   shapeFileName.c_str() );
	  painCave.severity = OOPSE_ERROR;
	  painCave.isFatal = 1;
	  simError();
	} else {
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
    shapeFile.getline(inLine, MAXLEN); 
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
      nTokens = countTokens(inLine, delim);
      if (nTokens != 0) {
	if (nTokens < 4) {
	  sprintf( painCave.errMsg,
		   "Not enough data on a ContactFunctions line in file: %s\n",
		   shapeFileName.c_str() );
	  painCave.severity = OOPSE_ERROR;
	  painCave.isFatal = 1;
	  simError();
	} else {
	  // read in a spherical harmonic function
	  token = strtok(inLine, delim);
          rsh = new RealSphericalHarmonic();
	  rsh->setL(atoi(token));
	  token = strtok(NULL, delim);
	  rsh->setM(atoi(token));
	  token = strtok(NULL, delim);
	  if (!strcasecmp("sin",token))
	    rsh->makeSinFunction();
	  else
	    rsh->makeCosFunction();		
	  token = strtok(NULL, delim);
	  rsh->setCoefficient(atof(token));
	  
	  functionVector.push_back(rsh);
	}
      }
    }
    shapeFile.getline(inLine, MAXLEN);
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
      
      nTokens = countTokens(inLine, delim);
      if (nTokens != 0) {
	if (nTokens < 4) {
	  sprintf( painCave.errMsg,
		   "Not enough data on a RangeFunctions line in file: %s\n",
		   shapeFileName.c_str() );
	  painCave.severity = OOPSE_ERROR;
	  painCave.isFatal = 1;
	  simError();
	} else {
	  
	  // read in a spherical harmonic function
	  token = strtok(inLine, delim);

          rsh = new RealSphericalHarmonic();
	  rsh->setL(atoi(token));
	  token = strtok(NULL, delim);
	  rsh->setM(atoi(token));
	  token = strtok(NULL, delim);
	  if (!strcasecmp("sin",token))
	    rsh->makeSinFunction();
	  else
	    rsh->makeCosFunction();		
	  token = strtok(NULL, delim);
	  rsh->setCoefficient(atof(token));
	  
	  functionVector.push_back(rsh);
	}
      }
    }
    shapeFile.getline(inLine, MAXLEN);
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
      
      nTokens = countTokens(inLine, delim);
      if (nTokens != 0) {
	if (nTokens < 4) {
	  sprintf( painCave.errMsg,
		   "Not enough data on a StrengthFunctions line in file: %s\n",
		   shapeFileName.c_str() );
	  painCave.severity = OOPSE_ERROR;
	  painCave.isFatal = 1;
	  simError();
	} else {
	  
	  // read in a spherical harmonic function
	  token = strtok(inLine, delim);
          rsh = new RealSphericalHarmonic();
	  rsh->setL(atoi(token));
	  token = strtok(NULL, delim);
	  rsh->setM(atoi(token));
	  token = strtok(NULL, delim);
	  if (!strcasecmp("sin",token))
	    rsh->makeSinFunction();
	  else
	    rsh->makeCosFunction();		
	  token = strtok(NULL, delim);
	  rsh->setCoefficient(atof(token));
	  
	  functionVector.push_back(rsh);
	}
      }
    }
    shapeFile.getline(inLine, MAXLEN);
  }

  // pass strength functions to ShapeType
  st->setStrengthFuncs(functionVector);

  // we're done reading from this file
  shapeFile.close();
}

double Shapes_FF::findLargestContactDistance(ShapeAtomType* st) {
  int i, j,  nSteps;
  double theta, thetaStep, thetaMin, costheta;
  double phi, phiStep;
  double sigma, bs;
  
  nSteps = 16;

  thetaStep = M_PI / nSteps;
  thetaMin = thetaStep / 2.0;
  phiStep = thetaStep * 2.0;
  bs = 0.0;
  
  for (i = 0; i < nSteps; i++) {
    
    theta = thetaMin + i * thetaStep;
    costheta = cos(theta);

    for (j = 0; j < nSteps; j++) {

      phi = j*phiStep;

      sigma = st->getContactValueAt(costheta, phi);
      
      if (sigma > bs) bs = sigma;
    }
  }

  return bs;  
}


double Shapes_FF::findCutoffDistance(ShapeAtomType* st) {
  int i, j,  nSteps;
  double theta, thetaStep, thetaMin, costheta;
  double phi, phiStep;
  double sigma, range;
  double bigCut, tempCut;

  nSteps = 16;

  thetaStep = M_PI / nSteps;
  thetaMin = thetaStep / 2.0;
  phiStep = thetaStep * 2.0;
  bigCut = 0.0;
  
  for (i = 0; i < nSteps; i++) {
    
    theta = thetaMin + i * thetaStep;
    costheta = cos(theta);

    for (j = 0; j < nSteps; j++) {

      phi = j*phiStep;

      sigma = st->getContactValueAt(costheta, phi);
      range = st->getRangeValueAt(costheta, phi);

       // cutoff for a shape is taken to be (2.5*rangeVal + contactVal)
      tempCut = 1.5*range + sigma;

      if (tempCut > bigCut) bigCut = tempCut; 
    }
  }
 
  return bigCut;  
}
