#include <iostream>

using namespace std;

#include "utils/simError.h"
#include "primitives/Atom.hpp"

Atom::Atom(int theIndex, SimState* theConfig) {

  objType = OT_ATOM;
  myConfig = theConfig;
  hasCoords = false;

  has_dipole = 0;
  has_charge = 0;
  
  index = theIndex;
  offset = 0;
  offsetX = offset;
  offsetY = offset+1;
  offsetZ = offset+2;
  
  Axx = 0;
  Axy = Axx+1;
  Axz = Axx+2;
  
  Ayx = Axx+3;
  Ayy = Axx+4;
  Ayz = Axx+5;
  
  Azx = Axx+6;
  Azy = Axx+7;
  Azz = Axx+8;
}

void Atom::setIndex(int theIndex) {
  index = theIndex; 
}

void Atom::setCoords(void){

  if( myConfig->isAllocated() ){

    myConfig->getAtomPointers( index,
		     &pos, 
		     &vel, 
		     &frc, 
		     &trq, 
		     &Amat,
		     &mu,  
		     &ul);
  }
  else{
    sprintf( painCave.errMsg,
	     "Attempted to set Atom %d  coordinates with an unallocated "
	     "SimState object.\n", index );
    painCave.isFatal = 1;
    simError();
  }
  
  hasCoords = true;
  
}

void Atom::getPos( double theP[3] ){
  
  if( hasCoords ){
    theP[0] = pos[offsetX];
    theP[1] = pos[offsetY];
    theP[2] = pos[offsetZ];
  }
  else{

    sprintf( painCave.errMsg,
	     "Attempt to get Pos for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}

void Atom::setPos( double theP[3] ){

  if( hasCoords ){
    pos[offsetX] = theP[0];
    pos[offsetY] = theP[1];
    pos[offsetZ] = theP[2];
  }
  else{

    sprintf( painCave.errMsg,
	     "Attempt to set Pos for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}

void Atom::getVel( double theV[3] ){
  
  if( hasCoords ){
    theV[0] = vel[offsetX];
    theV[1] = vel[offsetY];
    theV[2] = vel[offsetZ];
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to get vel for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }

}

void Atom::setVel( double theV[3] ){
  
  if( hasCoords ){
    vel[offsetX] = theV[0];
    vel[offsetY] = theV[1];
    vel[offsetZ] = theV[2];
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to set vel for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}

void Atom::getFrc( double theF[3] ){
  
  if( hasCoords ){
    theF[0] = frc[offsetX];
    theF[1] = frc[offsetY];
    theF[2] = frc[offsetZ];
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to get frc for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}

void Atom::addFrc( double theF[3] ){
  
  if( hasCoords ){
    frc[offsetX] += theF[0];
    frc[offsetY] += theF[1];
    frc[offsetZ] += theF[2];
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to add frc for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}


void Atom::zeroForces( void ){
  
  if( hasCoords ){
    frc[offsetX] = 0.0; 
    frc[offsetY] = 0.0; 
    frc[offsetZ] = 0.0;
  }
  else{
    
    sprintf( painCave.errMsg,
	     "Attempt to zero frc for atom %d before coords set.\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}

