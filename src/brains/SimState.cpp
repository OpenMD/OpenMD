#include <iostream>

#include <stdlib.h>
#include <stdio.h>

#include "utils/simError.h"
#include "brains/SimState.hpp"


SimState::SimState(){

  nElements = 0;
  arraysAllocated = false;

  pos  = NULL;
  vel  = NULL;
  frc  = NULL;
  trq  = NULL;
  Amat = NULL;
  mu   = NULL;
  ul   = NULL;
}


SimState::~SimState(){

  if(pos  != NULL) delete[] pos;
  if(vel  != NULL) delete[] vel;
  if(frc  != NULL) delete[] frc;
  if(trq  != NULL) delete[] trq;
  if(Amat != NULL) delete[] Amat;
  if(mu   != NULL) delete[] mu;
  if(ul   != NULL) delete[] ul;  

}


void SimState::createArrays (int the_nElements) {
  int i, index3, index9;

  if( arraysAllocated ){

    sprintf( painCave.errMsg,
	     "Someone has attempted to allocate SimState arrays before "
	     "deallocating the previous arrays.\n" );
    painCave.isFatal = 1;
    simError();
  }

  nElements = the_nElements;

  pos  = new double[nElements*3];
  vel  = new double[nElements*3];
  frc  = new double[nElements*3];
  trq  = new double[nElements*3];
  Amat = new double[nElements*9];
  mu   = new double[nElements];
  ul   = new double[nElements*3];
  
  // init directional values to zero
  
  for( i=0; i<nElements; i++){
    index3 = i*3;
    index9 = i*9;
    
    trq[index3+0] = 0.0;
    trq[index3+1] = 0.0;
    trq[index3+2] = 0.0;
    
    Amat[index9+0] = 1.0;
    Amat[index9+1] = 0.0;
    Amat[index9+2] = 0.0;
    
    Amat[index9+3] = 0.0;
    Amat[index9+4] = 1.0;
    Amat[index9+5] = 0.0;
    
    Amat[index9+6] = 0.0;
    Amat[index9+7] = 0.0;
    Amat[index9+8] = 1.0;
    
    mu[i] = 0.0;    
    
    ul[index3+0] = 1.0;
    ul[index3+1] = 0.0;
    ul[index3+2] = 0.0;
  }

  arraysAllocated = true;
}

void SimState::destroyArrays( void ){

  if(pos  != NULL) delete[] pos;
  if(vel  != NULL) delete[] vel;
  if(frc  != NULL) delete[] frc;
  if(trq  != NULL) delete[] trq;
  if(Amat != NULL) delete[] Amat;
  if(mu   != NULL) delete[] mu;
  if(ul   != NULL) delete[] ul;  

  pos  = NULL;
  vel  = NULL;
  frc  = NULL;
  trq  = NULL;
  Amat = NULL;
  mu   = NULL;
  ul   = NULL;


  arraysAllocated = false;
  nElements = 0;
}

void SimState::getAtomPointers( int index,
				double** the_pos,
				double** the_vel,
				double** the_frc,
				double** the_trq,
				double** the_Amat,
				double** the_mu, 
				double** the_ul){
  int index3, index9;

  if( arraysAllocated ){
    
    index3 = index*3;
    index9 = index*9;
    
    *the_pos  = &(pos[index3]);
    *the_vel  = &(vel[index3]); 
    *the_frc  = &(frc[index3]); 
    *the_trq  = &(trq[index3]); 
    *the_Amat = &(Amat[index9]);
    *the_mu   = &(mu[index]); 
    *the_ul   = &(ul[index3]);
  }
  else{

    sprintf( painCave.errMsg,
	     "Atom %d attempted to access its arrays before they had been"
	     " allocated\n",
	     index );
    painCave.isFatal = 1;
    simError();
  }
}
