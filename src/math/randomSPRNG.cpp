#include <iostream>
#include <math.h>

#include "randomSPRNG.hpp"
#include "simError.h"
#include <sprng.h>

#ifdef IS_MPI
#include "mpiSimulation.hpp"
#endif

using namespace std;

/* randomStreamSPRNF creates a new SPRNG stream for random numbers
 */

int randomSPRNG::nStreamsInitialized = 0;

randomSPRNG::randomSPRNG(int iseed){
  int newSeed;
  nStreamsInitialized++;
  newSeed = abs(iseed) + nStreamsInitialized;  
  if( newSeed < 0 ) newSeed = abs( newSeed );

#ifdef IS_MPI

  nSPRNGStreams = mpiSim->getNProcessors();
 
  myStreamNumber = mpiSim->getMyNode();

  

#else

  nSPRNGStreams = 1;
  myStreamNumber = 0;

#endif


  thisStream = init_sprng(GTYPE,myStreamNumber,nSPRNGStreams,
			 newSeed,SPRNG_DEFAULT);
}

randomSPRNG::~randomSPRNG(){
  if ( thisStream != NULL){
    free_sprng(thisStream);
    nStreamsInitialized--;
  }
}


double randomSPRNG::getRandom(){
  return sprng(thisStream);
}


// Gaussian SPRNG class...

double gaussianSPRNG::getGaussian(){
  double ranNum1;
  double ranNum2;
  double gaussianNumber;

  ranNum1 = getRandom();
  ranNum2 = getRandom();

  gaussianNumber = sqrt(-2.0 * log(ranNum1)) * cos(2 * M_PI * ranNum2);

  return gaussianNumber;
}
