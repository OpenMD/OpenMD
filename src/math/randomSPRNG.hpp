#ifndef __randomSPRNG_H
#define __randomSPRNG_H

/* Define the random number generator used by SPRNG
   to be type 3 = Combined Multiple Recursive Generator.
*/
#define GTYPE 3
#ifdef IS_MPI
#define USE_MPI
#endif

class randomSPRNG{
public:
  randomSPRNG(int myseed);
  virtual ~randomSPRNG();

  double getRandom();

protected:
  int *thisStream;
  int myStreamNumber;
  int nSPRNGStreams;
  static int nStreamsInitialized;

};


class gaussianSPRNG : protected randomSPRNG{

public:
  gaussianSPRNG(int iseed):randomSPRNG(iseed){}
  ~gaussianSPRNG(){}

  double getGaussian();

protected:

};



#endif
