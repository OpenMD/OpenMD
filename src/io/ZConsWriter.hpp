#ifndef __ZCONSWRITER_H__
#define __ZCONSWRITER_H__

#define _LARGEFILE_SOURCE64
# ifndef _FILE_OFFSET_BITS
#   define _FILE_OFFSET_BITS 64
# endif

#include <iostream>
#include <fstream>

#include "GenericData.hpp"

#ifdef IS_MPI
#include <mpi.h>
#include "mpiSimulation.hpp"
#endif

using namespace std;

class ZConsWriter {

public:
  ZConsWriter(const char* filename, vector<ZConsParaItem>* thePara);
  ~ZConsWriter();  
  
  void writeFZ(double time, int num, int* index, double* fz, double* curZPos, double* zpos);
  
private:
  void writeZPos();
  ofstream output;
  vector<ZConsParaItem>* parameters;
};

#endif
