#ifndef __READWRITE_H__
#define __READWRITE_H__
#define _LARGEFILE_SOURCE64
#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64
#endif
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>


#include "Atom.hpp"
#include "SimInfo.hpp"
#include "Thermo.hpp"
#include "StuntDouble.hpp"

using namespace std;

class DumpWriter{

public:
  DumpWriter( SimInfo* the_entry_plug );
  ~DumpWriter();

  void writeDump( double currentTime );
  void writeFinal( double currentTime);
  void writeFrame( vector<ofstream*>& outFile, double finalTime );

#ifdef IS_MPI
  void update();
#endif

private:

#ifdef IS_MPI
  void sortByGlobalIndex();
#endif

  SimInfo* entry_plug;
  ofstream dumpFile;

#ifdef IS_MPI 
  vector<pair<int, int> > indexArray;
#endif
};

class StatWriter{

public:
  StatWriter( SimInfo* the_entry_plug );
  ~StatWriter();

  void writeStat( double currentTime );

private:

  SimInfo* entry_plug;
  ofstream outFile;
  string outName;
  Thermo* tStats;

};

class InitializeFromFile{

public:
  InitializeFromFile( char *in_name );
  ~InitializeFromFile();

  void readInit( SimInfo* the_entry_plug );

private:
  char* parseDumpLine(char* line, StuntDouble* sd);
  char* parseCommentLine(char* line, SimInfo* entry_plug);
  FILE *c_in_file;
  char c_in_name[500];
  SimInfo* simnfo;
};

class DumpReader{

public:
  DumpReader(const char *in_name );
  ~DumpReader();

  int getNframes();
  void scanFile( void );

  void getNextFrame() {}
  void readFrame(SimInfo* the_simnfo, int whichFrame);

#ifdef IS_MPI
  void anonymousNodeDie( void );
  void nodeZeroError( void );
#endif
private:

  void readSet( int whichFrame );
  char* parseDumpLine(char* readLine, StuntDouble* sd);
  char* parseCommentLine(char* readLine, SimInfo* entry_plug);
  FILE *inFile;
  string inFileName;
  bool isScanned;

  vector<fpos_t*> framePos;
  SimInfo *simnfo;
};



#endif
