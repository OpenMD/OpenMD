#ifndef _ZCONSREADER_H_
#define _ZCONSREADER_H_

#include <fstream>
#include <string>
#include <vector>

#include "integrators/Integrator.hpp"
#include "brains/SimInfo.hpp"

using namespace std;

class ZConsReader{
  public:
    
    ZConsReader(SimInfo* info);
    ZConsReader(const string& filename);
    ~ZConsReader();
    
    void readHeader();
    void readNextFrame();
    bool hasNextFrame();
    int getNumZMol();
    vector<int> getZConsIndex();
    vector<double> getZConsPos();
    //vector<double> getKRatio();
    
    vector<double> getCurZPos(); 
    vector<double> getCurFZ();
    double getCurTime();

  private:
    ifstream* istream;    
    SimInfo* info;
    string zconsFileName;

    int nZMol;
    vector<int> index;
    vector<double> zconsPos;
    //vector<double> kRatio;

    double curTime;    
    vector<double> curFZ;
    vector<double> curZPos;

};

#endif
