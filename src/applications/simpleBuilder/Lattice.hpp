#ifndef _LATTICE_H_
#define _LATTICE_H_
#include "applications/simpleBuilder/BaseLattice.hpp"
#include <string>
#include <vector>
using namespace std;

const string FCCLatticeType = "FCC";
const string BCCLatticeType = "BCC";
const string HCPCLatticeType = "HCP";
const string OrthorhombicLatticeType = "ORTHORHOMBIC";


class CubicLattice : public BaseLattice{
  protected:
    CubicLattice();
  public:
    //get lattice constant of unit cell
    virtual vector<double> getLatticeConstant();

    //set lattice constant of unit cell
    virtual void setLatticeConstant(const vector<double>& lc);
  protected:
    double latticeParam;
};


class FCCLattice : public CubicLattice{
  public:
    FCCLattice();
    virtual const string getLatticeType() {return FCCLatticeType;}    
    virtual void update();
    
};

/*
class BCCLattice : public CubicLattice{
  public:
    BCCLattice();
    virtual const string getLatticeType() {return BCCLatticeType;}
};


class HCPLattice : public BaseLattice{
  public:
    HCPLattice();
    virtual const string getLatticeType() {return HCPCLatticeType;}
};

class OrthorhombicLattice : public BaseLattice{
  public:
    OrthorhombicLattice();
    virtual const string getLatticeType() {return OrthorhombicLatticeType;}
};
*/


#endif
