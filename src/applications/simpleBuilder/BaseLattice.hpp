#ifndef _BASELATTICE_H_
#define _BASELATTICE_H_

#include <vector>
#include "math/Vector3.hpp"

using namespace std;

class BaseLattice{
  protected:
    BaseLattice(){
      Vector3d zeroVector3d(0.0, 0.0, 0.0);
      
      setOrigin(zeroVector3d);
    }
    
  public:

    //virtual destructor of BaseLattice
    virtual ~BaseLattice() {}

    //get lattice type
    virtual const string getLatticeType() = 0;
    
    int getNumSitesPerCell() {return nCellSites;}

    void getLatticePointsPos(vector<Vector3d>& latticePos, int nx, int ny, int nz);

    vector<Vector3d> getLatticePointsOrt() {return cellSitesOrt;}
    
    //get lattice constant of unit cell
    virtual vector<double> getLatticeConstant() =0;

    //set lattice constant of unit cell
    virtual void setLatticeConstant(const vector<double>& lc)=0;

    //get origin of unit cell
    Vector3d getOrigin( ) {return origin;} 

    //set origin of unit cell
    void setOrigin(const Vector3d& newOrigin){
      this->origin = newOrigin;
    }

  protected:
    virtual void update() =0;
    
    int nCellSites;
    Vector3d origin;    
    vector<Vector3d> cellSitesPos;
    vector<Vector3d> cellSitesOrt;
    Vector3d cellLen;
};


#endif
