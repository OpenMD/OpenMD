#include "Lattice.hpp"
#include "LatticeFactory.hpp"
#include "LatticeCreator.hpp"

static LatticeCreator<FCCLattice> *FCCLatticeCreator = new LatticeCreator<FCCLattice>(FCCLatticeType);
//static LatticeCreator<FCCLattice> *BCCLatticeCreator = new LatticeCreator<FCCLattice>(BCCLatticeType);
//static LatticeCreator<FCCLattice> *HCPLatticeCreator = new LatticeCreator<FCCLattice>(HCPCLatticeType);
//static LatticeCreator<FCCLattice> *OrthorhombicLattice = new LatticeCreator<FCCLattice>(OrthorhombicLatticeType);

CubicLattice::CubicLattice(){
  latticeParam = 1.0;
  
  cellLen.x = latticeParam;
  cellLen.y = latticeParam;
  cellLen.z = latticeParam;
  
}

vector<double> CubicLattice::getLatticeConstant(){
  vector<double> lc;
  
  lc.push_back(cellLen.x);
  return lc;
}

void CubicLattice::setLatticeConstant(const vector<double>& lc){
  
  if(lc.size() < 1){
    cerr << "CubicLattice::setLatticeConstant Error: the size of lattice constant vector  is 0" << endl;
    exit(1);
  }
  else if (lc.size() > 1){
    cerr << "CubicLattice::setLatticeConstant Warning: the size of lattice constant vector  is " << lc.size() << endl;
  }
  
  latticeParam = lc[0];
  
  cellLen.x = latticeParam;
  cellLen.y = latticeParam;
  cellLen.z = latticeParam;
  
  update();
}

FCCLattice::FCCLattice() : CubicLattice(){
  nCellSites = 4;
  cellSitesPos.resize(nCellSites);
  cellSitesOrt.resize(nCellSites);
  update();

}

void FCCLattice::update(){

  double cellLenOver2;
  double oneOverRoot3;

  cellLenOver2  = 0.5 * latticeParam;
  oneOverRoot3 = 1.0 / sqrt(3.0);

  // Molecule 1
  cellSitesPos[0].x = 0.0; 
  cellSitesPos[0].y = 0.0;
  cellSitesPos[0].z = 0.0;
  
   cellSitesOrt[0].x = oneOverRoot3;
   cellSitesOrt[0].y = oneOverRoot3;
   cellSitesOrt[0].z = oneOverRoot3;

  // Molecule 2  
  cellSitesPos[1].x   = 0.0;
  cellSitesPos[1].y   = cellLenOver2;
  cellSitesPos[1].z   = cellLenOver2;

  cellSitesOrt[1].x = -oneOverRoot3;
  cellSitesOrt[1].y = oneOverRoot3;
  cellSitesOrt[1].z = -oneOverRoot3;
   
  // Molecule 3
  cellSitesPos[2].x   = cellLenOver2;
  cellSitesPos[2].y   = cellLenOver2;
  cellSitesPos[2].z   = 0.0;

  cellSitesOrt[2].x = oneOverRoot3;
  cellSitesOrt[2].y = -oneOverRoot3;
  cellSitesOrt[2].z = -oneOverRoot3;

  // Molecule 4

  cellSitesPos[3].x   = cellLenOver2;
  cellSitesPos[3].y   = 0.0;
  cellSitesPos[3].z   = cellLenOver2;

  cellSitesOrt[3].x = -oneOverRoot3;
  cellSitesOrt[3].y = oneOverRoot3;
  cellSitesOrt[3].z = oneOverRoot3;
}

