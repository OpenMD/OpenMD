#include "BaseLattice.hpp"
void BaseLattice::getLatticePointsPos(vector<Vector3d>& latticePos, int nx, int ny, int nz){

  latticePos.resize(nCellSites);
                                                           
  for( int i=0;i < nCellSites;i++){

    latticePos[i].x = origin.x + cellSitesPos[i].x + cellLen.x * (double(nx) - 0.5);
    latticePos[i].y = origin.y + cellSitesPos[i].y + cellLen.y * (double(ny) - 0.5);
    latticePos[i].z = origin.z + cellSitesPos[i].z + cellLen.z * (double(nz) - 0.5);    
  }

}

