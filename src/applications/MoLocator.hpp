#ifndef __MOLOCATOR_H__
#define __MOLOCATOR_H__

#include <vector>
#include "primitives/Atom.hpp"
#include "primitives/DirectionalAtom.hpp"
#include "types/MoleculeStamp.hpp"
#include "primitives/Molecule.hpp"
#include "Vector3d.hpp"
#include "UseTheForce/ForceFields.hpp"
using namespace std;

//convert lattice vector to rotation matrix
void latVec2RotMat(const Vector3d& lv, double rotMat[3][3]);

class MoLocator{
  
public:
  
  MoLocator( MoleculeStamp* theStamp, ForceFields* theFF);

  void placeMol( const Vector3d& offset, const Vector3d& ort, Molecule* mol);

private:
  
  void calcRefCoords( void );
  
  MoleculeStamp* myStamp;
  ForceFields* myFF;

  vector<Vector3d> refCoords;
  int nIntegrableObjects;
};

#endif
