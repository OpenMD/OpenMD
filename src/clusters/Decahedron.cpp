/*
 * Copyright (c) 2013 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *                                                                      
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include "clusters/Decahedron.hpp"

using namespace std;

namespace OpenMD {

  Decahedron::Decahedron(int columnAtoms, int shells, int twinAtoms) : 
    N_(columnAtoms), M_(shells), K_(twinAtoms) {
    
    Basis.clear();
    Points.clear();
    
    //
    // Initialize Basis vectors.
    //
    const RealType phi = 2.0 * M_PI / 5.0;  // 72 degrees
    const RealType r3o2 = 0.5 * sqrt(3.0);
    
    Basis.push_back( Vector3d(  r3o2*cos(0.0*phi), r3o2*sin(0.0*phi),  0.0 ));
    Basis.push_back( Vector3d(  r3o2*cos(1.0*phi), r3o2*sin(1.0*phi),  0.0 ));
    Basis.push_back( Vector3d(  r3o2*cos(2.0*phi), r3o2*sin(2.0*phi),  0.0 ));
    Basis.push_back( Vector3d(  r3o2*cos(3.0*phi), r3o2*sin(3.0*phi),  0.0 ));
    Basis.push_back( Vector3d(  r3o2*cos(4.0*phi), r3o2*sin(4.0*phi),  0.0 ));
  }
  
  Decahedron::~Decahedron() {
    Basis.clear();
    Points.clear();
  }
  
  vector<Vector3d> Decahedron::getPoints() {
    // Generate central column of Decahedron

    for (int i = 0; i < N_; i++) {
      Points.push_back( Vector3d( 0.0, 0.0, RealType(i) - 0.5 * (N_ - 1) ) );
    }
    
    for (int i = 1; i < M_ + 1; i++) {
      // generate the shells of the decahedron:
      
      vector<Vector3d> ring;
        
      if (i > K_ - 1) {
        ring = truncatedRing(i, i - K_ + 1);
      } else {
        ring = truncatedRing(i, 0);
      }

      // shift the rings in the z-direction (along the shell)
      
      for (int j = 0; j < N_ - i; j++) {       
        Vector3d shift = Vector3d(0, 0, -0.5 * RealType((N_-i)-1) + RealType(j));
        
        for (vector<Vector3d>::iterator k = ring.begin(); 
             k != ring.end(); ++k) {

          Points.push_back( (*k) + shift);

        }
      }
    }
    return Points;
  }

  vector<Vector3d> Decahedron::truncatedRing( int n, int k ) {
    // This function generates the rings of a Decahedron
    // n: index of shell (order of ring)
    // k: how many atoms are missing from both ends of one side of
    //    pentagon ring

    vector<Vector3d> ring;
    
    // Generate atomic coordinates along each side of pentagonal ring
    for (int i = 0; i < 5; i++) {

      Vector3d b1 = Basis[i];
      Vector3d b2 = Basis[(i + 1) % 5];

      if (k == 0) {
        // without truncation 
        for (int j = 0; j < n; j++) {
          ring.push_back( RealType(n) * b1 + RealType(j) * (b2-b1));
        }
        
      } else {
        for (int j = k; j <= n - k; j++) {
          // with truncation        
          ring.push_back( RealType(n) * b1 + RealType(j) * (b2-b1));
        }
      }
    }
    return ring;
  }

  CurlingStoneDecahedron::CurlingStoneDecahedron(int columnAtoms, int shells,
                                                 int twinAtoms, 
                                                 int truncatedPlanes) :
    Decahedron(columnAtoms, shells, twinAtoms), T_(truncatedPlanes) {}
    
  vector<Vector3d> CurlingStoneDecahedron::getPoints() {

    vector<Vector3d> raw = Decahedron::getPoints();
    vector<Vector3d> snipped;
    RealType maxZ, minZ;   

    maxZ = raw.begin()->z();
    minZ = raw.begin()->z();

    for (vector<Vector3d>::iterator i = raw.begin(); i != raw.end(); ++i) {
      maxZ = max(maxZ, (*i).z());
      minZ = min(minZ, (*i).z());
    }
    
    for (vector<Vector3d>::iterator i = raw.begin(); i != raw.end(); ++i) {
      if ( ((*i).z() < maxZ - 0.995 * (T_ / 2.0) ) && 
           ((*i).z() > minZ + 0.995 * (T_ / 2.0) ) ){
        snipped.push_back( (*i) );
      }
    }
    return snipped;
  }

}
