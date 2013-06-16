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

#include "clusters/Icosahedron.hpp"

using namespace std;

namespace OpenMD {

  Icosahedron::Icosahedron() {
    Basis.clear();
    Edges.clear();
    Facets.clear();
    Points.clear();
    
    //
    // Initialize Basis vectors.
    //
    const RealType HGR = ( sqrt(5.0) + 1.0 ) / 4.0; // half of the golden ratio
    
    Basis.push_back( Vector3d(  HGR,  0.0,  0.5 ));
    Basis.push_back( Vector3d(  HGR,  0.0, -0.5 ));
    Basis.push_back( Vector3d(  0.5,  HGR,  0.0 ));
    Basis.push_back( Vector3d( -0.5,  HGR,  0.0 ));
    Basis.push_back( Vector3d(  0.0,  0.5,  HGR ));
    Basis.push_back( Vector3d(  0.0, -0.5,  HGR ));
    Basis.push_back( Vector3d(  0.5, -HGR,  0.0 ));
    Basis.push_back( Vector3d(  0.0,  0.5, -HGR ));
    Basis.push_back( Vector3d( -HGR,  0.0,  0.5 ));
    Basis.push_back( Vector3d(  0.0, -0.5, -HGR ));
    Basis.push_back( Vector3d( -HGR,  0.0, -0.5 ));
    Basis.push_back( Vector3d( -0.5, -HGR,  0.0 ));
    
    //
    // Initialize 30 edges
    //
    
    Edges.push_back(std::make_pair(0, 1));
    Edges.push_back(std::make_pair(0, 2));
    Edges.push_back(std::make_pair(0, 4));
    Edges.push_back(std::make_pair(0, 5));
    Edges.push_back(std::make_pair(0, 6));
    
    Edges.push_back(std::make_pair(10, 3));
    Edges.push_back(std::make_pair(10, 7));
    Edges.push_back(std::make_pair(10, 8));
    Edges.push_back(std::make_pair(10, 9));
    Edges.push_back(std::make_pair(10, 11));
  
    Edges.push_back(std::make_pair(1, 2));
    Edges.push_back(std::make_pair(1, 6));
    Edges.push_back(std::make_pair(1, 7));
    Edges.push_back(std::make_pair(1, 9));
    
    Edges.push_back(std::make_pair(8, 3));
    Edges.push_back(std::make_pair(8, 4));
    Edges.push_back(std::make_pair(8, 5));
    Edges.push_back(std::make_pair(8, 11));
    
    Edges.push_back(std::make_pair(2, 3));
    Edges.push_back(std::make_pair(2, 4));
    Edges.push_back(std::make_pair(2, 7));
    
    Edges.push_back(std::make_pair(11, 5));
    Edges.push_back(std::make_pair(11, 6));
    Edges.push_back(std::make_pair(11, 9));
    
    Edges.push_back(std::make_pair(6, 5));
    Edges.push_back(std::make_pair(6, 9));
    
    Edges.push_back(std::make_pair(3, 4));
    Edges.push_back(std::make_pair(3, 7));
    
    Edges.push_back(std::make_pair(7, 9));
    
    Edges.push_back(std::make_pair(5, 4));
    
    //
    // Initialize 20 facets
    //
    
    Facets.push_back(make_tuple3(0, 1, 2));
    Facets.push_back(make_tuple3(0, 2, 4));
    Facets.push_back(make_tuple3(0, 4, 5));
    Facets.push_back(make_tuple3(0, 5, 6));
    Facets.push_back(make_tuple3(0, 1, 6));
    
    Facets.push_back(make_tuple3(10, 3, 7));
    Facets.push_back(make_tuple3(10, 3, 8));
    Facets.push_back(make_tuple3(10, 8, 11));
    Facets.push_back(make_tuple3(10, 9, 11));
    Facets.push_back(make_tuple3(10, 7, 9));
    
    Facets.push_back(make_tuple3(1, 2, 7));
    Facets.push_back(make_tuple3(1, 7, 9));
    Facets.push_back(make_tuple3(1, 6, 9));
    
    Facets.push_back(make_tuple3(8, 5, 11));
    Facets.push_back(make_tuple3(8, 4, 5));
    Facets.push_back(make_tuple3(8, 3, 4));
    
    Facets.push_back(make_tuple3(2, 3, 7));
    Facets.push_back(make_tuple3(2, 3, 4));
    
    Facets.push_back(make_tuple3(11, 5, 6));
    Facets.push_back(make_tuple3(11, 6, 9));
  }
  
  Icosahedron::~Icosahedron() {
    Facets.clear();
    Edges.clear();
    Basis.clear();
    Points.clear();
  }
  
  int Icosahedron::getNpoints( int n ) {
    int count=0;
    for( int i = 0; i <= n; i++ ) count += np( i );
    return count;
  }

  int Icosahedron::np( int n ) {
    if( n<0 ) return -1;
    else if( n==0 ) return 1;
    else if( n==1 ) return 12;
    else if( n==2 ) return 42;
    else {
      int count = 0;    
      count += 12;   // edge particles
      count += (n-1)*30;   // side particles
      for( int i = 1; i <= n-2; i++ ) count += i*20;   // body particles
      return count;
    }
  }

  vector<Vector3d> Icosahedron::ih( int n ) {
    
    if( n < 0 ) return Points; 
    
    if( n==0 ) {

      // center particle only

      Points.push_back(Vector3d( 0.0, 0.0, 0.0 ));
      return Points;
    } 
    
    //
    // Generate edge particles
    //

    for( vector<Vector3d>::iterator i = Basis.begin(); i != Basis.end(); ++i ) {
      
      Points.push_back( (*i) * RealType(n) );
    }
    
    //
    // Generate side particles
    //
    
    if( n<2 ) return Points;
    
    for( vector<pair<int,int> >::iterator i=Edges.begin(); 
         i != Edges.end(); ++i ) {
      
      Vector3d e1 = Basis[ (*i).first  ] * RealType(n);
      Vector3d e2 = Basis[ (*i).second ] * RealType(n);
      
      for( int j = 1; j <= n-1; j++ ) {
        Points.push_back( e1 + (e2-e1) * RealType(j) / RealType(n));
      }      
    }
    
    //
    // Generate body particles
    //
    
    if( n<3 ) return Points;
    
    for( vector<tuple3<int,int,int> >::iterator i = Facets.begin();
         i != Facets.end(); ++i) {
      
      Vector3d e1 = Basis[ (*i).first  ] * RealType(n);
      Vector3d e2 = Basis[ (*i).second ] * RealType(n);
      Vector3d e3 = Basis[ (*i).third  ] * RealType(n);
      
      for( int j=1; j<=n-2; j++ ) {
        
        Vector3d v1 = e1 + (e2-e1) * RealType(j+1) / RealType(n);
        Vector3d v2 = e1 + (e3-e1) * RealType(j+1) / RealType(n);
        
        for( int k=1; k<=j; k++ ) {
          Points.push_back(v1 + (v2-v1) * RealType(k) / RealType(j+1));
        }
      }
    }
    return Points;
  }
  
  vector<Vector3d> Icosahedron::getPoints(int nshells) {
    //generate the coordinates
    for( int i = 0; i <= nshells; i++ ) ih( i );
    return Points;
  }
}
