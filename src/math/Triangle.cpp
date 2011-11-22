/* Copyright (c) 2008, 2010 The University of Notre Dame. All Rights Reserved.
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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [4] , Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). *
 *
 *  Triangle.cpp
 *
 *  Purpose: Provide basic triangle object for OpenMD
 *
 *  Created by Charles F. Vardeman II on 29 July 2008.
 *  @author  Charles F. Vardeman II
 *  @version $Id$
 *
 */

#include "math/Triangle.hpp"

using namespace OpenMD;


Triangle::Triangle() : HaveNormal_(false), HaveUnitNormal_(false), HaveCentroid_(false),
                       HaveArea_(false), area_(0.0), normal_(V3Zero), 
		       centroid_(V3Zero), facetVelocity_(V3Zero), mass_(0.0),
		       a_(V3Zero), b_(V3Zero), c_(V3Zero){
}

void Triangle::addVertices(Vector3d P1, Vector3d P2, Vector3d P3){
  vertices_[0] = P1;
  vertices_[1] = P2;
  vertices_[2] = P3;

  // Compute some quantites like a,b,c
  a_ = P1-P2;
  b_ = P1-P3;
  c_ = P2-P3;  
}


RealType Triangle::computeArea(){
  HaveArea_ = true;
  area_ = getNormal().length() * 0.5;
  return area_;
}
// This should return the normal for our calculations.
Vector3d Triangle::computeNormal(){
  HaveNormal_ = true;
  normal_ = cross(a_,b_);
  return normal_;
}
// This should return the normal for our calculations.
Vector3d Triangle::computeUnitNormal(){
  HaveUnitNormal_ = true;
  unitnormal_ = cross(a_,b_);
  unitnormal_.normalize();
  return unitnormal_;
}

Vector3d Triangle::computeCentroid(){
  HaveCentroid_ = true;
  centroid_ = (vertices_[0] + vertices_[1] + vertices_[2])/3.0;
  return centroid_;
}


Mat3x3d Triangle::computeHydrodynamicTensor(RealType viscosity) {
  
  Vector3d u0 = -a_;
  Vector3d v0 = centroid_ - vertices_[0];
  RealType s0 = 0.5*cross(u0,v0).length();
  
  Vector3d u1 = -c_;
  Vector3d v1 = centroid_ - vertices_[1];
  RealType s1 = 0.5*cross(u1,v1).length();
  
  Vector3d u2 = b_;
  Vector3d v2 = centroid_ - vertices_[2];
  RealType s2 = 0.5*cross(u2,v2).length();
  
  Mat3x3d H;
  H = hydro_tensor(centroid_,centroid_,vertices_[1],vertices_[0],s0,viscosity)+
    hydro_tensor(centroid_,centroid_,vertices_[1],vertices_[2],s1,viscosity)+
    hydro_tensor(centroid_,centroid_,vertices_[2],vertices_[0],s2,viscosity);

  return H.inverse();
}

Mat3x3d Triangle::hydro_tensor(
                               const Vector3d& ri,
                               const Vector3d& rj0,
                               const Vector3d& rj1,
                               const Vector3d& rj2,
                               RealType s, RealType viscosity){
  
  Vector3d v2 = (rj0 + rj1 + rj2)/3.0;  // sub-centroid
  Vector3d dr = ri - v2;                // real centroid to sub-centroid
  RealType l2 = 1.0/dr.lengthSquare();  
 
  Mat3x3d G;
  G = (SquareMatrix3<RealType>::identity() + outProduct(dr,dr)*l2)*sqrt(l2);

  G *= 0.125/3.14159285358979;
  G *= s/viscosity;
  return G;
}
