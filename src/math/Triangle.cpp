/* Copyright (c) 2008 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 *
 *  Triangle.cpp
 *
 *  Purpose: Provide basic triangle object for OOPSE
 *
 *  Created by Charles F. Vardeman II on 29 July 2008.
 *  @author  Charles F. Vardeman II
 *  @version $Id: Triangle.cpp,v 1.2 2008-12-05 16:20:38 chuckv Exp $
 *
 */

#include "math/Triangle.hpp"

using namespace oopse;


Triangle::Triangle() : HaveNormal_(false), HaveCentroid_(false),HaveArea_(false), area_(0.0), normal_(V3Zero), 
		       centroid_(V3Zero),facetVelocity_(V3Zero), mass_(0.0),
		       a_(V3Zero),b_(V3Zero),c_(V3Zero){
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

Vector3d Triangle::computeNormal(){
  HaveNormal_ = true;
  normal_ = cross(a_,b_);
  return normal_;
}

Vector3d Triangle::computeCentroid(){
  HaveCentroid_ = true;
  centroid_ = (vertices_[0] + vertices_[1] + vertices_[2])/3.0;
  return centroid_;
}
