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
 *  Triangle.hpp
 *
 *  Purpose: Provide basic triangle class for OpenMD. Hates Particle class.
 *
 *  Created by Charles F. Vardeman II on 29 July 2008.
 *  @author  Charles F. Vardeman II
 *  @version $Id$
 *
 */


#ifndef MATH_FACET_HPP
#define MATH_FACET_HPP

#include "math/Vector3.hpp"
#include "math/SquareMatrix3.hpp"
#include "config.h"
#include "primitives/StuntDouble.hpp"

#include <vector>


namespace OpenMD {

/**
   * @class Triangle
   *
   * Triangle provides geometric data to OpenMD. Triangle includes
   * information about the normal, centroid and the atoms
   * that belong to this triangle.
   */
  class Triangle {

  public:
    Triangle();
    virtual ~Triangle() { };

    void setNormal(Vector3d normal) {
      normal_ = normal;
      HaveNormal_ = true;
    }
    void setUnitNormal(Vector3d normal) {
      unitnormal_ = normal;
      HaveUnitNormal_ = true;
    }

    void addVertices(Vector3d P1, Vector3d P2, Vector3d P3);

    void addVertexSD(StuntDouble* thisSD){
      vertexSD_.push_back(thisSD);
    }

    std::vector<StuntDouble*> getVertices(){return vertexSD_;}

    void setArea(RealType area) {
      area_ = area;
      HaveArea_ = true;
    }
    
    Vector3d getNormal() {
      if (HaveNormal_) {
	return normal_;
      } else {
	return computeNormal();
      }
    }
   Vector3d getUnitNormal() {
      if (HaveUnitNormal_) {
	return unitnormal_;
      } else {
	return computeUnitNormal();
      }
    }
    
    RealType getArea() { 
      if(HaveArea_){
	return area_;
      }else{
	return computeArea();
      }
    }

    RealType computeArea();
    Vector3d computeNormal();
    Vector3d computeCentroid();
    Vector3d computeUnitNormal();

    void setCentroid(Vector3d centroid) { 
      centroid_ = centroid;
      HaveCentroid_ = true;
    }

    Vector3d getCentroid() {
      if (HaveCentroid_) {
	return centroid_;
      } else {
	return computeCentroid();
      }
    }
    
    Vector3d getFacetVelocity(){
      return facetVelocity_;
    }
    
    void setFacetVelocity(Vector3d facetVelocity){
      facetVelocity_ = facetVelocity;
    }

    void setFacetMass(RealType mass){
      mass_ = mass;
    }

    RealType getFacetMass(){
      return mass_;
    }

    RealType a(){
      return a_.length();
    }

    RealType b(){
      return b_.length();
    }

    RealType c(){
      return c_.length();
    }

    RealType getHydroLength() {
      RealType a1 = a();
      RealType b1 = b();
      RealType c1 = c();
      RealType t1 =  a1 + b1 + c1;
      RealType t4 =  a1 + b1 - c1;

      return 32.0 * c1 / log(t1*t1/t4/t4);
    }


    RealType getIncircleRadius() {
      return 2.0 * getArea() / (a() + b() + c());
    }

    RealType getCircumcircleRadius() {
      RealType a1 = a();
      RealType b1 = b();
      RealType c1 = c();
      RealType t1 =  a1 + b1 + c1;
      RealType t2 = -a1 + b1 + c1;
      RealType t3 =  a1 - b1 + c1;
      RealType t4 =  a1 + b1 - c1;
      RealType junk = t1*t2*t3*t4;
      return a1 * b1 * c1 / sqrt(t1 * t2 * t3 * t4);
    }
      
    Mat3x3d computeHydrodynamicTensor(RealType viscosity);


  private:
    Mat3x3d hydro_tensor(const Vector3d& ri, const Vector3d& rj0, const Vector3d& rj1, const Vector3d& rj2,RealType s, RealType viscosity);
    
    /* Local Indentity of vertex atoms in pos array*/
    std::vector <StuntDouble*> vertexSD_;
    Vector3d normal_;
    Vector3d unitnormal_;
    Vector3d centroid_;
    Vector3d vertices_[3];
    RealType area_;    
    RealType mass_;
    Vector3d facetVelocity_;
    //Length of triangle sides
    Vector3d a_,b_,c_;
    RealType alpha_,beta_,gamma_;
    bool HaveArea_;
    bool HaveNormal_;
    bool HaveUnitNormal_;
    bool HaveCentroid_;
    
  }; // End class Triangle
    
    

} //end namespace OpenMD



#endif // MATH_FACET_HPP





