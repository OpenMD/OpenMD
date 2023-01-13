/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include "hydrodynamics/BoundaryElementModel.hpp"
#include "hydrodynamics/Mesh.hpp"
#include "math/integration/StrangFixCowperTriangleQuadrature.hpp"
#include "math/integration/TriangleQuadrature.hpp"
#include "math/DynamicRectMatrix.hpp"
#include "math/LU.hpp"
#include "math/SquareMatrix3.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  BoundaryElementModel::BoundaryElementModel() : ApproximateModel() { }

  std::size_t BoundaryElementModel::assignElements() {
    if (shape_ != NULL ) {
      if (shape_->isMesh()) {
        createTriangles( dynamic_cast<Mesh*>(shape_) );
      } else {
	snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, 
		 "BoundaryElementModel::assignElements Error: No mesh was given as the shape\n");
	painCave.severity = OPENMD_ERROR;
	painCave.isFatal  = 1;
	simError();
      }
      return elements_.size();
    }
    return 0;
  }

  void BoundaryElementModel::createTriangles(Mesh* m) {
    if (m != NULL ) {
      std::string name = m->getName();
      std::vector<Triangle> facets = m->getFacets();
      for (std::vector<Triangle>::iterator i = facets.begin();
	   i != facets.end(); ++i) {
	HydrodynamicsElement currTri;
	currTri.name = name;
	currTri.pos = (*i).getCentroid();
	currTri.t = (*i);
	elements_.push_back(currTri);	
      }
    }
  }

  
  void BoundaryElementModel::checkElement(std::size_t i) {
  }

  void BoundaryElementModel::writeElements(std::ostream& os) {
    std::vector<HydrodynamicsElement>::iterator iter;
    std::string name = shape_->getName();
    os << "solid" << " " << name << std::endl;
    for (iter = elements_.begin(); iter != elements_.end(); ++iter) {
      Triangle t = iter->t;
      os << "\t" << "facet normal" << " " << t.getUnitNormal() << std::endl;
      os << "\t\t" << "outer loop" << std::endl;
      os << "\t\t\t" << " " << "vertex" << " " << t.vertex1() << std::endl;
      os << "\t\t\t" << " " << "vertex" << " " << t.vertex2() << std::endl;
      os << "\t\t\t" << " " << "vertex" << " " << t.vertex3() << std::endl;
      os << "\t\t" << "endloop" << std::endl;
      os << "\t" << "endfacet" << std::endl;
    }
    os << "endsolid" << " " << name << std::endl;
  }
  
  Mat3x3d BoundaryElementModel::interactionTensor(const std::size_t i,
						  const std::size_t j,
						  const RealType viscosity) {
    
    Mat3x3d B;
    Mat3x3d I = SquareMatrix3<RealType>::identity();

    StrangFixCowperTriangleQuadratureRule rule(6);
        
    Vector3d centroid = elements_[i].pos;
    Triangle t = elements_[j].t;

    auto Tij = [&t, &centroid, &I, &viscosity](const Vector3d& p)  {
      // p are in barycentric coordinates
      Vector3d r = t.barycentricToCartesian(p);
      Vector3d ab = centroid - r;
      RealType abl = ab.length();
      Mat3x3d T;
      T = (I + outProduct(ab, ab) / (abl*abl) );
      T /= (8.0 * Constants::PI * viscosity * abl);
      return T;
    };

    B = TriangleQuadrature<RectMatrix<RealType,3,3>, RealType>::Integrate(Tij,
									  rule,
									  1.0);
    
    centroid = elements_[j].pos;
    t = elements_[i].t;

    auto Tji = [&t, &centroid, &I, &viscosity](const Vector3d& p) {
      // p are in barycentric coordinates
      Vector3d r = t.barycentricToCartesian(p);
      Vector3d ab = centroid - r;
      RealType abl = ab.length();
      Mat3x3d T;
      T = (I + outProduct(ab, ab) / (abl*abl) );
      T /= (8.0 * Constants::PI * viscosity * abl);
      return T;
    };

    B += TriangleQuadrature<RectMatrix<RealType,3,3>, RealType>::Integrate(Tji,
									   rule,
									   1.0);
    B *= 0.5;
    return B;
  }
}  // namespace OpenMD
