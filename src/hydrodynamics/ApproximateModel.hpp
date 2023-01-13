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

#ifndef HYDRODYNAMICS_APPROXIMATEMODEL_HPP
#define HYDRODYNAMICS_APPROXIMATEMODEL_HPP

#include <vector>

#include "hydrodynamics/HydrodynamicsModel.hpp"
#include "math/DynamicRectMatrix.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "math/Triangle.hpp"

namespace OpenMD {

  struct HydrodynamicsElement {
    Vector3d pos;        // minimally, this is all we need
    std::string name;
    RealType radius;     // for Bead models
    Triangle t;          // for BoundaryElement model
    RealType mass;       // for center of mass calculation only
  };

  class Shape;
  class ApproximateModel : public HydrodynamicsModel {
  public:
    ApproximateModel();
    
    virtual void setShape(Shape* shape);
    virtual HydroProp* calcHydroProps(RealType viscosity);

    virtual std::size_t assignElements() = 0;
    virtual void checkElement(std::size_t i) = 0;
    virtual void writeElements(std::ostream& os) = 0;
    virtual Mat3x3d interactionTensor(std::size_t i, std::size_t j,
				      RealType viscosity) = 0;
    
    virtual RealType volumeCorrection() { return 0.0; }

  protected:
    std::vector<HydrodynamicsElement> elements_;

  };
}  // namespace OpenMD

#endif
