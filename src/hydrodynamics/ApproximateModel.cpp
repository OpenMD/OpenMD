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

#include "hydrodynamics/ApproximateModel.hpp"
#include "hydrodynamics/CompositeShape.hpp"
#include "hydrodynamics/Ellipsoid.hpp"
#include "hydrodynamics/Sphere.hpp"
#include "math/DynamicRectMatrix.hpp"
#include "math/LU.hpp"
#include "math/SquareMatrix3.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  /**
   * References:
   *
   * For the General Hydro Framework:  
   * Beatriz Carrasco and Jose Gracia de la Torre; "Hydrodynamic
   * Properties of Rigid Particles: Comparison of Different Modeling
   * and Computational Procedures", Biophysical Journal, 75(6), 3044,
   * 1999
   *
   * Xiuquan Sun, Teng Lin, and J. Daniel Gezelter; "Langevin dynamics
   * for rigid bodies of arbitrary shape", J. Chem. Phys. 128, 234107
   * (2008)
   *
   * For overlapping beads and overlapping volume: 
   *
   * Beatriz Carrasco and Jose Garcia de la Torre and Peter Zipper;
   * "Calculation of hydrodynamic properties of macromolecular bead
   * models with overlapping spheres", Eur Biophys J (1999) 28:
   * 510-515
   *
   * For overlapping volume between two spherical beads:
   * http://mathworld.wolfram.com/Sphere-SphereIntersection.html 
   *
   * For non-overlapping and overlapping translation-translation
   * mobility tensors: 
   *
   * Zuk, P. J., E. Wajnryb, K. A. Mizerski, and P. Szymczak;
   * “Rotne–Prager–Yamakawa Approximation for Different-Sized
   * Particles in Application to Macromolecular Bead Models.”, Journal
   * of Fluid Mechanics, 741 (2014) 
   *
   * For distinctions between centers of resistance and diffusion:
   * Steven Harvey and Jose Garcia de la Torre; "Coordinate Systems
   * for Modeling the Hydrodynamic Resistance and Diffusion
   * Coefficients of Irregularly Shaped Rigid Macromolecules",
   * Macromolecules 1980 13 (4), 960-964
   **/
  ApproximateModel::ApproximateModel() : HydrodynamicsModel() {    
    elements_.clear();
  }

  void ApproximateModel::setShape(Shape* shape) {
    elements_.clear();
    shape_ = shape;
    std::cerr << "Assigning Elements...";
    assignElements();
    std::cerr << " done.\n";
  }

  HydroProp* ApproximateModel::calcHydroProps(RealType viscosity) {
    
    HydroProp* hp = new HydroProp();
    hp->setName( shape_->getName() );
    std::size_t nelements = elements_.size();
    if (nelements == 0)
      nelements = assignElements();   
    
    DynamicRectMatrix<RealType> B(3 * nelements, 3 * nelements);
    DynamicRectMatrix<RealType> C(3 * nelements, 3 * nelements);
    Mat3x3d I = Mat3x3d::identity();

    std::cerr << "Checking " << nelements << " elements...";
    for (std::size_t i = 0; i < nelements; ++i) 
      checkElement(i);
    std::cerr << " done.\n";

    std::cerr << "Calculating B matrix...";
    for (std::size_t i = 0; i < nelements; ++i) {
      for (std::size_t j = 0; j < nelements; ++j) {
        Mat3x3d Tij = interactionTensor(i, j, viscosity);
        B.setSubMatrix(i * 3, j * 3, Tij);
      }
    }
    std::cerr << " done.\n";
        
    // invert B Matrix
    std::cerr << "Inverting B matrix...";
    invertMatrix(B, C);  // B is modified during the inversion
    std::cerr << " done.\n";

    // prepare U Matrix relative to arbitrary origin O(0.0, 0.0, 0.0);
    // relative to the original geometry
    
    std::vector<Mat3x3d> U;
    for (unsigned int i = 0; i < nelements; ++i) {
      Mat3x3d currU;
      currU.setupSkewMat(elements_[i].pos);
      U.push_back(currU);
    }

    // calculate Xi matrix at arbitrary origin O
    Mat3x3d Xiott;
    Mat3x3d Xiorr;
    Mat3x3d Xiotr;

    std::cerr << "Assembling resistance tensor...";
    for (std::size_t i = 0; i < nelements; ++i) {
      for (std::size_t j = 0; j < nelements; ++j) {
        Mat3x3d Cij;
        C.getSubMatrix(i * 3, j * 3, Cij);

        Xiott += Cij;
        Xiotr += U[i] * Cij;
        // Uncorrected here. Volume correction is added after we
        // assemble Xiorr
        Xiorr += -U[i] * Cij * U[j];
      }
    }

    // If the model requires it, calculate the total volume,
    // discounting overlap:
    
    RealType volume = volumeCorrection();

    // Add the volume correction
    Xiorr += (6.0 * viscosity * volume) * I;

    Xiott *= Constants::viscoConvert;
    Xiotr *= Constants::viscoConvert;
    Xiorr *= Constants::viscoConvert;

    // calculate center of resistance
    Mat3x3d tmp;
    Mat3x3d tmpInv;
    Vector3d tmpVec;
    tmp(0, 0) =  Xiott(1, 1) + Xiott(2, 2);
    tmp(0, 1) = -Xiott(0, 1);
    tmp(0, 2) = -Xiott(0, 2);
    tmp(1, 0) = -Xiott(0, 1);
    tmp(1, 1) =  Xiott(0, 0) + Xiott(2, 2);
    tmp(1, 2) = -Xiott(1, 2);
    tmp(2, 0) = -Xiott(0, 2);
    tmp(2, 1) = -Xiott(1, 2);
    tmp(2, 2) =  Xiott(1, 1) + Xiott(0, 0);
    
    tmpVec[0] =  Xiotr(2, 1) - Xiotr(1, 2);
    tmpVec[1] =  Xiotr(0, 2) - Xiotr(2, 0);
    tmpVec[2] =  Xiotr(1, 0) - Xiotr(0, 1);

    // invert tmp Matrix
    invertMatrix(tmp, tmpInv);

    // center of resistance
    Vector3d ror = tmpInv * tmpVec;

    // calculate Resistance Tensor at center of resistance
    Mat3x3d Uor;
    Uor.setupSkewMat(ror);

    Mat3x3d Xirtt;
    Mat3x3d Xirrr;
    Mat3x3d Xirtr;

    // Resistance tensors at the center of resistance
    Xirtt = Xiott;
    Xirtr = (Xiotr - Uor * Xiott);
    Xirrr = Xiorr - Uor * Xiott * Uor + Xiotr * Uor - Uor * Xiotr.transpose();

    // calculate Diffusion tensors at center of resistance
    Mat6x6d Xi;
    Mat6x6d Dr;

    hp->setCenterOfResistance(ror);

    Xi.setSubMatrix(0, 0, Xirtt);
    Xi.setSubMatrix(0, 3, Xirtr.transpose());
    Xi.setSubMatrix(3, 0, Xirtr);
    Xi.setSubMatrix(3, 3, Xirrr);
    std::cerr << " done.\n";

    hp->setResistanceTensor(Xi);
    
    return hp;
  }
}  // namespace OpenMD
