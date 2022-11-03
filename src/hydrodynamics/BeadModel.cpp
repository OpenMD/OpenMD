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

#include "hydrodynamics/BeadModel.hpp"
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
  BeadModel::BeadModel() : ApproximateModel() { volumeOverlap_ = 0.0; }

  void BeadModel::checkElement(std::size_t i) {
    // checking if the radius is a non-negative value.
    if (elements_[i].radius < 0) {
      sprintf(painCave.errMsg,
	      "BeadModel::checkElement: There is a bead with a negative radius.\n"
	      "\tStarting from index 0, check bead (%lu).\n", i);
      painCave.isFatal = 1;
      simError();
    }
    // if the bead's radius is below 1.0e-14, substitute by 1.0e-14;
    // to avoid problem in the self-interaction part (i.e., to not divide by
    // zero)
    if (elements_[i].radius < 1.0e-14) elements_[i].radius = 1.0e-14;
  }


  RealType BeadModel::volumeCorrection() {
    
    RealType volume(0.0); 
    for (std::vector<HydrodynamicsElement>::iterator iter = elements_.begin();
         iter != elements_.end(); ++iter) {
      volume += 4.0 / 3.0 * Constants::PI * pow((*iter).radius, 3);
    }
    // double loop double counts overlap volume Vij = Vji
    volume -= 0.5 * volumeOverlap_;

    return volume;
  }

  void BeadModel::writeElements(std::ostream& os) {
    std::vector<HydrodynamicsElement>::iterator iter;
    os << elements_.size() << std::endl;
    os << "Generated by Hydro" << std::endl;
    for (iter = elements_.begin(); iter != elements_.end(); ++iter) {
      os << iter->name << "\t" << iter->pos[0] << "\t" << iter->pos[1]
         << "\t" << iter->pos[2] << std::endl;
    }
  }
  Mat3x3d BeadModel::interactionTensor(const std::size_t i, const std::size_t j,
				       const RealType viscosity) {

    Mat3x3d Tij;
    Mat3x3d I = SquareMatrix3<RealType>::identity();
    
    if (i == j) {

      // self interaction, there is no overlapping volume
      RealType c = 1.0 / (6.0 * Constants::PI * viscosity * elements_[i].radius);
      
      Tij(0, 0) = c;
      Tij(1, 1) = c;
      Tij(2, 2) = c;

    } else {      

      // non-self interaction: divided into overlapping and
      // non-overlapping beads; the transitions between these cases
      // are continuous
      
      Vector3d Rij  = elements_[i].pos - elements_[j].pos;
      RealType rij  = Rij.length();
      RealType rij2 = rij * rij;
      
      if (rij >= (elements_[i].radius + elements_[j].radius)) {

	// non-overlapping beads
	RealType sumSigma2OverRij2 = ((elements_[i].radius * elements_[i].radius) +
				      (elements_[j].radius * elements_[j].radius)) / rij2;
	Mat3x3d tmpMat;
	tmpMat = outProduct(Rij, Rij) / rij2;
	RealType constant = 8.0 * Constants::PI * viscosity * rij;
	RealType tmp1     = 1.0 + sumSigma2OverRij2 / 3.0;
	RealType tmp2     = 1.0 - sumSigma2OverRij2;
	Tij               = (tmp1 * I + tmp2 * tmpMat) / constant;
	
      } else if (rij > fabs(elements_[i].radius - elements_[j].radius) &&
		 rij < (elements_[i].radius + elements_[j].radius)) {

	// overlapping beads, part I
	std::cout << "There is overlap between beads: (" << i
		  << ") and (" << j << ")" << std::endl;
	
	RealType sum_sigma       = (elements_[i].radius + elements_[j].radius);
	RealType subtr_sigma     = (elements_[i].radius - elements_[j].radius);
	RealType subtr_sigma_sqr = subtr_sigma * subtr_sigma;
	
	RealType constant = 6.0 * Constants::PI * viscosity *
	  (elements_[i].radius * elements_[j].radius);
	
	RealType rij3 = rij2 * rij;	
	Mat3x3d tmpMat;
	tmpMat = outProduct(Rij, Rij) / rij2;

	RealType tmp_var1 = subtr_sigma_sqr + 3.0 * rij2;
	RealType tmp1overlap =
	  (16.0 * rij3 * sum_sigma - tmp_var1 * tmp_var1) / (32.0 * rij3);
	
	RealType tmp_var2    = subtr_sigma_sqr - rij2;
	RealType tmp2overlap = (3.0 * tmp_var2 * tmp_var2) / (32.0 * rij3);
	
	Tij = (tmp1overlap * I + tmp2overlap * tmpMat) / constant;
	
	RealType volumetmp1 = (-rij + sum_sigma) * (-rij + sum_sigma);
	RealType volumetmp2 = (rij2 + 2.0 * (rij * elements_[i].radius) -
			       3.0 * (elements_[i].radius * elements_[i].radius) +
			       2 * (rij * elements_[j].radius) +
			       6 * (elements_[i].radius * elements_[j].radius) -
			       3 * (elements_[j].radius * elements_[j].radius));
	volumeOverlap_ += (Constants::PI / (12.0 * rij)) * volumetmp1 * volumetmp2;
	
      } else {

	// overlapping beads, part II: one bead inside the other

	if ((elements_[i].radius >= elements_[j].radius)) {

	  std::cout << "Bead: (" << j << ") is inside bead (" << i << ")"
		    << std::endl;

	  RealType constant =
	    1.0 / (6.0 * Constants::PI * viscosity * elements_[i].radius);
	  Tij(0, 0) = constant;
	  Tij(1, 1) = constant;
	  Tij(2, 2) = constant;

	  volumeOverlap_ += (4.0 / 3.0) * Constants::PI * pow(elements_[j].radius, 3);
	  
	} else {
	  
	  std::cout << "Bead: (" << i << ") is inside bead (" << j << ")"
		    << std::endl;

	  RealType constant =
	    1.0 / (6.0 * Constants::PI * viscosity * elements_[j].radius);
	  Tij(0, 0) = constant;
	  Tij(1, 1) = constant;
	  Tij(2, 2) = constant;
	  	  
	  volumeOverlap_ += (4.0 / 3.0) * Constants::PI * pow(elements_[i].radius, 3);
	}
      }
    }

    return Tij;
  }

}  // namespace OpenMD
