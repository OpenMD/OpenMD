/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
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
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "BeadModel::checkElement: There is a bead with a negative radius.\n"
          "\tStarting from index 0, check bead (%zu).\n",
          i);
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
      os << iter->name << "\t" << iter->pos[0] << "\t" << iter->pos[1] << "\t"
         << iter->pos[2] << std::endl;
    }
  }

  Mat3x3d BeadModel::interactionTensor(const std::size_t i, const std::size_t j,
                                       const RealType viscosity) {
    Mat3x3d Tij;
    Mat3x3d I = SquareMatrix3<RealType>::identity();
    RealType c {};

    if (i == j) {
      // self interaction, there is no overlapping volume
      c = 1.0 / (6.0 * Constants::PI * viscosity * elements_[i].radius);

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
        RealType a = ((elements_[i].radius * elements_[i].radius) +
                      (elements_[j].radius * elements_[j].radius)) /
                     rij2;
        Mat3x3d op;
        op          = outProduct(Rij, Rij) / rij2;
        c           = 1.0 / (8.0 * Constants::PI * viscosity * rij);
        RealType b1 = 1.0 + a / 3.0;
        RealType b2 = 1.0 - a;
        Tij         = (b1 * I + b2 * op) * c;

      } else if (rij > fabs(elements_[i].radius - elements_[j].radius) &&
                 rij < (elements_[i].radius + elements_[j].radius)) {
        // overlapping beads, part I

        RealType sum   = (elements_[i].radius + elements_[j].radius);
        RealType diff  = (elements_[i].radius - elements_[j].radius);
        RealType diff2 = diff * diff;

        c = 1.0 / (6.0 * Constants::PI * viscosity *
                   (elements_[i].radius * elements_[j].radius));

        RealType rij3 = rij2 * rij;
        Mat3x3d op;
        op = outProduct(Rij, Rij) / rij2;

        RealType a  = diff2 + 3.0 * rij2;
        RealType ao = (16.0 * rij3 * sum - a * a) / (32.0 * rij3);

        RealType b  = diff2 - rij2;
        RealType bo = (3.0 * b * b) / (32.0 * rij3);

        Tij = (ao * I + bo * op) * c;

        RealType v1 = (-rij + sum) * (-rij + sum);
        RealType v2 = (rij2 + 2.0 * (rij * elements_[i].radius) -
                       3.0 * (elements_[i].radius * elements_[i].radius) +
                       2.0 * (rij * elements_[j].radius) +
                       6.0 * (elements_[i].radius * elements_[j].radius) -
                       3.0 * (elements_[j].radius * elements_[j].radius));

        volumeOverlap_ += (Constants::PI / (12.0 * rij)) * v1 * v2;

      } else {
        // overlapping beads, part II: one bead inside the other

        RealType rmin = std::min(elements_[i].radius, elements_[j].radius);
        RealType rmax = std::max(elements_[i].radius, elements_[j].radius);

        c = 1.0 / (6.0 * Constants::PI * viscosity * rmax);

        Tij(0, 0) = c;
        Tij(1, 1) = c;
        Tij(2, 2) = c;

        volumeOverlap_ += (4.0 / 3.0) * Constants::PI * pow(rmin, 3);
      }
    }

    return Tij;
  }

}  // namespace OpenMD
