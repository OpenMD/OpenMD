// Note that this code is mostly from the Drake project: https://drake.mit.edu/
// and is governed by their BSD 3-Clause license:
// https://github.com/RobotLocomotion/drake/blob/master/LICENSE.TXT
// Only minor modifications have been made to make it work in OpenMD.

#ifndef MATH_INTEGRATION_GAUSSIANTRIANGLEQUADRATURERULE_HPP
#define MATH_INTEGRATION_GAUSSIANTRIANGLEQUADRATURERULE_HPP

#include <stdexcept>
#include <vector>

#include "math/integration/TriangleQuadratureRule.hpp"

namespace OpenMD {
  
  class GaussianTriangleQuadratureRule final : public TriangleQuadratureRule {
  public:
    /// Constructs the Gaussian quadrature rule of the specified order, which
    /// must be between 1 and 5.
    explicit GaussianTriangleQuadratureRule(int order) : order_(order) {
      assert(order >= 1);
      if (order > 5) {
	snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
		 "GaussianTriangleQuadrature does not implement a %d point algorithm\n",
		 order);
	painCave.isFatal = 1;;
	simError();
      }
      SetWeightsAndQuadraturePoints();
    }

  private:
    // Sets the weights and quadrature points depending on the order of the
    // quadrature rule. Weights and quadrature points for the rational numbers
    // were taken from pp. 8-9 of:
    // http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
    // Weights and quadrature points for decimal numbers were taken from p. 1140
    // of:
    // D. A. Dunavant. High degree efficient symmetrical Gaussian quadrature
    // rules for the triangle. Intl. J. Num. Meth. Eng. pp. 1129-1148, 1985.
    void SetWeightsAndQuadraturePoints() {
      switch (order_) {
      case 1:
        weights_.resize(1);
        quadrature_points_.resize(1);
        weights_[0] = 1.0;
        quadrature_points_[0] = Vector2d(1.0/3.0, 1.0/3.0);
        break;
	
      case 2:
        weights_.resize(3);
        quadrature_points_.resize(3);
        weights_[0] = weights_[1] = weights_[2] = 1.0/3.0;
        quadrature_points_[0] = Vector2d(1.0/6.0, 1.0/6.0);
        quadrature_points_[1] = Vector2d(1.0/6.0, 2.0/3.0);
        quadrature_points_[2] = Vector2d(2.0/3.0, 1.0/6.0);
        break;
	
      case 3:
        weights_.resize(4);
        quadrature_points_.resize(4);
        weights_[0] = -27.0/48.0;
        weights_[1] = weights_[2] = weights_[3] = 25.0/48.0;
        quadrature_points_[0] = Vector2d(1.0/3.0, 1.0/3.0);
        quadrature_points_[1] = Vector2d(1.0/5.0, 1.0/5.0);
        quadrature_points_[2] = Vector2d(1.0/5.0, 3.0/5.0);
        quadrature_points_[3] = Vector2d(3.0/5.0, 1.0/5.0);
        break;
	
      case 4:
        weights_.resize(6);
        quadrature_points_.resize(6);
        weights_[0] = weights_[1] = weights_[2] = 0.223381589678011;
        weights_[3] = weights_[4] = weights_[5] = 0.109951743655322;
        quadrature_points_[0] = Vector2d(0.445948490915965, 0.445948490915965);
        quadrature_points_[1] = Vector2d(0.445948490915965, 0.108103018168070);
        quadrature_points_[2] = Vector2d(0.108103018168070, 0.445948490915965);
        quadrature_points_[3] = Vector2d(0.091576213509771, 0.091576213509771);
        quadrature_points_[4] = Vector2d(0.091576213509771, 0.816847572980459);
        quadrature_points_[5] = Vector2d(0.816847572980459, 0.091576213509771);
        break;
	
      case 5:
        weights_.resize(7);
        quadrature_points_.resize(7);
        weights_[0] = 0.225;
        weights_[1] = weights_[2] = weights_[3] = 0.132394152788506;
        weights_[4] = weights_[5] = weights_[6] = 0.125939180544827;
        quadrature_points_[0] = Vector2d(1.0/3.0, 1.0/3.0);
        quadrature_points_[1] = Vector2d(0.470142064105115, 0.470142064105115);
        quadrature_points_[2] = Vector2d(0.470142064105115, 0.059715871789770);
        quadrature_points_[3] = Vector2d(0.059715871789770, 0.470142064105115);
        quadrature_points_[4] = Vector2d(0.101286507323456, 0.101286507323456);
        quadrature_points_[5] = Vector2d(0.101286507323456, 0.797426985353087);
        quadrature_points_[6] = Vector2d(0.797426985353087, 0.101286507323456);
        break;
	
      default:
	snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
		 "GaussianTriangleQuadrature does not implement a %d point algorithm\n",
		 order_);
	painCave.isFatal = 1;;
	simError();
      }
    }
    
    int do_order() const final { return order_; }
    
    const std::vector<RealType>& do_weights() const final {
      return weights_;
    }
    
    const std::vector<Vector2d>& do_quadrature_points() const final {
      return quadrature_points_;
    }
    
    const int order_{-1};
    std::vector<RealType> weights_;
    std::vector<Vector2d> quadrature_points_;
  };
  
}
#endif
