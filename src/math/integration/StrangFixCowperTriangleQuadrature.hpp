#ifndef MATH_INTEGRATION_STRANGFIXCOWPERTRIANGLEQUADRATURE_HPP
#define MATH_INTEGRATION_STRANGFIXCOWPERTRIANGLEQUADRATURE_HPP

#include <stdexcept>
#include <vector>

#include "math/integration/TriangleQuadratureRule.hpp"

namespace OpenMD {
  
  class StrangFixCowperTriangleQuadratureRule final : public TriangleQuadratureRule {
  public:
    /// Constructs the StrangFixCowper quadrature rule of the specified order, which
    /// must be between 1 and 5.
    explicit StrangFixCowperTriangleQuadratureRule(int order) : order_(order) {
      assert(order >= 1);
      if (order > 7) {
	snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
		 "StrangFixCowperTriangleQuadrature does not implement a %d point algorithm\n",
		 order);
	painCave.isFatal = 1;;
	simError();
      }
      SetWeightsAndQuadraturePoints();
    }

  private:
    // A partial implementation of the Triangular quadrature schemes
    // in: Gilbert Strang & George Fix, An Analysis of the Finite
    // Element Method, (Wellesley-Cambridge Press, 1973),
    // https://bookstore.siam.org/wc08/
    //
    // and G.R. Cowper, "Gaussian quadrature formulas for triangles",
    // Numerical Methods in Engineering, 7(3), pp. 405-408 (1973).
    // https://doi.org/10.1002/nme.1620070316
    void SetWeightsAndQuadraturePoints() {
      switch (order_) {
      case 2:
	//Strang-Fix-Cowper scheme 2
	weights_.resize(3);
        quadrature_points_.resize(3);
        weights_[0] = weights_[1] = weights_[2] = 1.0/3.0;
        quadrature_points_[0]  = Vector2d(1.0/2.0, 1.0/2.0);
        quadrature_points_[1]  = Vector2d(0.0,     1.0/2.0);
        quadrature_points_[2]  = Vector2d(1.0/2.0, 0.0);
        break;
      case 3:
	//Strang-Fix-Cowper Scheme 3
        weights_.resize(4);
        quadrature_points_.resize(4);
        weights_[0] = -27.0/48.0;
        weights_[1] = weights_[2] = weights_[3] = 25.0/48.0;
        quadrature_points_[0]  = Vector2d(1.0/3.0, 1.0/3.0);
        quadrature_points_[1]  = Vector2d(1.0/5.0, 1.0/5.0);
        quadrature_points_[2]  = Vector2d(1.0/5.0, 3.0/5.0);
        quadrature_points_[3]  = Vector2d(3.0/5.0, 1.0/5.0);
        break;	
      case 4:
        weights_.resize(6);
        quadrature_points_.resize(6);
        weights_[0] = weights_[1] = weights_[2] = 0.223381589678011;
        weights_[3] = weights_[4] = weights_[5] = 0.109951743655322;
        quadrature_points_[0]  = Vector2d(0.445948490915965, 0.445948490915965);
        quadrature_points_[1]  = Vector2d(0.445948490915965, 0.108103018168070);
        quadrature_points_[2]  = Vector2d(0.108103018168070, 0.445948490915965);
        quadrature_points_[3]  = Vector2d(0.091576213509771, 0.091576213509771);
        quadrature_points_[4]  = Vector2d(0.091576213509771, 0.816847572980459);
        quadrature_points_[5]  = Vector2d(0.816847572980459, 0.091576213509771);
        break;	
      case 6:
	//Strang-Fix-Cowper scheme 4
	weights_.resize(6);
	quadrature_points_.resize(6);
	weights_[0] = weights_[1] = weights_[2] = 1.0/6.0;
	weights_[3] = weights_[4] = weights_[5] = 1.0/6.0;
	quadrature_points_[0]  = Vector2d(0.659027622374092, 0.231933368553031);
	quadrature_points_[1]  = Vector2d(0.109039009072877, 0.659027622374092);
	quadrature_points_[2]  = Vector2d(0.231933368553031, 0.109039009072877);
	quadrature_points_[3]  = Vector2d(0.231933368553031, 0.659027622374092);
	quadrature_points_[4]  = Vector2d(0.109039009072877, 0.231933368553031);
	quadrature_points_[5]  = Vector2d(0.659027622374092, 0.109039009072877);
	break;
      case 7:
	//Strang-Fix-Cowper scheme 7
	weights_.resize(7);
	quadrature_points_.resize(7);
	weights_[0] = 0.225;
	weights_[1] = weights_[2] = weights_[3] = 0.125939180544827;	
	weights_[4] = weights_[5] = weights_[6] = 0.132394152788506 ;
	quadrature_points_[0]  = Vector2d(1.0 / 3.0, 1.0 / 3.0);
	quadrature_points_[1]  = Vector2d(0.797426985353087, 0.101286507323456);
	quadrature_points_[2]  = Vector2d(0.101286507323456, 0.797426985353087);
	quadrature_points_[3]  = Vector2d(0.101286507323456, 0.101286507323456);
	quadrature_points_[4]  = Vector2d(0.059715871789770, 0.470142064105115);
	quadrature_points_[5]  = Vector2d(0.470142064105115, 0.059715871789770);
	quadrature_points_[6]  = Vector2d(0.470142064105115, 0.470142064105115);
	break;
	
      default:
	snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
		 "StrangFixCowperTriangleQuadrature does not implement a %d point algorithm\n",
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
