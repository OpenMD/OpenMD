// Note that this code is mostly from the Drake project: https://drake.mit.edu/
// and is governed by their BSD 3-Clause license:
// https://github.com/RobotLocomotion/drake/blob/master/LICENSE.TXT
// Only minor modifications have been made to make it work in OpenMD.

#ifndef MATH_INTEGRATION_TRIANGLEQUADRATURERULE_HPP
#define MATH_INTEGRATION_TRIANGLEQUADRATURERULE_HPP

#include <vector>

#include "math/Vector2.hpp"
#include "math/Vector3.hpp"

namespace OpenMD {

  /// A "rule" (weights and quadrature points) for computing quadrature over
  /// triangular domains.
  class TriangleQuadratureRule {
  public:
    TriangleQuadratureRule(const TriangleQuadratureRule&)            = default;
    TriangleQuadratureRule& operator=(const TriangleQuadratureRule&) = default;
    TriangleQuadratureRule(TriangleQuadratureRule&&)                 = default;
    TriangleQuadratureRule& operator=(TriangleQuadratureRule&&)      = default;
    TriangleQuadratureRule()                                         = default;
    virtual ~TriangleQuadratureRule() {}

    /// Returns the order of this rule.
    int order() const {
      int rule_order = do_order();
      assert(rule_order >= 1);
      return rule_order;
    }

    /// Returns the vector of quadrature points. These are returned as
    /// the first two barycentric coordinates b0 b1; the third is just
    /// b2 = 1 - b0 - b1.  Each of these has a corresponding weight
    /// returned by weights().
    const std::vector<Vector2d>& quadrature_points() const {
      return do_quadrature_points();
    }

    /// Returns the vector of weights. These sum to 1 and there is one
    /// weight for each point returned by quadrature_points().
    const std::vector<RealType>& weights() const { return do_weights(); }

  protected:
    /// Derived classes shall return the order (>= 1) of this rule.
    virtual int do_order() const = 0;

    /// Derived classes shall return the vector of quadrature
    /// points. Each of these Vector2d objects represents
    /// the barycentric coordinates of a triangle (the third
    /// barycentric coordinate is implicit: it is the difference
    /// between unity and the sum of the other two coordinates).
    virtual const std::vector<Vector2d>& do_quadrature_points() const = 0;

    /// Derived classes shall return the vector of weights. The sum of
    /// all weights must equal 1.0.
    virtual const std::vector<RealType>& do_weights() const = 0;
  };
}  // namespace OpenMD
#endif
