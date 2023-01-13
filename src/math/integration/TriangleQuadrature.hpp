// Note that this code is mostly from the Drake project: https://drake.mit.edu/
// and is governed by their BSD 3-Clause license:
// https://github.com/RobotLocomotion/drake/blob/master/LICENSE.TXT
// Only minor modifications have been made to make it work in OpenMD.

#ifndef MATH_INTEGRATION_TRIANGLEQUADRATURE_HPP
#define MATH_INTEGRATION_TRIANGLEQUADRATURE_HPP

#include <functional>
#include <vector>

#include "math/Vector2.hpp"
#include "math/Vector3.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/integration/TriangleQuadratureRule.hpp"

namespace OpenMD {

  /// A class for integrating a function using numerical quadrature
  /// over triangular domains.
  ///
  /// @tparam NumericReturnType the output type of the function being
  ///         integrated.  Commonly will be a IEEE floating point
  ///         scalar (e.g., `double`), but could be a
  ///         VectorXd, or any other
  ///         numeric type that supports both scalar multiplication
  ///         (i.e., operator*(const NumericReturnType&, double) and
  ///         addition with another of the same type (i.e.,
  ///         operator+(const NumericReturnType&, const
  ///         NumericReturnType&)).  
  /// @tparam T the scalar type of the function being integrated
  ///         over. Supported types are currently only IEEE floating
  ///         point scalars.
  template <typename NumericReturnType, typename T>
  class TriangleQuadrature {
  public:
    /// Numerically integrates the function f over a triangle using the given
    /// quadrature rule and the initial value.
    /// @param f(p) a function that returns a numerical value for point p in the
    ///        domain of the triangle, specified in barycentric coordinates.
    ///        The barycentric coordinates are given by
    ///        (p[0], p[1], 1 - p[0] - p[1]).
    /// @param area the area of the triangle.
    static NumericReturnType Integrate(
				       const std::function<NumericReturnType(const Vector<T,2>&)>& f,
				       const TriangleQuadratureRule& rule,
				       const T& area);

    /// Alternative signature for Integrate() that uses three-dimensional
    /// barycentric coordinates.
    /// @param f(p) a function that returns a numerical value for point p in the
    ///        domain of the triangle, specified in barycentric coordinates.
    ///        The barycentric coordinates are given by
    ///        (p[0], p[1], p[2]).
    /// @param area the area of the triangle.
    static NumericReturnType Integrate(
				       const std::function<NumericReturnType(const Vector<T,3>&)>& f,
				       const TriangleQuadratureRule& rule,
				       const T& area);
  };

  template <typename NumericReturnType, typename T>
  NumericReturnType TriangleQuadrature<NumericReturnType, T>::Integrate(
									const std::function<NumericReturnType(const Vector<T,2>&)>& f,
									const TriangleQuadratureRule& rule,
									const T& area) {
    // Get the quadrature points and weights.
    const std::vector<Vector2<T>>& barycentric_coordinates =
      rule.quadrature_points();
    const std::vector<T>& weights = rule.weights();
    assert(barycentric_coordinates.size() == weights.size());
    assert(weights.size() >= 1);

    // Sum the weighted function evaluated at the transformed
    // quadrature points.  The looping is done in this particular way
    // so that the return type can be either a traditional scalar
    // (e.g., `double`), a Vector, or some other numeric type,
    // without having to worry about the numerous possible ways to
    // initialize to zero (e.g., `double integral = 0.0` vs.
    // VectorXd = VectorXd:zZero())` or having to know
    // the dimension of the NumericReturnType (i.e., scalar or vector
    // dimension).
    NumericReturnType integral = f(barycentric_coordinates[0]) * weights[0];
    for (int i = 1; i < static_cast<int>(weights.size()); ++i)
      integral += f(barycentric_coordinates[i]) * weights[i];
    
    return integral * area;
  }
  
  template <typename NumericReturnType, typename T>
  NumericReturnType TriangleQuadrature<NumericReturnType, T>::Integrate(
									const std::function<NumericReturnType(const Vector<T,3>&)>& f,
									const TriangleQuadratureRule& rule,
									const T& area) {
    // Create a function fprime accepting a 2D barycentric representation but
    // using it to call the given 3D function f.
    std::function<NumericReturnType(const Vector<T,2>&)> fprime = [&f](
								       const Vector<T,2>& barycentric_2D) {
      return f(Vector3<T>(
			   barycentric_2D[0],
			   barycentric_2D[1],
			   T(1.0) - barycentric_2D[0] - barycentric_2D[1]));
    };

    return Integrate(fprime, rule, area);
  }
}
#endif
