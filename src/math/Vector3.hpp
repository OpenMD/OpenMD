/*
 * Copyright (C) 2000-2004  Object Oriented Parallel Simulation Engine (OOPSE) project
 * 
 * Contact: oopse@oopse.org
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */

/**
 * @file Vector3.hpp
 * @author Teng Lin
 * @date 09/14/2004
 * @version 1.0
 */
 
#ifndef MATH_VECTOR3_HPP
#define MATH_VECTOR3_HPP

#include <cassert>
#include <cmath>

#include "Vector.hpp"

namespace oopse {
  
  /**
   * @class Vector3 Vector3.hpp "math/Vector3.hpp"
   * @brief 
   */
  
  template<typename Real>
  class Vector3 : public Vector<Real, 3>{
  public:

    Vector3() : Vector<Real, 3>(){}
    
    /** Constructs and initializes a Vector3 from x, y, z coordinates */
    inline Vector3( double x, double y, double z) {
      data_[0] = x;
      data_[1] = y;
      data_[2] = z;
    }
    
    inline Vector3(const Vector<Real, 3>& v) : Vector<Real, 3>(v) {}
    
    inline Vector3<Real>& operator = (const Vector<Real, 3>& v) {
      if (this == &v) { return *this; }
      Vector<Real, 3>::operator=(v);
      return *this;
    }
    
    /**
     * Retunrs reference of the first element of Vector3.
     * @return reference of the first element of Vector3
     */
    inline Real& x() {  return data_[0];}
    
    /**
     * Retunrs the first element of Vector3.
     * @return  the first element of Vector3
     */
    inline Real x() const {  return data_[0];}
    
    /**
     * Retunrs reference of the second element of Vector3.
     * @return reference of the second element of Vector3
     */
    inline Real& y() {  return data_[1];}
    
    /**
     * Retunrs  the second element of Vector3.
     * @return c the second element of Vector3
     */
    inline Real y() const {  return data_[1];}
    
    /**
     * Retunrs reference of the third element of Vector3.
     * @return reference of the third element of Vector3
     */
    inline Real& z() {  return data_[2];}
    
    /**
     * Retunrs  the third element of Vector3.
     * @return f the third element of Vector3
     */
    inline Real z() const {  return data_[2];}
    
  };
  
  /**
   * Returns the cross product of two Vectors
   * @param v1 first vector
   * @param v2 second vector
   * @return the cross product  of v1 and v2
   * @see #vector::dot
   */
  template<typename Real>
  Vector3<Real> cross( const Vector3<Real>& v1, const Vector3<Real>& v2 ) {
    Vector3<Real> result;
    
    result.x() = v1.y() * v2.z() - v1.z() * v2.y();
    result.y() = v1.z() * v2.x() - v1.x() * v2.z();
    result.z() = v1.x() * v2.y() - v1.y() * v2.x();
    
    return result;
  }
    
  typedef template Vector3<double> Vector3d;    
  
}

#endif
