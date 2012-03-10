/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
/**
 * @file Quaternion.hpp
 * @author Teng Lin
 * @date 10/11/2004
 * @version 1.0
 */

#ifndef MATH_QUATERNION_HPP
#define MATH_QUATERNION_HPP

#include "math/Vector3.hpp"
#include "math/SquareMatrix.hpp"
#define ISZERO(a,eps) ( (a)>-(eps) && (a)<(eps) )
const RealType tiny=1.0e-6;     

namespace OpenMD{

  /**
   * @class Quaternion Quaternion.hpp "math/Quaternion.hpp"
   * Quaternion is a sort of a higher-level complex number.
   * It is defined as Q = w + x*i + y*j + z*k,
   * where w, x, y, and z are numbers of type T (e.g. RealType), and
   * i*i = -1; j*j = -1; k*k = -1;
   * i*j = k; j*k = i; k*i = j;
   */
  template<typename Real>
  class Quaternion : public Vector<Real, 4> {

  public:
    Quaternion() : Vector<Real, 4>() {}

    /** Constructs and initializes a Quaternion from w, x, y, z values */     
    Quaternion(Real w, Real x, Real y, Real z) {
      this->data_[0] = w;
      this->data_[1] = x;
      this->data_[2] = y;
      this->data_[3] = z;                
    }
            
    /** Constructs and initializes a Quaternion from a  Vector<Real,4> */     
    Quaternion(const Vector<Real,4>& v) 
      : Vector<Real, 4>(v){
    }

    /** copy assignment */
    Quaternion& operator =(const Vector<Real, 4>& v){
      if (this == & v)
	return *this;
      
      Vector<Real, 4>::operator=(v);
      
      return *this;
    }
    
    /**
     * Returns the value of the first element of this quaternion.
     * @return the value of the first element of this quaternion
     */
    Real w() const {
      return this->data_[0];
    }

    /**
     * Returns the reference of the first element of this quaternion.
     * @return the reference of the first element of this quaternion
     */
    Real& w() {
      return this->data_[0];    
    }

    /**
     * Returns the value of the first element of this quaternion.
     * @return the value of the first element of this quaternion
     */
    Real x() const {
      return this->data_[1];
    }

    /**
     * Returns the reference of the second element of this quaternion.
     * @return the reference of the second element of this quaternion
     */
    Real& x() {
      return this->data_[1];    
    }

    /**
     * Returns the value of the thirf element of this quaternion.
     * @return the value of the third element of this quaternion
     */
    Real y() const {
      return this->data_[2];
    }

    /**
     * Returns the reference of the third element of this quaternion.
     * @return the reference of the third element of this quaternion
     */           
    Real& y() {
      return this->data_[2];    
    }

    /**
     * Returns the value of the fourth element of this quaternion.
     * @return the value of the fourth element of this quaternion
     */
    Real z() const {
      return this->data_[3];
    }
    /**
     * Returns the reference of the fourth element of this quaternion.
     * @return the reference of the fourth element of this quaternion
     */
    Real& z() {
      return this->data_[3];    
    }

    /**
     * Tests if this quaternion is equal to other quaternion
     * @return true if equal, otherwise return false
     * @param q quaternion to be compared
     */
    inline bool operator ==(const Quaternion<Real>& q) {

      for (unsigned int i = 0; i < 4; i ++) {
	if (!equal(this->data_[i], q[i])) {
	  return false;
	}
      }
                
      return true;
    }
            
    /**
     * Returns the inverse of this quaternion
     * @return inverse
     * @note since quaternion is a complex number, the inverse of quaternion
     * q = w + xi + yj+ zk is inv_q = (w -xi - yj - zk)/(|q|^2)
     */
    Quaternion<Real> inverse() {
      Quaternion<Real> q;
      Real d = this->lengthSquare();
                
      q.w() = w() / d;
      q.x() = -x() / d;
      q.y() = -y() / d;
      q.z() = -z() / d;
                
      return q;
    }

    /**
     * Sets the value to the multiplication of itself and another quaternion
     * @param q the other quaternion
     */
    void mul(const Quaternion<Real>& q) {
      Quaternion<Real> tmp(*this);

      this->data_[0] = (tmp[0]*q[0]) -(tmp[1]*q[1]) - (tmp[2]*q[2]) - (tmp[3]*q[3]);
      this->data_[1] = (tmp[0]*q[1]) + (tmp[1]*q[0]) + (tmp[2]*q[3]) - (tmp[3]*q[2]);
      this->data_[2] = (tmp[0]*q[2]) + (tmp[2]*q[0]) + (tmp[3]*q[1]) - (tmp[1]*q[3]);
      this->data_[3] = (tmp[0]*q[3]) + (tmp[3]*q[0]) + (tmp[1]*q[2]) - (tmp[2]*q[1]);                
    }

    void mul(const Real& s) {
      this->data_[0] *= s;
      this->data_[1] *= s;
      this->data_[2] *= s;
      this->data_[3] *= s;
    }

    /** Set the value of this quaternion to the division of itself by another quaternion */
    void div(Quaternion<Real>& q) {
      mul(q.inverse());
    }

    void div(const Real& s) {
      this->data_[0] /= s;
      this->data_[1] /= s;
      this->data_[2] /= s;
      this->data_[3] /= s;
    }
            
    Quaternion<Real>& operator *=(const Quaternion<Real>& q) {
      mul(q);
      return *this;
    }

    Quaternion<Real>& operator *=(const Real& s) {
      mul(s);
      return *this;
    }
            
    Quaternion<Real>& operator /=(Quaternion<Real>& q) {                
      *this *= q.inverse();
      return *this;
    }

    Quaternion<Real>& operator /=(const Real& s) {
      div(s);
      return *this;
    }            
    /**
     * Returns the conjugate quaternion of this quaternion
     * @return the conjugate quaternion of this quaternion
     */
    Quaternion<Real> conjugate() const {
      return Quaternion<Real>(w(), -x(), -y(), -z());            
    }


    /**
       return rotation angle from -PI to PI 
    */
    inline Real get_rotation_angle() const{
      if( w() < (Real)0.0 )
        return 2.0*atan2(-sqrt( x()*x() + y()*y() + z()*z() ), -w() );
      else
        return 2.0*atan2( sqrt( x()*x() + y()*y() + z()*z() ),  w() );
    }

    /**
       create a unit quaternion from axis angle representation
    */
    Quaternion<Real> fromAxisAngle(const Vector3<Real>& axis, 
                                   const Real& angle){
      Vector3<Real> v(axis);
      v.normalize();
      Real half_angle = angle*0.5;
      Real sin_a = sin(half_angle);
      *this = Quaternion<Real>(cos(half_angle), 
                               v.x()*sin_a, 
                               v.y()*sin_a, 
                               v.z()*sin_a);
      return *this;
    }
    
    /**
       convert a quaternion to axis angle representation, 
       preserve the axis direction and angle from -PI to +PI
    */
    void toAxisAngle(Vector3<Real>& axis, Real& angle)const {
      Real vl = sqrt( x()*x() + y()*y() + z()*z() );
      if( vl > tiny ) {
        Real ivl = 1.0/vl;
        axis.x() = x() * ivl;
        axis.y() = y() * ivl;
        axis.z() = z() * ivl;

        if( w() < 0 )
          angle = 2.0*atan2(-vl, -w()); //-PI,0 
        else
          angle = 2.0*atan2( vl,  w()); //0,PI 
      } else {
        axis = Vector3<Real>(0.0,0.0,0.0);
        angle = 0.0;
      }
    }

    /**
       shortest arc quaternion rotate one vector to another by shortest path.
       create rotation from -> to, for any length vectors.
    */
    Quaternion<Real> fromShortestArc(const Vector3d& from, 
                                     const Vector3d& to ) {
      
      Vector3d c( cross(from,to) );
      *this = Quaternion<Real>(dot(from,to), 
                               c.x(), 
                               c.y(),
                               c.z());

      this->normalize();    // if "from" or "to" not unit, normalize quat
      w() += 1.0f;            // reducing angle to halfangle
      if( w() <= 1e-6 ) {     // angle close to PI
        if( ( from.z()*from.z() ) > ( from.x()*from.x() ) ) {
          this->data_[0] =  w();    
          this->data_[1] =  0.0;       //cross(from , Vector3d(1,0,0))
          this->data_[2] =  from.z();
          this->data_[3] = -from.y();
        } else {
          this->data_[0] =  w();
          this->data_[1] =  from.y();  //cross(from, Vector3d(0,0,1))
          this->data_[2] = -from.x();
          this->data_[3] =  0.0;
        }
      }
      this->normalize(); 
    }

    Real ComputeTwist(const Quaternion& q) {
      return  (Real)2.0 * atan2(q.z(), q.w());
    }

    void RemoveTwist(Quaternion& q) {
      Real t = ComputeTwist(q);
      Quaternion rt = fromAxisAngle(V3Z, t);
      
      q *= rt.inverse();
    }

    void getTwistSwingAxisAngle(Real& twistAngle, Real& swingAngle, 
                                Vector3<Real>& swingAxis) {
      
      twistAngle = (Real)2.0 * atan2(z(), w());
      Quaternion rt, rs;
      rt.fromAxisAngle(V3Z, twistAngle);
      rs = *this * rt.inverse();
      
      Real vl = sqrt( rs.x()*rs.x() + rs.y()*rs.y() + rs.z()*rs.z() );
      if( vl > tiny ) {
        Real ivl = 1.0 / vl;
        swingAxis.x() = rs.x() * ivl;
        swingAxis.y() = rs.y() * ivl;
        swingAxis.z() = rs.z() * ivl;

        if( rs.w() < 0.0 )
          swingAngle = 2.0*atan2(-vl, -rs.w()); //-PI,0 
        else
          swingAngle = 2.0*atan2( vl,  rs.w()); //0,PI 
      } else {
        swingAxis = Vector3<Real>(1.0,0.0,0.0);
        swingAngle = 0.0;
      }           
    }


    Vector3<Real> rotate(const Vector3<Real>& v) {

      Quaternion<Real> q(v.x() * w() + v.z() * y() - v.y() * z(),
                         v.y() * w() + v.x() * z() - v.z() * x(),
                         v.z() * w() + v.y() * x() - v.x() * y(),
                         v.x() * x() + v.y() * y() + v.z() * z());

      return Vector3<Real>(w()*q.x() + x()*q.w() + y()*q.z() - z()*q.y(),
                           w()*q.y() + y()*q.w() + z()*q.x() - x()*q.z(),
                           w()*q.z() + z()*q.w() + x()*q.y() - y()*q.x())*
        ( 1.0/this->lengthSquare() );      
    }   

    Quaternion<Real>& align (const Vector3<Real>& V1,
                             const Vector3<Real>& V2) {

      // If V1 and V2 are not parallel, the axis of rotation is the unit-length
      // vector U = Cross(V1,V2)/Length(Cross(V1,V2)).  The angle of rotation,
      // A, is the angle between V1 and V2.  The quaternion for the rotation is
      // q = cos(A/2) + sin(A/2)*(ux*i+uy*j+uz*k) where U = (ux,uy,uz).
      //
      // (1) Rather than extract A = acos(Dot(V1,V2)), multiply by 1/2, then
      //     compute sin(A/2) and cos(A/2), we reduce the computational costs
      //     by computing the bisector B = (V1+V2)/Length(V1+V2), so cos(A/2) =
      //     Dot(V1,B).
      //
      // (2) The rotation axis is U = Cross(V1,B)/Length(Cross(V1,B)), but
      //     Length(Cross(V1,B)) = Length(V1)*Length(B)*sin(A/2) = sin(A/2), in
      //     which case sin(A/2)*(ux*i+uy*j+uz*k) = (cx*i+cy*j+cz*k) where
      //     C = Cross(V1,B).
      //
      // If V1 = V2, then B = V1, cos(A/2) = 1, and U = (0,0,0).  If V1 = -V2,
      // then B = 0.  This can happen even if V1 is approximately -V2 using
      // floating point arithmetic, since Vector3::Normalize checks for
      // closeness to zero and returns the zero vector accordingly.  The test
      // for exactly zero is usually not recommend for floating point
      // arithmetic, but the implementation of Vector3::Normalize guarantees
      // the comparison is robust.  In this case, the A = pi and any axis
      // perpendicular to V1 may be used as the rotation axis.

      Vector3<Real> Bisector = V1 + V2;
      Bisector.normalize();

      Real CosHalfAngle = dot(V1,Bisector);

      this->data_[0] = CosHalfAngle;

      if (CosHalfAngle != (Real)0.0) {
        Vector3<Real> Cross = cross(V1, Bisector);
        this->data_[1] = Cross.x();
        this->data_[2] = Cross.y();
        this->data_[3] = Cross.z();
      } else {
        Real InvLength;
        if (fabs(V1[0]) >= fabs(V1[1])) {
          // V1.x or V1.z is the largest magnitude component
          InvLength = (Real)1.0/sqrt(V1[0]*V1[0] + V1[2]*V1[2]);

          this->data_[1] = -V1[2]*InvLength;
          this->data_[2] = (Real)0.0;
          this->data_[3] = +V1[0]*InvLength;
        } else {
          // V1.y or V1.z is the largest magnitude component
          InvLength = (Real)1.0 / sqrt(V1[1]*V1[1] + V1[2]*V1[2]);
          
          this->data_[1] = (Real)0.0;
          this->data_[2] = +V1[2]*InvLength;
          this->data_[3] = -V1[1]*InvLength;
        }
      }
      return *this;
    }

    void toTwistSwing ( Real& tw, Real& sx, Real& sy ) {
      
      // First test if the swing is in the singularity:

      if ( ISZERO(z(),tiny) && ISZERO(w(),tiny) ) { sx=sy=M_PI; tw=0; return; }

      // Decompose into twist-swing by solving the equation:
      //
      //                       Qtwist(t*2) * Qswing(s*2) = q
      //
      // note: (x,y) is the normalized swing axis (x*x+y*y=1)
      //
      //          ( Ct 0 0 St ) * ( Cs xSs ySs 0 ) = ( qw qx qy qz )
      //  ( CtCs  xSsCt-yStSs  xStSs+ySsCt  StCs ) = ( qw qx qy qz )      (1)
      // From (1): CtCs / StCs = qw/qz => Ct/St = qw/qz => tan(t) = qz/qw (2)
      //
      // The swing rotation/2 s comes from:
      //
      // From (1): (CtCs)^2 + (StCs)^2 = qw^2 + qz^2 =>  
      //                                       Cs = sqrt ( qw^2 + qz^2 ) (3)
      //
      // From (1): (xSsCt-yStSs)^2 + (xStSs+ySsCt)^2 = qx^2 + qy^2 => 
      //                                       Ss = sqrt ( qx^2 + qy^2 ) (4)
      // From (1):  |SsCt -StSs| |x| = |qx|
      //            |StSs +SsCt| |y|   |qy|                              (5)

      Real qw, qx, qy, qz;
      
      if ( w()<0 ) {
        qw=-w(); 
        qx=-x(); 
        qy=-y(); 
        qz=-z();
      } else {
        qw=w(); 
        qx=x(); 
        qy=y(); 
        qz=z();
      }
      
      Real t = atan2 ( qz, qw ); // from (2)
      Real s = atan2( sqrt(qx*qx+qy*qy), sqrt(qz*qz+qw*qw) ); // from (3)
                                                              // and (4)

      Real x=0.0, y=0.0, sins=sin(s);

      if ( !ISZERO(sins,tiny) ) {
        Real sint = sin(t);
        Real cost = cos(t);
        
        // by solving the linear system in (5):
        y = (-qx*sint + qy*cost)/sins;
        x = ( qx*cost + qy*sint)/sins;
      }

      tw = (Real)2.0*t;
      sx = (Real)2.0*x*s;
      sy = (Real)2.0*y*s;
    }

    void toSwingTwist(Real& sx, Real& sy, Real& tw ) {

      // Decompose q into swing-twist using a similar development as
      // in function toTwistSwing

      if ( ISZERO(z(),tiny) && ISZERO(w(),tiny) ) { sx=sy=M_PI; tw=0; }
      
      Real qw, qx, qy, qz;
      if ( w() < 0 ){
        qw=-w(); 
        qx=-x(); 
        qy=-y(); 
        qz=-z();
      } else {
        qw=w(); 
        qx=x(); 
        qy=y(); 
        qz=z(); 
      }

      // Get the twist t:
      Real t = 2.0 * atan2(qz,qw);
      
      Real bet = atan2( sqrt(qx*qx+qy*qy), sqrt(qz*qz+qw*qw) );
      Real gam = t/2.0;
      Real sinc = ISZERO(bet,tiny)? 1.0 : sin(bet)/bet;
      Real singam = sin(gam);
      Real cosgam = cos(gam);

      sx = Real( (2.0/sinc) * (cosgam*qx - singam*qy) );
      sy = Real( (2.0/sinc) * (singam*qx + cosgam*qy) );
      tw = Real( t );
    }
      
    
    
    /**
     * Returns the corresponding rotation matrix (3x3)
     * @return a 3x3 rotation matrix
     */
    SquareMatrix<Real, 3> toRotationMatrix3() {
      SquareMatrix<Real, 3> rotMat3;
      
      Real w2;
      Real x2;
      Real y2;
      Real z2;

      if (!this->isNormalized())
	this->normalize();
                
      w2 = w() * w();
      x2 = x() * x();
      y2 = y() * y();
      z2 = z() * z();

      rotMat3(0, 0) = w2 + x2 - y2 - z2;
      rotMat3(0, 1) = 2.0 * ( x() * y() + w() * z() );
      rotMat3(0, 2) = 2.0 * ( x() * z() - w() * y() );

      rotMat3(1, 0) = 2.0 * ( x() * y() - w() * z() );
      rotMat3(1, 1) = w2 - x2 + y2 - z2;
      rotMat3(1, 2) = 2.0 * ( y() * z() + w() * x() );

      rotMat3(2, 0) = 2.0 * ( x() * z() + w() * y() );
      rotMat3(2, 1) = 2.0 * ( y() * z() - w() * x() );
      rotMat3(2, 2) = w2 - x2 -y2 +z2;

      return rotMat3;
    }

  };//end Quaternion


    /**
     * Returns the vaule of scalar multiplication of this quaterion q (q * s). 
     * @return  the vaule of scalar multiplication of this vector
     * @param q the source quaternion
     * @param s the scalar value
     */
  template<typename Real, unsigned int Dim>                 
  Quaternion<Real> operator * ( const Quaternion<Real>& q, Real s) {       
    Quaternion<Real> result(q);
    result.mul(s);
    return result;           
  }
    
  /**
   * Returns the vaule of scalar multiplication of this quaterion q (q * s). 
   * @return  the vaule of scalar multiplication of this vector
   * @param s the scalar value
   * @param q the source quaternion
   */  
  template<typename Real, unsigned int Dim>
  Quaternion<Real> operator * ( const Real& s, const Quaternion<Real>& q ) {
    Quaternion<Real> result(q);
    result.mul(s);
    return result;           
  }    

  /**
   * Returns the multiplication of two quaternion
   * @return the multiplication of two quaternion
   * @param q1 the first quaternion
   * @param q2 the second quaternion
   */
  template<typename Real>
  inline Quaternion<Real> operator *(const Quaternion<Real>& q1, const Quaternion<Real>& q2) {
    Quaternion<Real> result(q1);
    result *= q2;
    return result;
  }

  /**
   * Returns the division of two quaternion
   * @param q1 divisor
   * @param q2 dividen
   */

  template<typename Real>
  inline Quaternion<Real> operator /( Quaternion<Real>& q1,  Quaternion<Real>& q2) {
    return q1 * q2.inverse();
  }

  /**
   * Returns the value of the division of a scalar by a quaternion
   * @return the value of the division of a scalar by a quaternion
   * @param s scalar
   * @param q quaternion
   * @note for a quaternion q, 1/q = q.inverse()
   */
  template<typename Real>
  Quaternion<Real> operator /(const Real& s,  Quaternion<Real>& q) {

    Quaternion<Real> x;
    x = q.inverse();
    x *= s;
    return x;
  }
    
  template <class T>
  inline bool operator==(const Quaternion<T>& lhs, const Quaternion<T>& rhs) {
    return equal(lhs[0] ,rhs[0]) && equal(lhs[1] , rhs[1]) && equal(lhs[2], rhs[2]) && equal(lhs[3], rhs[3]);
  }
    
  typedef Quaternion<RealType> Quat4d;
}
#endif //MATH_QUATERNION_HPP 
