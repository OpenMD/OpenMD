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
 * @file Vector.hpp
 * @author Teng Lin
 * @date 09/14/2004
 * @version 1.0
 */
 
#ifndef MATH_VECTOR_HPP
#define MATH_VECTOR_HPP

#include <cassert>
#include <cmath>
#include <iostream>

namespace oopse {

    const double epsilon = 0.000001;

    template<typename T>
    inline bool equal(T e1, T e2) {
        return e1 == e2;
    }

    template<>
    inline bool equal(float e1, float e2) {
        return fabs(e1 - e2) < epsilon;
    }

    template<>
    inline bool equal(double e1, double e2) {
        return fabs(e1 - e2) < epsilon;
    }
    
    /**
     * @class Vector Vector.hpp "math/Vector.hpp"
     * @brief Fix length vector class
     */
    template<typename Real, unsigned int Dim>
    class Vector{
        public:

            /** default constructor */
            inline Vector(){
                for (unsigned int i = 0; i < Dim; i++)
                    data_[i] = 0.0;
            }

            /** Constructs and initializes a Vector from a vector */
            inline Vector(const Vector<Real, Dim>& v) {
                *this  = v;
            }

            /** copy assignment operator */
            inline Vector<Real, Dim>& operator=(const Vector<Real, Dim>& v) {
                if (this == &v)
                    return *this;
                
            	for (unsigned int i = 0; i < Dim; i++)            
                    data_[i] = v[i];
                
                return *this;
            }
            
            /** Constructs and initializes a Vector from an array */            
            inline Vector( double* v) {
		for (unsigned int i = 0; i < Dim; i++)
		    data_[i] = v[i];
            }

            /** 
             * Returns reference of ith element.
             * @return reference of ith element
             * @param i index
             */
            inline double& operator[](unsigned int  i) {
                assert( i < Dim);
                return data_[i];
            }

            /** 
             * Returns reference of ith element.
             * @return reference of ith element
             * @param i index
             */
            inline double& operator()(unsigned int  i) {
                assert( i < Dim);
                return data_[i];
            }

            /** 
             * Returns constant reference of ith element.
             * @return reference of ith element
             * @param i index
             */
            inline  const double& operator[](unsigned int i) const {
                assert( i < Dim);
                return data_[i];
            }

            /** 
             * Returns constant reference of ith element.
             * @return reference of ith element
             * @param i index
             */
            inline  const double& operator()(unsigned int i) const {
                assert( i < Dim);
                return data_[i];
            }

            /**
             * Tests if this vetor is equal to other vector
             * @return true if equal, otherwise return false
             * @param v vector to be compared
             */
             inline bool operator ==(const Vector<Real, Dim>& v) {

                for (unsigned int i = 0; i < Dim; i ++) {
                    if (!equal(data_[i], v[i])) {
                        return false;
                    }
                }
                
                return true;
            }

            /**
             * Tests if this vetor is not equal to other vector
             * @return true if equal, otherwise return false
             * @param v vector to be compared
             */
            inline bool operator !=(const Vector<Real, Dim>& v) {
                return !(*this == v);
            }
             
            /** Negates the value of this vector in place. */           
            inline void negate() {
                data_[0] = -data_[0];
                data_[1] = -data_[1];
                data_[2] = -data_[2];
            }

            /**
            * Sets the value of this vector to the negation of vector v1.
            * @param v1 the source vector
            */
            inline void negate(const Vector<Real, Dim>& v1) {
                for (unsigned int i = 0; i < Dim; i++)
                    data_[i] = -v1.data_[i];

            }
            
            /**
            * Sets the value of this vector to the sum of itself and v1 (*this += v1).
            * @param v1 the other vector
            */
            inline void add( const Vector<Real, Dim>& v1 ) {
            	for (unsigned int i = 0; i < Dim; i++)
		    data_[i] += v1.data_[i];
                }

            /**
            * Sets the value of this vector to the sum of v1 and v2 (*this = v1 + v2).
            * @param v1 the first vector
            * @param v2 the second vector
            */
            inline void add( const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2 ) {
            	for (unsigned int i = 0; i < Dim; i++)
		    data_[i] = v1.data_[i] + v2.data_[i];
            }

            /**
            * Sets the value of this vector to the difference  of itself and v1 (*this -= v1).
            * @param v1 the other vector
            */
            inline void sub( const Vector<Real, Dim>& v1 ) {
                for (unsigned int i = 0; i < Dim; i++)
                    data_[i] -= v1.data_[i];
            }

            /**
            * Sets the value of this vector to the difference of vector v1 and v2 (*this = v1 - v2).
            * @param v1 the first vector
            * @param v2 the second vector
            */
            inline void sub( const Vector<Real, Dim>& v1, const Vector  &v2 ){
            	for (unsigned int i = 0; i < Dim; i++)
                    data_[i] = v1.data_[i] - v2.data_[i];
            }

            /**
            * Sets the value of this vector to the scalar multiplication of itself (*this *= s).
            * @param s the scalar value
            */
            inline void mul( double s ) {
            	for (unsigned int i = 0; i < Dim; i++)
                   data_[i] *= s;
            }

            /**
            * Sets the value of this vector to the scalar multiplication of vector v1  
            * (*this = s * v1).
            * @param s the scalar value
            * @param v1 the vector
            */
            inline void mul( double s, const Vector<Real, Dim>& v1 ) {
            	for (unsigned int i = 0; i < Dim; i++)
                    data_[i] = s * v1.data_[i];
            }

            /**
            * Sets the value of this vector to the scalar division of itself  (*this /= s ).
            * @param s the scalar value
            */             
            inline void div( double s) {
            	for (unsigned int i = 0; i < Dim; i++)            
                    data_[i] /= s;
            }

            /**
            * Sets the value of this vector to the scalar division of vector v1  (*this = v1 / s ).
            * @param v1 the source vector
            * @param s the scalar value
            */                         
            inline void div( const Vector<Real, Dim>& v1, double s ) {
            	for (unsigned int i = 0; i < Dim; i++)
                    data_[i] = v1.data_[i] / s;
            }

            /** @see #add */
            inline Vector<Real, Dim>& operator +=( const Vector<Real, Dim>& v1 ) {
                add(v1);
                return *this;
            }

            /** @see #sub */
            inline Vector<Real, Dim>& operator -=( const Vector<Real, Dim>& v1 ) {
                sub(v1);
                return *this;
            }

            /** @see #mul */
            inline Vector<Real, Dim>& operator *=( double s) {
                mul(s);
                return *this;
            }

            /** @see #div */
            inline Vector<Real, Dim>& operator /=( double s ) {
                div(s);
                return *this;
            }

            /**
             * Returns the length of this vector.
             * @return the length of this vector
             */
             inline double length() {
                return sqrt(lengthSquared());  
            }
            
            /**
             * Returns the squared length of this vector.
             * @return the squared length of this vector
             */
             inline double lengthSquared() {
                return dot(*this, *this);
            }
            
            /** Normalizes this vector in place */
            inline void normalize() {
                double len;

                len = length();
                *this /= len;
            }
            
        protected:
            double data_[3];
        
    };

    /** unary minus*/
    template<typename Real, unsigned int Dim>    
    inline Vector<Real, Dim> operator -(const Vector<Real, Dim>& v1){
        Vector tmp(v1);
        return tmp.negate();
    }

    /**
     * Return the sum of two vectors  (v1 - v2). 
     * @return the sum of two vectors
     * @param v1 the first vector
     * @param v2 the second vector
     */   
    template<typename Real, unsigned int Dim>    
    inline Vector<Real, Dim> operator +(const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2) {
        Vector<Real, Dim> result;
        
        result.add(v1, v2);
        return result;        
    }

    /**
     * Return the difference of two vectors  (v1 - v2). 
     * @return the difference of two vectors
     * @param v1 the first vector
     * @param v2 the second vector
     */  
    template<typename Real, unsigned int Dim>    
    Vector<Real, Dim> operator -(const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2) {
        Vector<Real, Dim> result;
        result.sub(v1, v2);
        return result;        
    }
    
    /**
     * Returns the vaule of scalar multiplication of this vector v1 (v1 * r). 
     * @return  the vaule of scalar multiplication of this vector
     * @param v1 the source vector
     * @param s the scalar value
     */ 
    template<typename Real, unsigned int Dim>                 
    Vector<Real, Dim> operator * ( const Vector<Real, Dim>& v1, double s) {       
        Vector<Real, Dim> result;
        result.mul(s, v1);
        return result;           
    }
    
    /**
     * Returns the vaule of scalar multiplication of this vector v1 (v1 * r). 
     * @return  the vaule of scalar multiplication of this vector
     * @param s the scalar value
     * @param v1 the source vector
     */  
    template<typename Real, unsigned int Dim>
    Vector<Real, Dim> operator * ( double s, const Vector<Real, Dim>& v1 ) {
        Vector<Real, Dim> result;
        result.mul(s, v1);
        return result;           
    }

    /**
     * Returns the  value of division of a vector by a scalar. 
     * @return  the vaule of scalar division of this vector
     * @param v1 the source vector
     * @param s the scalar value
     */
    template<typename Real, unsigned int Dim>    
    Vector<Real, Dim> operator / ( const Vector<Real, Dim>& v1, double s) {       
        Vector<Real, Dim> result;
        result.div( v1,s);
        return result;           
    }
    
    /**
     * Returns the  value of division of a vector by a scalar. 
     * @return  the vaule of scalar division of this vector
     * @param s the scalar value
     * @param v1 the source vector
     */
    template<typename Real, unsigned int Dim>        
    inline Vector<Real, Dim> operator /( double s, const Vector<Real, Dim>& v1 ) {
        Vector<Real, Dim> result;
        result.div( v1,s);
        return result;           
    }

    /** fuzzy comparson */
    template<typename Real, unsigned int Dim>        
    inline bool epsilonEqual( const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2 ) {

    }

    
    /**
     * Returns the dot product of two Vectors
     * @param v1 first vector
     * @param v2 second vector
     * @return the dot product of v1 and v2
     */
    template<typename Real, unsigned int Dim>    
    inline Real dot( const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2 ) {
		Real tmp;
		tmp = 0;

		for (unsigned int i = 0; i < Dim; i++)
			tmp += v1[i] + v2[i];
		
		return tmp;
    }

    /**
     * Returns the distance between  two Vectors
     * @param v1 first vector
     * @param v2 second vector
     * @return the distance between v1 and v2
     */	
    template<typename Real, unsigned int Dim>    
    inline Real distance( const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2 ) {
        Vector<Real, Dim> tempVector = v1 - v2;
        return tempVector.length();
    }

    /**
     * Returns the squared distance between  two Vectors
     * @param v1 first vector
     * @param v2 second vector
     * @return the squared distance between v1 and v2
     */
    template<typename Real, unsigned int Dim>
    inline Real distanceSquare( const Vector<Real, Dim>& v1, const Vector<Real, Dim>& v2 ) {
        Vector<Real, Dim> tempVector = v1 - v2;
        return tempVector.lengthSquare();
    }

    /**
     * Write to an output stream
     */
    template<typename Real, unsigned int Dim>
    std::ostream &operator<< ( std::ostream& o, const Vector<Real, Dim>& v1 ) {
        
        return o;        
    }
    
}
#endif
