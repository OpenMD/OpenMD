/*
 * Matpack Wigner3jm special function imported and modified for use in
 * OpenMD
 *                                                                              
 * Matpack Library Release 1.9.0                                                
 * Copyright (C) 1991-2003 by Berndt M. Gammel. All rights reserved.            
 *
 * Permission to use, copy, and distribute Matpack in its entirety
 * and its documentation for non-commercial purpose and without fee
 * is hereby granted, provided that this license information and
 * copyright notice appear unmodified in all copies.  This software
 * is provided 'as is' without express or implied warranty.  In no
 * event will the author be held liable for any damages arising from
 * the use of this software.
 *
 * Note that distributing Matpack 'bundled' in with any product is
 * considered to be a 'commercial purpose'.
 *
 * The software may be modified for your own purposes, but modified
 * versions may not be distributed without prior consent of the
 * author.
 *                                                                  
 * Read the COPYRIGHT and README files in this distribution about
 * registration and installation of Matpack.
 */

#ifndef MATH_WIGNER3JM_HPP
#define MATH_WIGNER3JM_HPP

#include "config.h"

namespace MATPACK {
    
/// return sign of number
template <class T> inline int sign (T x)
{ return (x > 0) ? 1 : (x < 0) ? -1 : 0; }
  
  /*
   * MpMin(), MpMin(): min and max templates for 2 arguments (all that
   * is required for Wigner3jm )
   */
template <class T> inline T MpMin (T x, T y) 
{ return (x<y)?x:y; }
template <class T> inline T MpMax (T x, T y) 
{ return (x>y)?x:y; }

// even and odd

inline int even(int x){return!(x&1);}
inline int odd(int x){return(x&1);}

void Wigner3jm(RealType l1, RealType l2, RealType l3, RealType m1, 
               RealType &m2min, RealType &m2max, RealType *thrcof, int ndim, 
               int &errflag);
}

#endif
