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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).                        
 */
 
#ifndef UTILS_TUPLE_HPP
#define UTILS_TUPLE_HPP

namespace OpenMD {

  template <class T1, class T2, class T3>
  struct tuple3 {
    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;

    T1 first;
    T2 second;
    T3 third;
    tuple3() {}
    tuple3(const T1& a, const T2& b, const T3& c) : first(a), second(b), third(c) {}
  };

  template <class T1, class T2, class T3>
  tuple3<T1,T2,T3> make_tuple3( T1 t1, T2 t2, T3 t3 ) {
    return tuple3<T1,T2,T3>( t1, t2, t3 ); 
  }


  template <class T1, class T2, class T3, class T4>
  struct tuple4 {
    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;
    typedef T4 fourth_type;

    T1 first;
    T2 second;
    T3 third;
    T4 fourth;
    tuple4() {}
    tuple4(const T1& a, const T2& b, const T3& c, const T4& d)
      : first(a), second(b), third(c), fourth(d) {}
  };

  template <class T1, class T2, class T3, class T4>
  tuple4<T1,T2,T3,T4> make_tuple4( T1 t1, T2 t2, T3 t3, T4 t4 ) {
    return tuple4<T1,T2,T3,T4>( t1, t2, t3, t4 );
  }

template<class T1, class T2, class T3>
inline bool operator < (const tuple3<T1, T2, T3>& t1, const tuple3<T1, T2, T3>& t2) {

        return t1.first < t2.first
             || (!(t2.first < t1.first) && t1.second < t2.second)
             || (!(t2.first < t1.first) && !(t2.second < t2.second) && t1.third < t2.third);
}


inline bool operator < (const tuple3<int, int, std::vector<std::string> >& t1, const tuple3<int, int, std::vector<std::string> >& t2) {
  
  if (t1.first < t2.first)
    return true;
  else {
    if (t1.first > t2.first) 
      return false;
    
    if (t1.second < t2.second) 
      return true;
    else 
      if (t1.second > t2.second) 
	return false;
    
    return true;
  }  
}

template<class T1, class T2, class T3, class T4>
inline bool operator < (const tuple4<T1, T2, T3, T4>& t1, const tuple4<T1, T2, T3, T4>& t2) {

        return t1.first < t2.first
             || (!(t2.first < t1.first) && t1.second < t2.second)
             || (!(t2.first < t1.first) && !(t2.second < t2.second) && t1.third < t2.third)
             ||(!(t2.first < t1.first) && !(t2.second < t2.second) && !(t2.third < t1.third) && t1.fourth < t2.fourth);
}
typedef tuple3<int, int, int> IntTuple3;
typedef tuple4<int, int, int, int> IntTuple4;

}


#endif //UTILS_TUPLE_HPP

