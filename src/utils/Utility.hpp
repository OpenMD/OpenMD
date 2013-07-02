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
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#ifndef UTILS_UTILITY_HPP
#define UTILS_UTILITY_HPP
#include <vector>
#include <math.h>
#include "config.h"
#include "utils/next_combination.hpp"

using namespace std;
namespace OpenMD {
  inline RealType roundMe( const RealType &x ){
    return ( x >= 0 ) ? floor( x + 0.5 ) : ceil( x - 0.5 );
  }

  /**
   * @brief iteratively replace the sequence with wild cards
   * @return true if more combination sequence is avaliable, otherwise
   * return true
   * @param cont iterator container, if expect the whole series of
   * combinations, pass an empty iterator container. The user should
   * not modify this iterator container
   * @param sequence the whole sequence used to generate combination
   * @param result a possible combination sequence which is set on return
   * @param wildCard the wild card string. Its value is "X" by default
   * @note since next_combination never returns an empty sequence,
   * replaceWildCard will not generate one special combination, which
   * is n identical wild cards (n is equal to the size of the passing
   * sequence)
   * 
   * @code
   * vector<string> sv;
   * vector<vector<string>::iterator> sic;
   * vector<string> resultString;
   * sv.push_back("H");
   * sv.push_back("C");
   * sv.push_back("N");

   * while (replaceWithWildCard(sic, sv, resultString)) {   
   *     for(vector<string>::iterator i = resultString.begin(); 
   *          i != resultString.end(); ++i) {
   *         cout << *i << "\t";
   *     }
   *     cout << endl;
   * }
   * //output
   * //H X X
   * //X C X
   * //X X N
   * //H C X
   * //H X N
   * //X C N
   * //H C N
   * @endcode
   */
  bool replaceWithWildCard(vector<vector<string>::iterator>& cont,
			   vector<string>& sequence, vector<string>& result, 
                           const string& wildCard = "X");
}
#endif

