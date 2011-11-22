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
 * [4]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).                        
 */
 
#ifndef UTILS_TRIM_HPP
#define UTILS_TRIM_HPP

#include <string>
#include <cctype>
#include "utils/Predicate.hpp"

/**
 * @file Trim.hpp
 * Defines trim algorithms. Trim algorithms are used to remove trailing and leading spaces from a string. 
 */
namespace OpenMD {

  /**
   * Remove all leading spaces in-place. The supplied predicate is used to determine which 
   * characters are considered spaces
   * @param str An input sequence
   * @param IsSpace An unary predicate identifying spaces 
   *
   * @code
   * std::string str = "  acb  trimLeftIf test case"
   * trimLeftIf(str, pred() || isFromRange('a', 'c')); 
   * std::cout << str << std::endl; //print "trimLeft test case"
   * 
   * @endcode
   */
  template<typename Predict>     
  void trimLeftIf(std::string& str, Predict pred) {
    std::string::iterator i = str.begin();

    for (; i != str.end(); ++i) {
      if (!pred(*i)) {
	break;
      }
    }
        
    str.erase(str.begin(), i);
  }

  /**
   * Remove all trailing spaces in-place. The supplied predicate is used to determine which 
   * characters are considered spaces
   * @param str An input sequence 
   */
  template<typename Predict>     
  void trimRightIf(std::string& str, Predict pred) {
    std::string::iterator i = str.end();

    for (; i != str.begin();) {
      if (!pred(*(--i))) {
	++i;
	break;
      }
    }
        
    str.erase(i, str.end());
  }

  /**
   *Remove all leading and trailing spaces in-place. The supplied predicate is used to determine 
   * which characters are considered spaces
   * @param str An input sequence
   */
  template<typename Predict>     
  void trimIf(std::string& str, Predict pred) {
    trimLeftIf(str, pred);
    trimRightIf(str, pred);        
  }

  /**
   * Remove all leading spaces from the input. The supplied predicate is used to determine 
   * which characters are considered spaces
   * @return A trimmed copy of the input
   * @param input An input sequence 
   */
  template<typename Predict>
  std::string trimLeftCopyIf(const std::string& input, Predict pred) {
    std::string result(input);
    trimLeftIf(result, pred);
    return result;
  }

  /**
   * Remove all trailing spaces from the input. The supplied predicate is used to determine 
   * which characters are considered spaces
   * @return A trimmed copy of the input
   * @param input An input sequence
   */
  template<typename Predict>
  std::string trimRightCopyIf(const std::string& input, Predict pred) {
    std::string result(input);
    trimRightIf(result, pred);
    return result;
  }

  /**
   * Remove all leading and trailing spaces from the input. The supplied predicate is used to 
   * determine which characters are considered spaces
   * @return A trimmed copy of the input
   * @param input An input sequence
   */
  template<typename Predict>     
  std::string trimCopyIf(const std::string& input, Predict pred) {
    std::string result(input);
    trimIf(result, pred);
    return result;
  }

    
  /**
   * Remove all leading spaces in-place.
   * @param str An input sequence
   */
  void trimLeft(std::string& str);

  /**
   * Remove all trailing spaces in-place.
   * @param str An input sequence 
   */
  void trimRight(std::string& str);

  /**
   *Remove all leading and trailing spaces in-place
   * @param str An input sequence
   */
  void trim(std::string& str);

  /**
   * Remove all leading spaces from the input.
   * @return A trimmed copy of the input
   * @param input An input sequence 
   */
  std::string trimLeftCopy(const std::string& input);

  /**
   * Remove all trailing spaces from the input.
   * @return A trimmed copy of the input
   * @param input An input sequence
   */
  std::string trimRightCopy(const std::string& input);

  /**
   *Remove all leading and trailing spaces from the input.
   * @return A trimmed copy of the input
   * @param input An input sequence
   */
  std::string trimCopy(const std::string& input);

}//end namespace OpenMD
#endif //UTILS_TRIM_HPP    
