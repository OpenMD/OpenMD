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
 * @file StringTokenizer.hpp
 * @author tlin
 * @date 09/20/2004
 * @time 11:30am
 * @version 1.0
 */

#ifndef UTIL_STRINGTOKENIZER_HPP
#define UTIL_STRINGTOKENIZER_HPP

#include <string>
#include <stdlib.h>
#include <vector>
#include "config.h"
namespace OpenMD {

  /**
   * @class StringTokenizer.hpp "util/StringTokenizer.hpp"
   * @brief The string tokenizer class allows an application to break a string into tokens
   * The set of delimiters (the characters that separate tokens) may be specified either 
   * at creation time or on a per-token basis. 
   * An instance of StringTokenizer behaves in one of two ways, depending on whether it was 
   * created with the returnTokens flag having the value true or false.
   */
  class StringTokenizer {
  public:

    /**
     * Constructs a string tokenizer for the specified string. The characters in the delim argument
     * are the delimiters for separating tokens. characters are skipped and only serve as 
     * separators between tokens.
     * @param str a string to be parsed.
     * @param delim the delimiters, default value is " ;\t\n\r".
     * @note this is still a little bit java like implementation. Pure c++ one should use TokenIterator.
     * Boost's tokenizer class is one of them 
     */
    StringTokenizer(const std::string & str,
		    const std::string & delim = " ;\t\n\r");

    /**
     * Constructs a string tokenizer for an iterator range [first, last). The characters in the delim argument
     * are the delimiters for separating tokens. characters are skipped and only serve as 
     * separators between tokens.
     * @param first begin iterator
     * @param last end iterator
     * @param delim the delimiters, default value is " ;\t\n\r".
     * @note this is still a little bit java like implementation. Pure c++ one should use TokenIterator.
     * Boost's tokenizer class is one of them 
     */
    StringTokenizer(std::string::const_iterator& first, std::string::const_iterator& last,
		    const std::string & delim = " ;\t\n\r");

    /**
     * Constructs a string tokenizer for the specified string. The characters in the delim argument
     * are the delimiters for separating tokens. 
     * If the returnTokens flag is true, then the delimiter characters are also returned as tokens. 
     * Each delimiter is returned as a string of length one. If the flag is false, the delimiter 
     * characters are skipped and only serve as separators between tokens.
     * @param str a string to be parsed. 
     * @param delim the delimiters. 
     * @param returnTokens flag indicating whether to return the delimiters as tokens.
     */
    StringTokenizer(const std::string&str, const std::string&delim,
		    bool returnTokens);

    /**
     * Calculates the number of times that this tokenizer's nextToken method can be called 
     * before it generates an exception.
     * @return the number of tokens remaining in the string using the current delimiter set.
     */
    int countTokens();

    /**
     * Tests if there are more tokens available from this tokenizer's string.
     * @return true if there are more tokens available from this tokenizer's string, false otherwise
     */
    bool hasMoreTokens();

    /**
     * Returns the next token from this string tokenizer.
     * @return the next token from this string tokenizer.
     * @exception NoSuchElementException if there are no more tokens in this tokenizer's string
     */
    std::string nextToken();

    //actually, nextToken Can be template function
    //template <typename ReturnType>
    //ReturnType nextToken();
        
    /**
     * Returns the next token from this string tokenizer as a bool.
     * @return the next token from this string tokenizer  as a bool.
     */
    bool nextTokenAsBool();

    /**
     * Returns the next token from this string tokenizer as an integer.
     * @return the next token from this string tokenizer  as an integer.
     */
    int nextTokenAsInt();

    /**
     * Returns the next token from this string tokenizer as a float.
     * @return the next token from this string tokenizer as a float.
     */
    float nextTokenAsFloat();

    /**
     * Returns the next token from this string tokenizer as a RealType.
     * @return the next token from this string tokenizer as a RealType.
     */
    RealType nextTokenAsDouble();

    /**
     * Returns the next token without advancing the position of the StringTokenizer.
     * @return the next token
     */
    std::string  peekNextToken();

    /**
     * Returns the current delimiter set of this string tokenizer
     * @return the current delimiter set
     */
    const std::string& getDelimiters() {
      return delim_;
    }

    /** 
     * Returns the original string before tokenizing.
     * @return the original string before tokenizing 
     */
    const std::string& getOriginal() {
      return tokenString_;
    }

    /** 
     * Returns all of the tokens
     * @return all of the tokens
     */
    std::vector<std::string> getAllTokens();
  private:

    /**
     * Test if character is in current delimiter set.
     * @param c character to be tested
     * @return true if character is in current delimiter set, flase otherwise.
     */
    bool isDelimiter(const char c);

    /** convert a fortran number to a c/c++ number */
    void convertFortranNumber(std::string& fortranNumber);
         

    std::string tokenString_;

    std::string delim_;         /**< current delimiter set of this string tokenizer */

    bool returnTokens_; /**< flag indicating whether to return the delimiters as tokens */

    std::string::const_iterator currentPos_;
    std::string::const_iterator end_;
  };

}                               //namespace OpenMD

#endif                          //UTIL_STRINGTOKENIZER_HPP
