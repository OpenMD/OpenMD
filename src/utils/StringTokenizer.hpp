/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

/**
 * @file StringTokenizer.hpp
 * @author tlin
 * @date 09/20/2004
 * @version 1.0
 */

#ifndef UTIL_STRINGTOKENIZER_HPP
#define UTIL_STRINGTOKENIZER_HPP

#include <config.h>

#include <cstdlib>
#include <string>
#include <vector>

namespace OpenMD {

  /**
   * @class StringTokenizer
   * @brief The string tokenizer class allows an application to break a string
   * into tokens The set of delimiters (the characters that separate tokens) may
   * be specified either at creation time or on a per-token basis. An instance
   * of StringTokenizer behaves in one of two ways, depending on whether it was
   * created with the returnTokens flag having the value true or false.
   */
  class StringTokenizer {
  public:
    /**
     * Constructs a string tokenizer for the specified string. The characters in
     * the delim argument are the delimiters for separating tokens. characters
     * are skipped and only serve as separators between tokens.
     * @param str a string to be parsed.
     * @param delim the delimiters, default value is " ;\t\n\r".
     * @note this is still a little bit java like implementation. Pure c++ one
     * should use TokenIterator. Boost's tokenizer class is one of them
     */
    StringTokenizer(const std::string& str,
                    const std::string& delim = " ;\t\n\r");

    /**
     * Constructs a string tokenizer for an iterator range [first, last). The
     * characters in the delim argument are the delimiters for separating
     * tokens. characters are skipped and only serve as separators between
     * tokens.
     * @param first begin iterator
     * @param last end iterator
     * @param delim the delimiters, default value is " ;\t\n\r".
     * @note this is still a little bit java like implementation. Pure c++ one
     * should use TokenIterator. Boost's tokenizer class is one of them
     */
    StringTokenizer(std::string::const_iterator& first,
                    std::string::const_iterator& last,
                    const std::string& delim = " ;\t\n\r");

    /**
     * Constructs a string tokenizer for the specified string. The characters in
     * the delim argument are the delimiters for separating tokens. If the
     * returnTokens flag is true, then the delimiter characters are also
     * returned as tokens. Each delimiter is returned as a string of length one.
     * If the flag is false, the delimiter characters are skipped and only serve
     * as separators between tokens.
     * @param str a string to be parsed.
     * @param delim the delimiters.
     * @param returnTokens flag indicating whether to return the delimiters as
     * tokens.
     */
    StringTokenizer(const std::string& str, const std::string& delim,
                    bool returnTokens);

    /**
     * Calculates the number of times that this tokenizer's nextToken method can
     * be called before it generates an exception.
     * @return the number of tokens remaining in the string using the current
     * delimiter set.
     */
    int countTokens();

    /**
     * Tests if there are more tokens available from this tokenizer's string.
     * @return true if there are more tokens available from this tokenizer's
     * string, false otherwise
     */
    bool hasMoreTokens();

    /**
     * Returns the next token from this string tokenizer.
     * @return the next token from this string tokenizer.
     * @exception NoSuchElementException if there are no more tokens in this
     * tokenizer's string
     */
    std::string nextToken();

    /**
     * Skips the next token from this string tokenizer.
     * @exception NoSuchElementException if there are no more tokens in this
     * tokenizer's string
     */
    void skipToken();

    // actually, nextToken Can be template function
    // template <typename ReturnType>
    // ReturnType nextToken();

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
     * Returns the next token without advancing the position of the
     * StringTokenizer.
     * @return the next token
     */
    std::string peekNextToken();

    /**
     * Returns the current delimiter set of this string tokenizer
     * @return the current delimiter set
     */
    const std::string& getDelimiters() { return delim_; }

    /**
     * Returns the original string before tokenizing.
     * @return the original string before tokenizing
     */
    const std::string& getOriginal() { return tokenString_; }

    /**
     * Returns all of the tokens
     * @return all of the tokens
     */
    std::vector<std::string> getAllTokens();

    /**
     * Returns the remaining unparsed string
     * @return the remaining unparsed string
     */
    std::string getRemainingString() const;

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

    std::string delim_; /**< current delimiter set of this string tokenizer */

    bool returnTokens_; /**< flag indicating whether to return the delimiters as
                           tokens */

    std::string::const_iterator currentPos_;
    std::string::const_iterator end_;
  };
}  // namespace OpenMD

#endif  // UTIL_STRINGTOKENIZER_HPP
