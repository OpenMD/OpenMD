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
  * @file StringTokenizer.hpp
  * @author tlin
  * @date 09/20/2004
  * @time 11:30am
  * @version 1.0
  */
 
#ifndef UTIL_STRINGTOKENIZER_HPP
#define UTIL_STRINGTOKENIZER_HPP

#include <vector>

#include "util/NoSuchElementException.hpp"

namespace oopse{

    /**
     * @class StringTokenizer.hpp "util/StringTokenizer.hpp"
     *
     * @brief The string tokenizer class allows an application to break a string into tokens
     *
     * The set of delimiters (the characters that separate tokens) may be specified either 
     * at creation time or on a per-token basis. 
     * An instance of StringTokenizer behaves in one of two ways, depending on whether it was 
     * created with the returnTokens flag having the value true or false.
     */
    class StringTokenizer{
    
        public:
            
            /**
             * Constructs a string tokenizer for the specified string. The characters in the delim argument
             * are the delimiters for separating tokens. characters are skipped and only serve as 
             * separators between tokens.
             * @param str a string to be parsed.
             * @param delim the delimiters, default value is "\t\n\r".
             */
            StringTokenizer(const std::string& str, const std::string& delim = "\t\n\r");

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
            StringTokenizer(const std::string& str, const std::string& delim, bool returnTokens);

            /**
             * Calculates the number of times that this tokenizer's nextToken method can be called 
             * before it generates an exception.
             *
             * @return the number of tokens remaining in the string using the current delimiter set.
             */
            int countTokens();

            /**
             * Tests if there are more tokens available from this tokenizer's string.
             *
             * @return true if there are more tokens available from this tokenizer's string, false otherwise
             */
            bool hasMoreTokens();

            /**
             * Returns the next token from this string tokenizer.
             *
             * @return the next token from this string tokenizer.
             *
             * @exception NoSuchElementException if there are no more tokens in this tokenizer's string
             */
            std::string nextToken();

            /**
             * Returns the next token in this string tokenizer's string. The new delimiter set remains the
             * default after this call.
             *
             * @param newDelim the new delimiters.
             *
             * @return the next token, after switching to the new delimiter set.
             *
             * @exception NoSuchElementException if there are no more tokens in this tokenizer's string.
             *
             */
            std::string nextToken(const std::string& newDelim); 

            /**
             * Returns the current delimiter set of this string tokenizer
             *
             * @return the current delimiter set
             */
            std::string getDelimiter();

        private:
            
            /**
             * Test if character is in current delimiter set.
             *
             * @param c character to be tested
             *
             * @return true if character is in current delimiter set, flase otherwise.
             */
            bool isDelimiter(char c);
            
            std::string delim_;  /**< current delimiter set of this string tokenizer */

            bool returnTokens_; /**< flag indicating whether to return the delimiters as tokens */
    };


} //namespace oopse
#endif //UTIL_STRINGTOKENIZER_HPP
