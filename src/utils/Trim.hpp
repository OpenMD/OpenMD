 /*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
 
#ifndef UTILS_TRIM_HPP
#define UTILS_TRIM_HPP

#include <string>
#include <cctype>

/**
 * @file Trim.hpp
 * Defines trim algorithms. Trim algorithms are used to remove trailing and leading spaces from a string. 
 */
namespace oopse {

    /**
     * Remove all leading spaces in-place. The supplied predicate is used to determine which 
     * characters are considered spaces
     * @param str An input sequence
     * @param IsSpace An unary predicate identifying spaces 
     */
    template<typename Predict>     
    void trimLeftIf(std::string& str, Predict isSpace) {
        std::string::iterator i = str.begin();

        for (; i != str.end(); ++i) {
            if (!isSpace(*i)) {
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
    void trimRightIf(std::string& str, Predict isSpace) {
        std::string::iterator i = str.end();

        for (; i != str.begin();) {
            if (!isSpace(*(--i))) {
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
    void trimIf(std::string& str, Predict isSpace) {
        trimLeftIf(str, isSpace);
        trimRightIf(str, isSpace);        
    }

    /**
     * Remove all leading spaces from the input. The supplied predicate is used to determine 
     * which characters are considered spaces
     * @return A trimmed copy of the input
     * @param input An input sequence 
     */
    template<typename Predict>
    std::string trimLeftCopyIf(const std::string& input, Predict isSpace) {
        std::string result(input);
        trimLeftIf(result, isSpace);
        return result;
    }

    /**
     * Remove all trailing spaces from the input. The supplied predicate is used to determine 
     * which characters are considered spaces
     * @return A trimmed copy of the input
     * @param input An input sequence
     */
    template<typename Predict>
    std::string trimRightCopyIf(const std::string& input, Predict isSpace) {
        std::string result(input);
        trimRightIf(result, isSpace);
        return result;
    }

    /**
     * Remove all leading and trailing spaces from the input. The supplied predicate is used to 
     * determine which characters are considered spaces
     * @return A trimmed copy of the input
     * @param input An input sequence
     */
    template<typename Predict>     
    std::string trimCopyIf(const std::string& input, Predict isSpace) {
        std::string result(input);
        trimIf(result, isSpace);
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

}//end namespace oopse
#endif //UTILS_TRIM_HPP    
