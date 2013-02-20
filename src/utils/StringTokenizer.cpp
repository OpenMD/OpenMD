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
 
#include <iostream>
#include <iterator>
#include <sstream>
#include "utils/StringTokenizer.hpp"

namespace OpenMD {


  StringTokenizer::StringTokenizer(const std::string & str, const std::string & delim) 
    : tokenString_(str), delim_(delim), returnTokens_(false),
      currentPos_(tokenString_.begin()), end_(tokenString_.end()){

    }

  StringTokenizer::StringTokenizer(std::string::const_iterator& first, std::string::const_iterator& last,
				   const std::string & delim)  
    : tokenString_(first, last) , delim_(delim), returnTokens_(false),
      currentPos_(tokenString_.begin()), end_(tokenString_.end()) {

    }

  StringTokenizer::StringTokenizer(const std::string&str, const std::string&delim,
				   bool returnTokens)
    : tokenString_(str), delim_(delim), returnTokens_(returnTokens),
      currentPos_(tokenString_.begin()), end_(tokenString_.end()) {

    }

  bool StringTokenizer::isDelimiter(const char c) {
    return delim_.find(c) == std::string::npos ? false : true;
  }

  int StringTokenizer::countTokens() {
    
    std::string::const_iterator tmpIter = currentPos_;    
    int numToken = 0;

    while (true) {

      //skip delimiter first
      while( tmpIter != end_ && isDelimiter(*tmpIter)) {
	++tmpIter;

	if (returnTokens_) {
	  //if delimiter is consider as token
	  ++numToken;
	}
      }
        
      if (tmpIter == end_) {
	break;
      }
        
      //encount a token here
      while ( tmpIter != end_ && !isDelimiter(*tmpIter) ) {
	++tmpIter;
      }

      ++numToken;

    }

    return numToken;
  }

  bool StringTokenizer::hasMoreTokens() {
    
    if (currentPos_ == end_) {
      return false;
    } else if (returnTokens_) {
      return true;
    } else {
      std::string::const_iterator i = currentPos_;

      //walk through the remaining string to check whether it contains non-delimeter or not
      while(i != end_ && isDelimiter(*i)) {
	++i;
      }

      return i != end_ ? true : false;
    }
  }

  std::string StringTokenizer::nextToken() {
    std::string result;
    
    if(currentPos_ != end_) {
      std::insert_iterator<std::string> insertIter(result, result.begin());

      while( currentPos_ != end_ && isDelimiter(*currentPos_)) {

	if (returnTokens_) {
	  *insertIter++ = *currentPos_++;
	  return result;
	}
            
	++currentPos_;
      }

      while (currentPos_ != end_ && !isDelimiter(*currentPos_)) {
	*insertIter++ = *currentPos_++;
      }
        
    }
    
    return result;
  }

  bool StringTokenizer::nextTokenAsBool() {
    std::string token = nextToken();
    std::istringstream iss(token);
    bool result;
    
    if (iss >> result) {
      return result;
    } else {
      std::cerr << "unable to convert " << token << " to a bool" << std::endl;
      return false;
    }
  }
 
  //Since libstdc++(GCC 3.2) has an i/ostream::operator>>/<<(streambuf*) bug (Bug 9318)
  //Instead of using iostream facility, we use C library
  int StringTokenizer::nextTokenAsInt() {
    std::string token = nextToken();
   
    return atoi(token.c_str());
  }

  float StringTokenizer::nextTokenAsFloat() {
    std::string token = nextToken();
    convertFortranNumber(token);
    return (float) (atof(token.c_str()));
  }

  RealType StringTokenizer::nextTokenAsDouble() {
    std::string token = nextToken();
    convertFortranNumber(token);
    return atof(token.c_str());
  }

  std::string  StringTokenizer::peekNextToken() {
    std::string result;
    std::string::const_iterator tmpIter = currentPos_;
    
    if(tmpIter != end_) {
      std::insert_iterator<std::string> insertIter(result, result.begin());

      while(tmpIter != end_ && isDelimiter(*tmpIter)) {

	if (returnTokens_) {
	  *insertIter++ = *tmpIter++;
	  return result;
	}
            
	++tmpIter;
      }

      while (tmpIter != end_ && !isDelimiter(*tmpIter)) {
	*insertIter++ = *tmpIter++;
      }
    }
    
    return result;    
  }

 std::vector<std::string>  StringTokenizer::getAllTokens() {
    std::vector<std::string> tokens;
    while (hasMoreTokens()) {
        tokens.push_back(nextToken());
    }
    return tokens;
 }
  void StringTokenizer::convertFortranNumber(std::string& fortranNumber) {
    std::string::iterator i;
    for(i = fortranNumber.begin(); i != fortranNumber.end(); ++i) {
      if (*i == 'd' || *i == 'D') {
	*i = 'E';
      }
    }
  }

}//end namespace OpenMD

