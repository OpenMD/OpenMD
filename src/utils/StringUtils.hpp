/**
 * @file StringUtils.hpp
 * @author Dan Gezelter
 * @date 10/18/2004
 * @version 1.0
 */

#ifndef UTILS_STRINGUTILS_HPP
#define UTILS_STRINGUTILS_HPP
#include <string>
#include <iostream>
#include <fstream>

namespace oopse {

  using namespace std;
 

  /**
   * Converts a string to UPPER CASE
   * @param S
   */
  string UpperCase(const string& S);

  /**
   * Converts a string to lower case
   * @param S
   */
  string LowerCase(const string& S);

  /**
   * Finds the location of the string "begin <startText>" in an input stream.
   * @param theStream
   * @param startText
   *
   * @return the line number of the block within the theStream
   */
  int findBegin(istream theStream, char* startText );

  /**
   * Counts the number of tokens on line which are delimited by the characters 
   * listed in delimiters
   * @param line
   * @param delimiters
   */
  int countTokens(char *line, char *delimiters);

  /**
   * Removes left and right spaces from a string
   *
   * @param str  String to trim
   *
   * @return  char* to the trimed string
   */
  char* TrimSpaces(char *str);

  /**
   * discovers whether or not the line contains the "end" token
   *
   * @param line  The line to test
   * 
   * @return int  (==1 if the line has "end", ==0 if not).
   */
  int isEndLine(char *line);
}

#endif
