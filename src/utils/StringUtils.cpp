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

#include "config.h"
#include <algorithm> 
#include <stdlib.h>
#include <cctype>
#include <cstdlib>
#include <string>
#include "utils/StringUtils.hpp"

#ifdef _MSC_VER
#define strcasecmp _stricmp
#define strdup _strdup
#define strtoull _strtoui64
#endif


namespace OpenMD {
  std::string UpperCase(const std::string& S) {
    std::string uc = S;
    unsigned int n = uc.size();
    for (unsigned int j = 0; j < n; j++) {   
      char sj = uc[j];
      if (sj >= 'a' && sj <= 'z') uc[j] = (char)(sj - ('a' - 'A'));
    }
    return uc;
  }

  std::string LowerCase(const std::string& S) {
    std::string lc = S;
    unsigned int n = lc.size();
    for (unsigned int j = 0; j < n; j++) {
      char sj = lc[j];
      if (sj >= 'A' && sj <= 'Z') lc[j] = (char)(sj + ('a' - 'A'));
    }
    return lc;
  }

  char* trimSpaces(char *str) {
    size_t len;
    char *right, *left;

    if (strlen(str) == 0) return(str);
  
    /* Trim whitespace from left side */
    for (left = str; isspace(*left); left++);

    /* Trim whitespace from right side */
    if ((len = strlen(left)))
      {
        right = left + (len - 1);
      
        while (isspace(*right))
          {
            *right = '\0';
            right--;
          }
      }

    /* Only do the str copy if there were spaces to the left */
    if (left != str)
      strcpy(str, left);
  
    return (str);
  }

  int findBegin(std::istream &theStream, const char* startText ){
    const int MAXLEN = 1024;
    char readLine[MAXLEN];   
    int foundText = 0;
    int lineNum;
    char* the_token;

    // rewind the stream
    theStream.seekg (0, std::ios::beg);
    lineNum = 0;
  
    if (!theStream.eof()) {
      theStream.getline(readLine, MAXLEN);
      lineNum++;
    } else {
      printf( "Error fast forwarding stream: stream is empty.\n");
      return -1;
    }
  
    while ( !foundText ) {
    
      if (theStream.eof()) {
        printf("Error fast forwarding stream at line %d: "
               "stream ended unexpectedly.\n", lineNum);
        return -1;
      }
    
      the_token = strtok( readLine, " ,;\t" );
      if ( the_token != NULL)
        if (!strcasecmp("begin", the_token)) {
          the_token = strtok( NULL, " ,;\t" );
          if ( the_token != NULL){
            foundText = !strcasecmp( startText, the_token );
          }
        }
    
      if (!foundText) {
        if (!theStream.eof()) {
          theStream.getline(readLine, MAXLEN);
          lineNum++;
        } else {
          printf( "Error fast forwarding stream at line %d: "
                  "stream ended unexpectedly.\n", lineNum);
          return -1;
        }
      }
    }
    return lineNum;
  }

  int countTokens(char *line, char *delimiters) {
    /* PURPOSE: RETURN A COUNT OF THE NUMBER OF TOKENS ON THE LINE. */
  
    char *working_line;   /* WORKING COPY OF LINE. */
    int ntokens;          /* NUMBER OF TOKENS FOUND IN LINE. */
    char *strtok_ptr;     /* POINTER FOR STRTOK. */

    strtok_ptr= working_line= strdup(line);

    ntokens=0;
    while (strtok(strtok_ptr,delimiters)!=NULL)
      {
        ntokens++;
        strtok_ptr=NULL;
      }

    free(working_line);
    return(ntokens);
  }

  int isEndLine(char *line) {
    char *working_line;
    char *foo;
  
    working_line = strdup(line);
  
    foo = strtok(working_line, " ,;\t");

    if (foo != NULL) {

      if (!strcasecmp(foo, "end")) return 1;

    }
 
    return 0;
  }
  
  std::string OpenMD_itoa(int value, unsigned int base) {    
    const char digitMap[] = "0123456789abcdef";	
    std::string buf;

    if (base == 0 || base > 16) {      
      return buf;
    }

    if (value == 0) {
      buf = "0";
      return buf;
    }
    
    // Take care negative int:
	
    std::string sign;    
    int _value = value;	
    if (value < 0) {
      _value = -value;
      sign = "-";
    }

    // Translating number to string with base:
    for (int i = 30; _value && i ; --i) {
      buf = digitMap[ _value % base ] + buf;
      _value /= base;
    }
    return sign.append(buf);
  }


  std::string getPrefix(const std::string& str) {
    return str.substr(0, str.rfind('.'));
  }

  std::string getSuffix(const std::string& str) {
    return str.substr(0, str.find('.'));
  }
  
  bool isInteger(const std::string& str) {
    
    bool result = false;
    
    std::string::const_iterator i = str.begin();    
    if (i != str.end() && (*i == '+' || *i == '-' || std::isdigit(*i) )) {
      ++i;        
      while (i != str.end() && std::isdigit(*i))
        ++i;
      if (i == str.end())
        result = true;
    }
    
    return result;
  }

  bool CaseInsensitiveEquals(const char ch1, const char ch2) {
    return std::toupper((unsigned char)ch1) == std::toupper((unsigned char)ch2);
  }
  
  size_t CaseInsensitiveFind(const std::string& str1, const std::string& str2) {
    std::string::const_iterator pos = std::search(str1.begin(), str1.end(), str2.begin(), str2.end(), CaseInsensitiveEquals);
    if (pos == str1.end())
      return std::string::npos;
    else
      return pos - str1.begin();
  }
  
  /**
   *    memparse - parse a string with mem suffixes into a number
   *    @param ptr: Where parse begins
   *    @param retptr: (output) Pointer to next char after parse completes
   *
   *    Parses a string into a number.  The number stored at @param ptr is
   *    potentially suffixed with %K (for kilobytes, or 1024 bytes),
   *    %M (for megabytes, or 1048576 bytes), or %G (for gigabytes, or
   *    1073741824).  If the number is suffixed with K, M, or G, then
   *    the return value is the number multiplied by one kilobyte, one
   *    megabyte, or one gigabyte, respectively.
   */  
  unsigned long long memparse (char *ptr,  char **retptr) {
    unsigned long long ret = strtoull (ptr, retptr, 0);
    
    switch (**retptr) {
    case 'G':
    case 'g':
      ret <<= 10;
    case 'M':
    case 'm':
      ret <<= 10;
    case 'K':
    case 'k':
      ret <<= 10;
      (*retptr)++;
    default:
      break;
    }
    return ret;
  }
  
}
