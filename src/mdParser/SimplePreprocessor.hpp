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

#ifndef MDPARSER_SIMPLEPREPROCESSOR_HPP
#define MDPARSER_SIMPLEPREPROCESSOR_HPP
#include <iostream>
#include <set>
#include <fstream>
#include <sstream>
#include "utils/StringTokenizer.hpp"
#include "utils/Trim.hpp"
#include "utils/OpenMDException.hpp"
#include "utils/simError.h"


/**
 * @class SimplePreprocessor
 * @brief A simple preprocessor.
 * @note only supports \#include \#ifdef, \#ifndef, \#endif, \#define and \#undef, c-like multiple line
 *  comment is not supported, macro substitute is not supported.
 */
namespace OpenMD { 
class SimplePreprocessor {
    public:
        bool preprocess(std::istream& myStream, const std::string& filename, int startingLine, ostream& os) {
            std::set<std::string> defineSet;
            std::stack<bool> ifStates;

            ifStates.push(true);
            return doPreprocess(myStream, filename, startingLine, os, defineSet, ifStates);
        }
        
    private:
        bool doPreprocess(std::istream& myStream, const std::string& filename, int startingLine, ostream& os, std::set<std::string>& defineSet, std::stack<bool>& ifStates) {
            //std::ifstream input(filename.c_str());
            //if (!input.is_open()) {
            //    std::stringstream ss;
            //    ss << "Can not open " << filename << " for preprocessing\n";
            //    
            //    sprintf(painCave.errMsg,
            //            "Can not open (%s) for processing. \n"
            //            "\tPlease check md file name syntax.\n", filename.c_str());
            //    
            //    painCave.isFatal = 1;
            //    simError();
            //    
            //    throw OpenMDException(ss.str());                
            //}
            int lineNo = startingLine;
            os << "#line " << lineNo << " \"" << filename << "\"\n";
            const int bufferSize = 1024;
            char buffer[bufferSize];
            while(myStream.getline(buffer, bufferSize)) {
              ++lineNo;
              std::string line = trimLeftCopy(buffer);
              if (!line.empty() && line[0] == '#') {
                    StringTokenizer tokenizer(line.substr(1, line.length()));
                    std::vector<std::string> tokens = tokenizer.getAllTokens();
                    if (tokens.size() < 1 ) {
                        return false;
                    }
                    std::string keyword = tokens[0];
                    if (tokens[0] == "endif") {
                        ifStates.pop();
                        if (ifStates.empty()) {
                            std::cout << "Error in preprocessing: endif \n";
                            return false;
                        }
                        os << std::endl;                        
                    } else if (tokens.size() == 2) {
                        if (tokens[0] == "include") {
                            SimplePreprocessor subPreprocessor;
                            std::string includeFilename = tokens[1];
                            includeFilename = includeFilename.substr(1, includeFilename.length() -2);
                            std::ifstream includeStream(includeFilename.c_str());
                            if (!includeStream.is_open()) {
                                std::stringstream ss;
                                ss << "Can not open " << includeFilename << " for preprocessing\n";
                                throw OpenMDException(ss.str()); 
                            }
                            
                            bool ret = subPreprocessor.doPreprocess(includeStream, includeFilename, 1, os, defineSet, ifStates);
                            if (!ret) {
                                std::cout << "Error in preprocessing\n";
                                return false;
                            }
                            os << "#line " << lineNo << " \"" << filename << "\"\n";
                        } else if (tokens[0] == "define") {
                           defineSet.insert(tokens[1]);
                           os << std::endl;
                        } else if (tokens[0] == "undef") {
                           defineSet.erase(tokens[1]);
                           os << std::endl;
                        } else if (tokens[0] == "ifdef") {
                           if (defineSet.find(tokens[1]) != defineSet.end() ) {
                                ifStates.push(true);
                           } else {
                              ifStates.push(false);
                           }
                           os << std::endl;
                        } else if (tokens[0] == "ifndef") {
                           if (defineSet.find(tokens[1]) == defineSet.end() ) {
                                ifStates.push(true);
                           } else {
                              ifStates.push(false);
                           }
                           os << std::endl;
                        } else {
                            std::cout << tokens[0] << " is not supported (yet)." << std::endl;
                            return false;
                        }
                    }else {
                        return false;
                    }
                    
              }else if (ifStates.top()){
                os << buffer << std::endl;
              }
              
            }

            return true;
        }
    private:
        
};

}
#endif
