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
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
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
            //    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
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
            const int bufferSize = 8192;
            char buffer[bufferSize];
            while(myStream.getline(buffer, bufferSize)) {
              ++lineNo;
              std::string line = Utils::trimLeftCopy(buffer);
              if (!line.empty() && line[0] == '#') {
                    StringTokenizer tokenizer(line.substr(1, line.length()));
                    std::vector<std::string> tokens = tokenizer.getAllTokens();
                    if (tokens.empty()) {
                        return false;
                    }
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
