/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * @file RestReader.hpp 
 */ 
 
#ifndef IO_RESTREADER_HPP 
#define IO_RESTREADER_HPP 
 
#include <cstdio>
#include <string>
#include "brains/SimInfo.hpp"
#include "primitives/StuntDouble.hpp"

namespace OpenMD { 
 
  /** 
   * @class RestReader RestReader.hpp "io/RestReader.hpp" 
   */ 
  class RestReader { 
  public: 
    RestReader(SimInfo* info, const std::string &filename, 
                           std::vector<int> stuntDoubleIndex) 
      : info_(info), stuntDoubleIndex_(stuntDoubleIndex), filename_(filename) {
#ifdef IS_MPI 
      
      if (worldRank == 0) { 
#endif 
        
        inFile_ = new std::ifstream(filename_.c_str(),
                                    ifstream::in | ifstream::binary); 
        
        if (inFile_->fail()) { 
          sprintf(painCave.errMsg, 
                  "RestReader: Cannot open file: %s\n", 
                  filename_.c_str()); 
          painCave.isFatal = 1; 
          simError(); 
        } 
        
#ifdef IS_MPI 
        
      } 
      
      strcpy(checkPointMsg, "Reference file opened for reading successfully."); 
      errorCheckPoint(); 
      
#endif 
      
      return;
    }

    void readReferenceStructure(void);

    ~RestReader(){
#ifdef IS_MPI
      if (worldRank == 0){
#endif 
        
        delete inFile_;
        
#ifdef IS_MPI       
      } 
      
      strcpy(checkPointMsg, "Reference file closed successfully."); 
      errorCheckPoint();
#endif
    }
    
  protected: 
    void parseDumpLine(const std::string& line);
    void readStuntDoubles(std::istream& inpuStream);
    void readFrameProperties(std::istream& inputStream);
    void scanFile(void);
    void readSet(void);
  private: 

    SimInfo* info_ {nullptr};

    std::vector<Vector3d> all_pos_;
    std::vector<int> stuntDoubleIndex_;

    std::istream* inFile_;
    std::string filename_;

    long long framePos_;

    const static int bufferSize = 4096;
    char buffer[bufferSize];
  };
 
}      //end namespace OpenMD 
 
#endif //IO_RESTREADER_HPP 
