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
 
#include "io/ZConsReader.hpp"
#include "utils/simError.h"
#include "utils/StringUtils.hpp"
namespace oopse {

  ZConsReader::ZConsReader(SimInfo* info) : info_(info){

    std::string zconsFileName_ = getPrefix(info_->getFinalConfigFileName()) + ".fz";
    istream_.open(zconsFileName_.c_str());

    if (!istream_){
      std::cerr << "open " << zconsFileName_ << "error" << std::endl;
      exit(1);
    }

    Globals* simParam = info_->getSimParams();
    int nZconstraints = simParam->getNZconsStamps();
    std::vector<ZConsStamp*> stamp = simParam->getZconsStamps();
    for (int i = 0; i < nZconstraints; i++){
      allZmols_.push_back(stamp[i]->getMolIndex());
    } 
  }

  ZConsReader::~ZConsReader(){
    istream_.close();
  }

  void ZConsReader::readNextFrame(){
    char line[MAXBUFFERSIZE];  
    int nFixedZmol;
    int sscanfCount;

    fixedZmolData_.clear();

    while(istream_.getline(line, MAXBUFFERSIZE) && line[0] == '/' && line[1] == '/');
    sscanfCount = sscanf(line, "%lf", &curTime_);
    if (sscanfCount != 1){
      std::cerr << "ZConsReader Error : reading file error" << std::endl;
      exit(1);
    }

    istream_.getline(line, MAXBUFFERSIZE);
    sscanfCount = sscanf(line, "%d", &nFixedZmol);
    if (sscanfCount != 1){
      std::cerr << "ZConsReader Error : reading file error" << std::endl;
      exit(1);
    }

    ZconsData data;
    for(int i = 0; i < nFixedZmol; i++){
      istream_.getline(line, MAXBUFFERSIZE);
      sscanfCount = sscanf(line, "%d\t%lf\t%lf\t%lf", &data.zmolIndex, &data.zforce, &data.zpos,&data.zconsPos);
      if (sscanfCount != 4){
	std::cerr << "ZConsReader Error : reading file error" << std::endl;
	exit(1);
      }

      fixedZmolData_.push_back(data);
    }

  }

  bool ZConsReader::hasNextFrame(){
    return istream_.peek() != EOF ? true : false;
  }

}
