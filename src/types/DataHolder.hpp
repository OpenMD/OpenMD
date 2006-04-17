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
 
#ifndef TYPES_DATAHOLDER_HPP
#define TYPES_DATAHOLDER_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <set>
#include "math/Vector3.hpp"
#include "utils/ParameterManager.hpp"
#include "io/ParamConstraint.hpp"
#include "utils/simError.h"
#include "utils/OOPSEException.hpp"
#include "utils/StringUtils.hpp"
namespace oopse {

class DataHolder {
    public:
        DataHolder() {}
        virtual ~DataHolder() {}
        
        template<class T>
        void assign(const std::string& keyword, T val) {
             ParamMap::iterator i =parameters_.find(keyword);
              if (i != parameters_.end()) {
                   bool result = i->second->setData(val);
                   if (!result ) {
                     std::stringstream ss;
              	  ss <<   "Error in parsing " << keyword << ": expected " << i->second->getParamType() <<"\n";
                      throw OOPSEException(ss.str());
                    }
              }else if (deprecatedKeywords_.find(keyword) != deprecatedKeywords_.end()){
                     std::cout << keyword << " has been deprecated in OOPSE 3.  Please update your .md file.\n";
              }else {
                     std::stringstream ss;
                     ss << keyword << " is not a recognized keyword.\n";
                     throw OOPSEException(ss.str());
              }
        }

        virtual void validate() {
          ParamMap::iterator i;
          for (i = parameters_.begin(); i != parameters_.end(); ++i) {
            if (!i->second->isOptional() && i->second->empty()) {
                std::stringstream ss;
                ss <<  i->second->getKeyword()  << " must be set.\n";
                throw OOPSEException(ss.str());
            }
          }
        }
    
    protected:
        typedef std::map<std::string, ParameterBase*> ParamMap;
        
        ParamMap parameters_;        
        std::set<std::string> deprecatedKeywords_;
        
    private:
        DataHolder(const DataHolder&);
        DataHolder& operator=(const DataHolder&);
};

}
#endif
