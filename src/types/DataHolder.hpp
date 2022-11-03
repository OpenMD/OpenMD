/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
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

#ifndef TYPES_DATAHOLDER_HPP
#define TYPES_DATAHOLDER_HPP

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

#include "io/ParamConstraint.hpp"
#include "math/Vector3.hpp"
#include "utils/OpenMDException.hpp"
#include "utils/ParameterManager.hpp"
#include "utils/StringUtils.hpp"
#include "utils/simError.h"
namespace OpenMD {

  class DataHolder {
  public:
    DataHolder() {}
    virtual ~DataHolder() = default;

    template<class T>
    void assign(const std::string& keyword, T val) {
      ParamMap::iterator i = parameters_.find(keyword);
      if (i != parameters_.end()) {
        bool result = i->second->setData(val);
        if (!result) {
          std::stringstream ss;
          ss << "Error in parsing " << keyword << ": expected "
             << i->second->getParamType() << " but got " << &val << "\n";
          throw OpenMDException(ss.str());
        }
      } else if (deprecatedKeywords_.find(keyword) !=
                 deprecatedKeywords_.end()) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "%s keyword has been deprecated in OpenMD. Please update your "
                 ".omd file.\n",
                keyword.c_str());
        painCave.isFatal  = 0;
        painCave.severity = OPENMD_WARNING;
        simError();
      } else {
        std::stringstream ss;
        ss << keyword << " is not a recognized keyword.\n";
        throw OpenMDException(ss.str());
      }
    }

    virtual void validate() {
      ParamMap::iterator i;
      for (i = parameters_.begin(); i != parameters_.end(); ++i) {
        if (!i->second->isOptional() && i->second->empty()) {
          std::stringstream ss;
          ss << i->second->getKeyword() << " must be set.\n";
          throw OpenMDException(ss.str());
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

}  // namespace OpenMD
#endif
