/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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
        sprintf(painCave.errMsg,
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
    using ParamMap = std::map<std::string, ParameterBase*>;

    ParamMap parameters_;
    std::set<std::string> deprecatedKeywords_;

  private:
    DataHolder(const DataHolder&);
    DataHolder& operator=(const DataHolder&);
  };
}  // namespace OpenMD

#endif
