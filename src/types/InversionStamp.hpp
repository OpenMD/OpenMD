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

#ifndef TYPES_INVERSIONSTAMP_HPP
#define TYPES_INVERSIONSTAMP_HPP

#include "types/DataHolder.hpp"

namespace OpenMD {
  class InversionStamp : public DataHolder {
    DeclareParameter(GhostVectorSource, int);

  public:
    InversionStamp();
    virtual ~InversionStamp();

    int getCenter() { return center_; }
    int getSatelliteAt(int index) { return satellites_.at(index); }
    int getNSatellites() { return satellites_.size(); }
    std::vector<int> getSatellites() { return satellites_; }
    void setCenter(int center) { center_ = center; }
    void addSatellite(int sat) {
      if (satellites_.size() > 3) {
        std::ostringstream oss;
        oss << "Too many satellites in inversion to add another!" << std::endl;
        throw OpenMDException(oss.str());
      } else {
        satellites_.push_back(sat);
      }
    }
    void setSatellites(const std::vector<int>& sats) {
      if (sats.size() == 3) {
        satellites_.push_back(sats.at(0));
        satellites_.push_back(sats.at(1));
        satellites_.push_back(sats.at(2));
      } else {
        std::ostringstream oss;
        oss << "Incorrect number of satellites to add to inversion!"
            << std::endl;
        throw OpenMDException(oss.str());
      }
    }
    void overrideType(std::string type, std::vector<RealType> pars) {
      orType_      = type;
      orPars_      = pars;
      hasOverride_ = true;
    }

    virtual void validate();
    bool hasOverride() { return hasOverride_; }
    std::string getOverrideType() { return orType_; }

    std::vector<RealType> getOverridePars() { return orPars_; }

  private:
    int center_;
    std::vector<int> satellites_;
    bool hasOverride_;
    std::string orType_;
    std::vector<RealType> orPars_;
  };
}  // namespace OpenMD

#endif
