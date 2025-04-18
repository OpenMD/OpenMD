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

/**
 * @file NonBondedInteractionType.hpp
 * @author    gezelter
 * @date  07/04/2007
 * @version 1.0
 */

#ifndef TYPES_NONBONDEDINTERACTIONTYPE_HPP
#define TYPES_NONBONDEDINTERACTIONTYPE_HPP

#include <memory>

#include "types/AtomType.hpp"

namespace OpenMD {

  typedef struct {
    bool is_LennardJones;
    bool is_Morse;
    bool is_MAW;
    bool is_EAMTable;
    bool is_EAMZhou;
    bool is_EAMOxides;
    bool is_SC;
    bool is_RepulsivePower;
    bool is_Mie;
    bool is_Buckingham;
    bool is_InversePowerSeries;
  } NonBondedInteractionTypeProperties;

  /**
   * @class NonBondedInteractionType
   *
   * NonBondedInteractionType class is responsible for keeping track
   * of static (unchanging) parameters for explicit non-bonded
   * interactions.
   */
  class NonBondedInteractionType {
  public:
    NonBondedInteractionType();
    virtual ~NonBondedInteractionType() {};

    void setLennardJones();
    bool isLennardJones();
    void setMorse();
    bool isMorse();
    void setMAW();
    bool isMAW();
    void setEAMTable();
    bool isEAMTable();
    void setEAMZhou();
    bool isEAMZhou();
    void setEAMOxides();
    bool isEAMOxides();
    bool isSC();
    void setSC();
    bool isMetal();
    void setRepulsivePower();
    bool isRepulsivePower();
    void setMie();
    bool isMie();
    void setBuckingham();
    bool isBuckingham();
    void setInversePowerSeries();
    bool isInversePowerSeries();

    void setAtomTypes(std::pair<AtomType*, AtomType*> ats);
    std::pair<AtomType*, AtomType*> getAtomTypes();

    // below functions are just forward functions
    /**
     * Adds property into property map
     * @param genData GenericData to be added into PropertyMap
     */
    void addProperty(std::shared_ptr<GenericData> genData);

    /**
     * Removes property from PropertyMap by name
     * @param propName the name of property to be removed
     */
    void removeProperty(const std::string& propName);

    /**
     * Returns all names of properties
     * @return all names of properties
     */
    std::vector<std::string> getPropertyNames();

    /**
     * Returns all of the properties in PropertyMap
     * @return all of the properties in PropertyMap
     */
    std::vector<std::shared_ptr<GenericData>> getProperties();

    /**
     * Returns property
     * @param propName name of property
     * @return a pointer point to property with propName. If no
     * property named propName exists, return NULL
     */
    std::shared_ptr<GenericData> getPropertyByName(const std::string& propName);

  protected:
    NonBondedInteractionTypeProperties nbitp;
    std::pair<AtomType*, AtomType*> atomTypes_;

  private:
    // prevent copy construction and copy assignment, since property
    // map contains pointers which can not be copied and managed
    // safely, except make the generic data at PropertyMap as copy on
    // write shared pointer
    NonBondedInteractionType(const NonBondedInteractionType&);
    NonBondedInteractionType& operator=(const NonBondedInteractionType& nbit);
    PropertyMap properties_;
  };
}  // namespace OpenMD

#endif  // TYPES_NONBONDEDINTERACTIONTYPE_HPP
