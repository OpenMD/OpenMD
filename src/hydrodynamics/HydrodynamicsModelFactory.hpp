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
 * @file HydrodynamicsModelFactory.hpp
 * @author Teng Lin
 * @date 10/24/2004
 * @version 1.0
 */

#ifndef HYDRODYNAMICS_HYDRODYNAMICSMODELFACTORY_HPP
#define HYDRODYNAMICS_HYDRODYNAMICSMODELFACTORY_HPP

#include <cassert>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace OpenMD {

  // forward declaration
  class HydrodynamicsModel;
  class HydrodynamicsModelCreator;

  /**
   * @class HydrodynamicsModelFactory
   * Factory pattern and Singleton Pattern are used to define an interface for
   * creating an HydrodynamicsModel.
   */
  class HydrodynamicsModelFactory {
  public:
    using CreatorMapType  = std::map<std::string, HydrodynamicsModelCreator*>;
    using IdentVectorType = std::vector<std::string>;
    using IdentVectorIterator = std::vector<std::string>::iterator;

    ~HydrodynamicsModelFactory();

    /**
     * Returns an instance of HydrodynamicsModel factory
     * @return an instance of HydrodynamicsModel factory
     */
    static HydrodynamicsModelFactory* getInstance() {
      if (instance_ == NULL) { instance_ = new HydrodynamicsModelFactory(); }
      return instance_;
    }

    /**
     * Registers a creator with a type identifier
     * @return true if registration is successful, otherwise return false
     * @param creator the object responsible to create the concrete object
     */
    bool registerHydrodynamicsModel(HydrodynamicsModelCreator* creator);

    /**
     * Unregisters the creator for the given type identifier. If the type
     * identifier was previously registered, the function returns true.
     * @return truethe type identifier was previously registered and the creator
     * is removed, otherwise return false
     * @param id the identification of the concrete object
     */
    bool unregisterHydrodynamicsModel(const std::string& id);
    /**
     * Looks up the type identifier in the internal map. If it is found, it
     * invokes the corresponding creator for the type identifier and returns its
     * result.
     * @return a pointer of the concrete object, return NULL if no creator is
     * registed for creating this concrete object
     * @param id the identification of the concrete object
     * @param sd a pointer to the StuntDouble being modeled
     * @param info a pointer to the SimInfo object
     */
    HydrodynamicsModel* createHydrodynamicsModel(const std::string& id);

    /**
     *  Returns all of the registed  type identifiers
     * @return all of the registed  type identifiers
     */
    IdentVectorType getIdents();

  private:
    HydrodynamicsModelFactory() {}
    static HydrodynamicsModelFactory* instance_;
    CreatorMapType creatorMap_;
  };

  /** write out all of the type identifiers to an output stream */
  std::ostream& operator<<(std::ostream& o, HydrodynamicsModelFactory& factory);

}  // namespace OpenMD

#endif
