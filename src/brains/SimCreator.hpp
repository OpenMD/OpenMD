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

/**
 * @file SimCreator.hpp
 * @author tlin
 * @date 11/02/2004
 * @version 1.0
 */

#ifndef BRAINS_SIMCREATOR_HPP
#define BRAINS_SIMCREATOR_HPP

#include <iostream>

#include "brains/ForceField.hpp"
#include "brains/SimInfo.hpp"
#include "io/Globals.hpp"
#include "primitives/Molecule.hpp"

namespace OpenMD {

  /**
   * @class SimCreator SimCreator.hpp "brains/SimCreator.hpp"
   *
   * The only responsibility of SimCreator is to parse the meta-data
   * file and create a SimInfo instance based on the information
   * returned by parser.
   */
  class SimCreator {
  public:
    virtual ~SimCreator() = default;

    /**
     * Setup Simulation
     * @return a pointer to SimInfo
     * @param mdFileName the meta-data file name
     * @param loadInitCoords should the initial coordinates be loaded from a
     * file?
     */
    SimInfo* createSim(const std::string& mdFileName,
                       bool loadInitCoords = true);

  private:
    /**
     * Parses the meta-data file
     * @param mdFileName the meta-data file name
     * @param rawMetaData the raw meta-data stream
     * @param mdFileVersion the version of code used to create the meta-data
     * file
     * @param metaDataStartingLine the starting line of the meta-data block
     * @return a pointer to the simulation parameters in a #Globals object
     */
    Globals* parseFile(std::istream& rawMetaData, const std::string& mdFileName,
                       int mdFileVersion, int metaDataStartingLine);

    /** create the molecules belong to current processor*/
    virtual void createMolecules(SimInfo* info);

    /**
     * Figure out the data storage layout based on what kinds of
     * objects are being simulated
     */
    void computeStorageLayouts(SimInfo* info);

    /**
     * Sets the global index for atoms, rigidbodies and cutoff groups
     * and fill up globalGroupMembership and globalMolMembership
     * arrays which map atoms' global index to the global index of the
     * groups (or molecules) they belong to.  These array are never
     * changed during the simulation.
     */
    void setGlobalIndex(SimInfo* info);

    void gatherParameters(SimInfo* info, const std::string& mdfile);

    /**
     * Divide the molecules among the processors
     */

    void divideMolecules(SimInfo* info);

    /** Load initial coordinates */
    void loadCoordinates(SimInfo* info, const std::string& mdFileName);

    std::string
        mdFileName_;  // save the meta-data file name which may be used later
  };
}  // namespace OpenMD

#endif  // BRAINS_SIMCREATOR_HPP
