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
 
 /**
  * @file SimCreatorr.hpp
  * @author tlin
  * @date 11/02/2004
  * @time 12:126am
  * @version 1.0
  */

#ifndef BRAINS_SIMCREATOR_HPP
#define BRAINS_SIMCREATOR_HPP


#include "primitives/Molecule.hpp"
#include "brains/SimInfo.hpp"
#include "types/MakeStamps.hpp"
#include "io/Globals.hpp"
#include "UseTheForce/ForceField.hpp"

// this routine is defined in BASS_interface.cpp
//another OOPS
extern void set_interface_stamps( MakeStamps* ms, Globals* g );

namespace oopse {

/**
 * @class SimCreator SimCreator.hpp "brains/SimCreator.hpp"
 * The only responsibility of SimCreator is to parse the meta-data file and create a SimInfo
 * instance based on the information returned by parser. 
 */
class SimCreator {
    public:

        /**
         * Setup Simulation
         * @return a pointer to SimInfo
         * @param mdfile the meta-data file name
         */
        SimInfo* createSim(const std::string & mdFileName, bool loadInitCoords = true);
        
    private:
        
        /**
         * Parses the meta-data file
         * @param mdfile
         * @param stamps
         * @param simParams
         */
        void parseFile(const std::string mdFileName,  MakeStamps* stamps, Globals* simParams);


        /** create the molecules belong to current processor*/
        virtual void createMolecules(SimInfo* info);

        /** 
         * Sets the global index for atoms, rigidbodies and cutoff groups and fill up
         * globalGroupMembership and globalMolMembership arrays which map atoms'
         * global index to the global index of the groups (or molecules) they belong to.
         * These array are never changed during the simulation.
         */
        void setGlobalIndex(SimInfo* info);

        void gatherParameters(SimInfo *info, const std::string& mdfile);             

        
        /** Extracts the molecules stamps and adds them into SimInfo class */
        void compList(MakeStamps* stamps,  Globals* simParams, 
                                     std::vector<std::pair<MoleculeStamp*, int> >& moleculeStamps) ;

        /**
         * Divide the molecules among the processors 
         */
         
        void divideMolecules(SimInfo* info);

        /** Load initial coordinates */
        void loadCoordinates(SimInfo* info);     

        std::string mdFileName_;  //save the meta-data file name which may be used later
};

} //end namespace oopse
#endif //BRAINS_SIMCREATOR_HPP

