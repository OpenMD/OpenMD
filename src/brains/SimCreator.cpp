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
 * @file SimCreator.cpp
 * @author tlin
 * @date 11/03/2004
 * @time 13:51am
 * @version 1.0
 */

#include <sprng.h>

#include "brains/MoleculeCreator.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimSnapshotManager.hpp"
#include "io/DumpReader.hpp"
#include "io/parse_me.h"
#include "UseTheForce/ForceFieldFactory.hpp"
#include "utils/simError.h"
#include "utils/StringUtils.hpp"
#ifdef IS_MPI
#include "io/mpiBASS.h"
#include "math/randomSPRNG.hpp"
#endif

namespace oopse {

void SimCreator::parseFile(const std::string mdFileName,  MakeStamps* stamps, Globals* simParams){

#ifdef IS_MPI

    if (worldRank == 0) {
#endif // is_mpi

        simParams->initalize();
        set_interface_stamps(stamps, simParams);

#ifdef IS_MPI

        mpiEventInit();

#endif

        yacc_BASS(mdFileName.c_str());

#ifdef IS_MPI

        throwMPIEvent(NULL);
    } else {
        set_interface_stamps(stamps, simParams);
        mpiEventInit();
        MPIcheckPoint();
        mpiEventLoop();
    }

#endif

}

SimInfo*  SimCreator::createSim(const std::string & mdFileName, bool loadInitCoords) {
    
    MakeStamps * stamps = new MakeStamps();

    Globals * simParams = new Globals();

    //parse meta-data file
    parseFile(mdFileName, stamps, simParams);

    //create the force field
    ForceField * ff = ForceFieldFactory::getInstance()->createForceField(
                          simParams->getForceField());
    
    if (ff == NULL) {
        sprintf(painCave.errMsg, "ForceField Factory can not create %s force field\n",
                simParams->getForceField());
        painCave.isFatal = 1;
        simError();
    }

    if (simParams->haveForceFieldFileName()) {
        ff->setForceFieldFileName(simParams->getForceFieldFileName());
    }
    
    std::string forcefieldFileName;
    forcefieldFileName = ff->getForceFieldFileName();

    if (simParams->haveForceFieldVariant()) {
        //If the force field has variant, the variant force field name will be
        //Base.variant.frc. For exampel EAM.u6.frc
        
        std::string variant = simParams->getForceFieldVariant();

        std::string::size_type pos = forcefieldFileName.rfind(".frc");
        variant = "." + variant;
        if (pos != std::string::npos) {
            forcefieldFileName.insert(pos, variant);
        } else {
            //If the default force field file name does not containt .frc suffix, just append the .variant
            forcefieldFileName.append(variant);
        }
    } 
    
    ff->parse(forcefieldFileName);
    
    //extract the molecule stamps
    std::vector < std::pair<MoleculeStamp *, int> > moleculeStampPairs;
    compList(stamps, simParams, moleculeStampPairs);

    //create SimInfo
    SimInfo * info = new SimInfo(moleculeStampPairs, ff, simParams);

    //gather parameters (SimCreator only retrieves part of the parameters)
    gatherParameters(info, mdFileName);

    //divide the molecules and determine the global index of molecules
#ifdef IS_MPI
    divideMolecules(info);
#endif 

    //create the molecules
    createMolecules(info);


    //allocate memory for DataStorage(circular reference, need to break it)
    info->setSnapshotManager(new SimSnapshotManager(info));
    
    //set the global index of atoms, rigidbodies and cutoffgroups (only need to be set once, the
    //global index will never change again). Local indices of atoms and rigidbodies are already set by 
    //MoleculeCreator class which actually delegates the responsibility to LocalIndexManager. 
    setGlobalIndex(info);

    //Alought addExculdePairs is called inside SimInfo's addMolecule method, at that point
    //atoms don't have the global index yet  (their global index are all initialized to -1).
    //Therefore we have to call addExcludePairs explicitly here. A way to work around is that
    //we can determine the beginning global indices of atoms before they get created.
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    for (mol= info->beginMolecule(mi); mol != NULL; mol = info->nextMolecule(mi)) {
        info->addExcludePairs(mol);
    }
    

    //load initial coordinates, some extra information are pushed into SimInfo's property map ( such as
    //eta, chi for NPT integrator)
    if (loadInitCoords)
        loadCoordinates(info);    
    
    return info;
}

void SimCreator::gatherParameters(SimInfo *info, const std::string& mdfile) {

    //setup seed for random number generator
    int seedValue;
    Globals * simParams = info->getSimParams();

    if (simParams->haveSeed()) {
        seedValue = simParams->getSeed();

        if (seedValue < 100000000 ) {
            sprintf(painCave.errMsg,
                    "Seed for sprng library should contain at least 9 digits\n"
                        "OOPSE will generate a seed for user\n");

            painCave.isFatal = 0;
            simError();

            //using seed generated by system instead of invalid seed set by user 

#ifndef IS_MPI

            seedValue = make_sprng_seed();

#else

            if (worldRank == 0) {
                seedValue = make_sprng_seed();
            }

            MPI_Bcast(&seedValue, 1, MPI_INT, 0, MPI_COMM_WORLD);

#endif

        } //end if (seedValue /1000000000 == 0)
    } else {

#ifndef IS_MPI

        seedValue = make_sprng_seed();

#else

        if (worldRank == 0) {
            seedValue = make_sprng_seed();
        }

        MPI_Bcast(&seedValue, 1, MPI_INT, 0, MPI_COMM_WORLD);

#endif

    } //end of simParams->haveSeed()

    info->setSeed(seedValue);


    //figure out the ouput file names
    std::string prefix;

#ifdef IS_MPI

    if (worldRank == 0) {
#endif // is_mpi

        if (simParams->haveFinalConfig()) {
            prefix = getPrefix(simParams->getFinalConfig());
        } else {
            prefix = getPrefix(mdfile);
        }

        info->setFinalConfigFileName(prefix + ".eor");
        info->setDumpFileName(prefix + ".dump");
        info->setStatFileName(prefix + ".stat");

#ifdef IS_MPI

    }

#endif

}

#ifdef IS_MPI
void SimCreator::divideMolecules(SimInfo *info) {
    double numerator;
    double denominator;
    double precast;
    double x;
    double y;
    double a;
    int old_atoms;
    int add_atoms;
    int new_atoms;
    int nTarget;
    int done;
    int i;
    int j;
    int loops;
    int which_proc;
    int nProcessors;
    std::vector<int> atomsPerProc;
    randomSPRNG myRandom(info->getSeed());
    int nGlobalMols = info->getNGlobalMolecules();
    std::vector<int> molToProcMap(nGlobalMols, -1); // default to an error condition:
    
    MPI_Comm_size(MPI_COMM_WORLD, &nProcessors);

    if (nProcessors > nGlobalMols) {
        sprintf(painCave.errMsg,
                "nProcessors (%d) > nMol (%d)\n"
                    "\tThe number of processors is larger than\n"
                    "\tthe number of molecules.  This will not result in a \n"
                    "\tusable division of atoms for force decomposition.\n"
                    "\tEither try a smaller number of processors, or run the\n"
                    "\tsingle-processor version of OOPSE.\n", nProcessors, nGlobalMols);

        painCave.isFatal = 1;
        simError();
    }

    a = 3.0 * nGlobalMols / info->getNGlobalAtoms();

    //initialize atomsPerProc
    atomsPerProc.insert(atomsPerProc.end(), nProcessors, 0);

    if (worldRank == 0) {
        numerator = info->getNGlobalAtoms();
        denominator = nProcessors;
        precast = numerator / denominator;
        nTarget = (int)(precast + 0.5);

        for(i = 0; i < nGlobalMols; i++) {
            done = 0;
            loops = 0;

            while (!done) {
                loops++;

                // Pick a processor at random

                which_proc = (int) (myRandom.getRandom() * nProcessors);

                //get the molecule stamp first
                int stampId = info->getMoleculeStampId(i);
                MoleculeStamp * moleculeStamp = info->getMoleculeStamp(stampId);

                // How many atoms does this processor have so far?
                old_atoms = atomsPerProc[which_proc];
                add_atoms = moleculeStamp->getNAtoms();
                new_atoms = old_atoms + add_atoms;

                // If we've been through this loop too many times, we need
                // to just give up and assign the molecule to this processor
                // and be done with it. 

                if (loops > 100) {
                    sprintf(painCave.errMsg,
                            "I've tried 100 times to assign molecule %d to a "
                                " processor, but can't find a good spot.\n"
                                "I'm assigning it at random to processor %d.\n",
                            i, which_proc);

                    painCave.isFatal = 0;
                    simError();

                    molToProcMap[i] = which_proc;
                    atomsPerProc[which_proc] += add_atoms;

                    done = 1;
                    continue;
                }

                // If we can add this molecule to this processor without sending
                // it above nTarget, then go ahead and do it:

                if (new_atoms <= nTarget) {
                    molToProcMap[i] = which_proc;
                    atomsPerProc[which_proc] += add_atoms;

                    done = 1;
                    continue;
                }

                // The only situation left is when new_atoms > nTarget.  We
                // want to accept this with some probability that dies off the
                // farther we are from nTarget

                // roughly:  x = new_atoms - nTarget
                //           Pacc(x) = exp(- a * x)
                // where a = penalty / (average atoms per molecule)

                x = (double)(new_atoms - nTarget);
                y = myRandom.getRandom();

                if (y < exp(- a * x)) {
                    molToProcMap[i] = which_proc;
                    atomsPerProc[which_proc] += add_atoms;

                    done = 1;
                    continue;
                } else {
                    continue;
                }
            }
        }

        // Spray out this nonsense to all other processors:

        MPI_Bcast(&molToProcMap[0], nGlobalMols, MPI_INT, 0, MPI_COMM_WORLD);
    } else {

        // Listen to your marching orders from processor 0:

        MPI_Bcast(&molToProcMap[0], nGlobalMols, MPI_INT, 0, MPI_COMM_WORLD);
    }

    info->setMolToProcMap(molToProcMap);
    sprintf(checkPointMsg,
            "Successfully divided the molecules among the processors.\n");
    MPIcheckPoint();
}

#endif

void SimCreator::createMolecules(SimInfo *info) {
    MoleculeCreator molCreator;
    int stampId;

    for(int i = 0; i < info->getNGlobalMolecules(); i++) {

#ifdef IS_MPI

        if (info->getMolToProc(i) == worldRank) {
#endif

            stampId = info->getMoleculeStampId(i);
            Molecule * mol = molCreator.createMolecule(info->getForceField(), info->getMoleculeStamp(stampId),
                                                                                    stampId, i, info->getLocalIndexManager());

            info->addMolecule(mol);

#ifdef IS_MPI

        }

#endif

    } //end for(int i=0)   
}

void SimCreator::compList(MakeStamps *stamps, Globals* simParams,
                        std::vector < std::pair<MoleculeStamp *, int> > &moleculeStampPairs) {
    int i;
    char * id;
    MoleculeStamp * currentStamp;
    Component** the_components = simParams->getComponents();
    int n_components = simParams->getNComponents();

    if (!simParams->haveNMol()) {
        // we don't have the total number of molecules, so we assume it is
        // given in each component

        for(i = 0; i < n_components; i++) {
            if (!the_components[i]->haveNMol()) {
                // we have a problem
                sprintf(painCave.errMsg,
                        "SimCreator Error. No global NMol or component NMol given.\n"
                            "\tCannot calculate the number of atoms.\n");

                painCave.isFatal = 1;
                simError();
            }

            id = the_components[i]->getType();
            currentStamp = (stamps->extractMolStamp(id))->getStamp();

            if (currentStamp == NULL) {
                sprintf(painCave.errMsg,
                        "SimCreator error: Component \"%s\" was not found in the "
                            "list of declared molecules\n", id);

                painCave.isFatal = 1;
                simError();
            }

            moleculeStampPairs.push_back(
                std::make_pair(currentStamp, the_components[i]->getNMol()));
        } //end for (i = 0; i < n_components; i++)
    } else {
        sprintf(painCave.errMsg, "SimSetup error.\n"
                                     "\tSorry, the ability to specify total"
                                     " nMols and then give molfractions in the components\n"
                                     "\tis not currently supported."
                                     " Please give nMol in the components.\n");

        painCave.isFatal = 1;
        simError();
    }

#ifdef IS_MPI

    strcpy(checkPointMsg, "Component stamps successfully extracted\n");
    MPIcheckPoint();

#endif // is_mpi

}

void SimCreator::setGlobalIndex(SimInfo *info) {
    SimInfo::MoleculeIterator mi;
    Molecule::AtomIterator ai;
    Molecule::RigidBodyIterator ri;
    Molecule::CutoffGroupIterator ci;
    Molecule * mol;
    Atom * atom;
    RigidBody * rb;
    CutoffGroup * cg;
    int beginAtomIndex;
    int beginRigidBodyIndex;
    int beginCutoffGroupIndex;
    int nGlobalAtoms = info->getNGlobalAtoms();
    
#ifndef IS_MPI

    beginAtomIndex = 0;
    beginRigidBodyIndex = 0;
    beginCutoffGroupIndex = 0;

#else

    int nproc;
    int myNode;

    myNode = worldRank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    std::vector < int > tmpAtomsInProc(nproc, 0);
    std::vector < int > tmpRigidBodiesInProc(nproc, 0);
    std::vector < int > tmpCutoffGroupsInProc(nproc, 0);
    std::vector < int > NumAtomsInProc(nproc, 0);
    std::vector < int > NumRigidBodiesInProc(nproc, 0);
    std::vector < int > NumCutoffGroupsInProc(nproc, 0);

    tmpAtomsInProc[myNode] = info->getNAtoms();
    tmpRigidBodiesInProc[myNode] = info->getNRigidBodies();
    tmpCutoffGroupsInProc[myNode] = info->getNCutoffGroups();

    //do MPI_ALLREDUCE to exchange the total number of atoms, rigidbodies and cutoff groups
    MPI_Allreduce(&tmpAtomsInProc[0], &NumAtomsInProc[0], nproc, MPI_INT,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&tmpRigidBodiesInProc[0], &NumRigidBodiesInProc[0], nproc,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&tmpCutoffGroupsInProc[0], &NumCutoffGroupsInProc[0], nproc,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    beginAtomIndex = 0;
    beginRigidBodyIndex = 0;
    beginCutoffGroupIndex = 0;

    for(int i = 0; i < myNode; i++) {
        beginAtomIndex += NumAtomsInProc[i];
        beginRigidBodyIndex += NumRigidBodiesInProc[i];
        beginCutoffGroupIndex += NumCutoffGroupsInProc[i];
    }

#endif

    for(mol = info->beginMolecule(mi); mol != NULL;
        mol = info->nextMolecule(mi)) {

        //local index(index in DataStorge) of atom is important
        for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
            atom->setGlobalIndex(beginAtomIndex++);
        }

        for(rb = mol->beginRigidBody(ri); rb != NULL;
            rb = mol->nextRigidBody(ri)) {
            rb->setGlobalIndex(beginRigidBodyIndex++);
        }

        //local index of cutoff group is trivial, it only depends on the order of travesing
        for(cg = mol->beginCutoffGroup(ci); cg != NULL;
            cg = mol->nextCutoffGroup(ci)) {
            cg->setGlobalIndex(beginCutoffGroupIndex++);
        }
    }

    //fill globalGroupMembership
    std::vector<int> globalGroupMembership(info->getNGlobalAtoms(), 0);
    for(mol = info->beginMolecule(mi); mol != NULL; mol = info->nextMolecule(mi)) {        
        for (cg = mol->beginCutoffGroup(ci); cg != NULL; cg = mol->nextCutoffGroup(ci)) {

            for(atom = cg->beginAtom(ai); atom != NULL; atom = cg->nextAtom(ai)) {
                globalGroupMembership[atom->getGlobalIndex()] = cg->getGlobalIndex();
            }

        }       
    }

#ifdef IS_MPI    
    // Since the globalGroupMembership has been zero filled and we've only
    // poked values into the atoms we know, we can do an Allreduce
    // to get the full globalGroupMembership array (We think).
    // This would be prettier if we could use MPI_IN_PLACE like the MPI-2
    // docs said we could.
    std::vector<int> tmpGroupMembership(nGlobalAtoms, 0);
    MPI_Allreduce(&globalGroupMembership[0], &tmpGroupMembership[0], nGlobalAtoms,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
     info->setGlobalGroupMembership(tmpGroupMembership);
#else
    info->setGlobalGroupMembership(globalGroupMembership);
#endif

    //fill molMembership
    std::vector<int> globalMolMembership(info->getNGlobalAtoms(), 0);
    
    for(mol = info->beginMolecule(mi); mol != NULL; mol = info->nextMolecule(mi)) {

        for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
            globalMolMembership[atom->getGlobalIndex()] = mol->getGlobalIndex();
        }
    }

#ifdef IS_MPI
    std::vector<int> tmpMolMembership(nGlobalAtoms, 0);

    MPI_Allreduce(&globalMolMembership[0], &tmpMolMembership[0], nGlobalAtoms,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    info->setGlobalMolMembership(tmpMolMembership);
#else
    info->setGlobalMolMembership(globalMolMembership);
#endif

}

void SimCreator::loadCoordinates(SimInfo* info) {
    Globals* simParams;
    simParams = info->getSimParams();
    
    if (!simParams->haveInitialConfig()) {
        sprintf(painCave.errMsg,
                "Cannot intialize a simulation without an initial configuration file.\n");
        painCave.isFatal = 1;;
        simError();
    }
        
    DumpReader reader(info, simParams->getInitialConfig());
    int nframes = reader.getNFrames();

    if (nframes > 0) {
        reader.readFrame(nframes - 1);
    } else {
        //invalid initial coordinate file
        sprintf(painCave.errMsg, "Initial configuration file %s should at least contain one frame\n",
                simParams->getInitialConfig());
        painCave.isFatal = 1;
        simError();
    }

    //copy the current snapshot to previous snapshot
    info->getSnapshotManager()->advance();
}

} //end namespace oopse


