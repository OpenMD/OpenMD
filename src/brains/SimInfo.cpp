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
 * @file SimInfo.cpp
 * @author    tlin
 * @date  11/02/2004
 * @version 1.0
 */

#include <algorithm>
#include <set>

#include "brains/SimInfo.hpp"
#include "math/Vector3.hpp"
#include "primitives/Molecule.hpp"
#include "UseTheForce/doForces_interface.h"
#include "UseTheForce/notifyCutoffs_interface.h"
#include "utils/MemoryUtils.hpp"
#include "utils/simError.h"
#include "selection/SelectionManager.hpp"

#ifdef IS_MPI
#include "UseTheForce/mpiComponentPlan.h"
#include "UseTheForce/DarkSide/simParallel_interface.h"
#endif 

namespace oopse {

SimInfo::SimInfo(std::vector<std::pair<MoleculeStamp*, int> >& molStampPairs, 
                                ForceField* ff, Globals* simParams) : 
                                forceField_(ff), simParams_(simParams), 
                                ndf_(0), ndfRaw_(0), ndfTrans_(0), nZconstraint_(0),
                                nGlobalMols_(0), nGlobalAtoms_(0), nGlobalCutoffGroups_(0), 
                                nGlobalIntegrableObjects_(0), nGlobalRigidBodies_(0),
                                nAtoms_(0), nBonds_(0),  nBends_(0), nTorsions_(0), nRigidBodies_(0),
                                nIntegrableObjects_(0),  nCutoffGroups_(0), nConstraints_(0),
                                sman_(NULL), fortranInitialized_(false), selectMan_(NULL) {

            
    std::vector<std::pair<MoleculeStamp*, int> >::iterator i;
    MoleculeStamp* molStamp;
    int nMolWithSameStamp;
    int nCutoffAtoms = 0; // number of atoms belong to cutoff groups
    int nGroups = 0;          //total cutoff groups defined in meta-data file
    CutoffGroupStamp* cgStamp;    
    RigidBodyStamp* rbStamp;
    int nRigidAtoms = 0;
    
    for (i = molStampPairs.begin(); i !=molStampPairs.end(); ++i) {
        molStamp = i->first;
        nMolWithSameStamp = i->second;
        
        addMoleculeStamp(molStamp, nMolWithSameStamp);

        //calculate atoms in molecules
        nGlobalAtoms_ += molStamp->getNAtoms() *nMolWithSameStamp;   


        //calculate atoms in cutoff groups
        int nAtomsInGroups = 0;
        int nCutoffGroupsInStamp = molStamp->getNCutoffGroups();
        
        for (int j=0; j < nCutoffGroupsInStamp; j++) {
            cgStamp = molStamp->getCutoffGroup(j);
            nAtomsInGroups += cgStamp->getNMembers();
        }

        nGroups += nCutoffGroupsInStamp * nMolWithSameStamp;
        nCutoffAtoms += nAtomsInGroups * nMolWithSameStamp;            

        //calculate atoms in rigid bodies
        int nAtomsInRigidBodies = 0;
        int nRigidBodiesInStamp = molStamp->getNRigidBodies();
        
        for (int j=0; j < nRigidBodiesInStamp; j++) {
            rbStamp = molStamp->getRigidBody(j);
            nAtomsInRigidBodies += rbStamp->getNMembers();
        }

        nGlobalRigidBodies_ += nRigidBodiesInStamp * nMolWithSameStamp;
        nRigidAtoms += nAtomsInRigidBodies * nMolWithSameStamp;            
        
    }

    //every free atom (atom does not belong to cutoff groups) is a cutoff group
    //therefore the total number of cutoff groups in the system is equal to 
    //the total number of atoms minus number of atoms belong to cutoff group defined in meta-data
    //file plus the number of cutoff groups defined in meta-data file
    nGlobalCutoffGroups_ = nGlobalAtoms_ - nCutoffAtoms + nGroups;

    //every free atom (atom does not belong to rigid bodies) is an integrable object
    //therefore the total number of  integrable objects in the system is equal to 
    //the total number of atoms minus number of atoms belong to  rigid body defined in meta-data
    //file plus the number of  rigid bodies defined in meta-data file
    nGlobalIntegrableObjects_ = nGlobalAtoms_ - nRigidAtoms + nGlobalRigidBodies_;

    nGlobalMols_ = molStampIds_.size();

#ifdef IS_MPI    
    molToProcMap_.resize(nGlobalMols_);
#endif

    selectMan_ = new SelectionManager(this);
    selectMan_->selectAll();
}

SimInfo::~SimInfo() {
    //MemoryUtils::deleteVectorOfPointer(molecules_);

    MemoryUtils::deleteVectorOfPointer(moleculeStamps_);
    
    delete sman_;
    delete simParams_;
    delete forceField_;
    delete selectMan_;
}

int SimInfo::getNGlobalConstraints() {
    int nGlobalConstraints;
#ifdef IS_MPI
    MPI_Allreduce(&nConstraints_, &nGlobalConstraints, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);    
#else
    nGlobalConstraints =  nConstraints_;
#endif
    return nGlobalConstraints;
}

bool SimInfo::addMolecule(Molecule* mol) {
    MoleculeIterator i;

    i = molecules_.find(mol->getGlobalIndex());
    if (i == molecules_.end() ) {

        molecules_.insert(std::make_pair(mol->getGlobalIndex(), mol));
        
        nAtoms_ += mol->getNAtoms();
        nBonds_ += mol->getNBonds();
        nBends_ += mol->getNBends();
        nTorsions_ += mol->getNTorsions();
        nRigidBodies_ += mol->getNRigidBodies();
        nIntegrableObjects_ += mol->getNIntegrableObjects();
        nCutoffGroups_ += mol->getNCutoffGroups();
        nConstraints_ += mol->getNConstraintPairs();

        addExcludePairs(mol);
        
        return true;
    } else {
        return false;
    }
}

bool SimInfo::removeMolecule(Molecule* mol) {
    MoleculeIterator i;
    i = molecules_.find(mol->getGlobalIndex());

    if (i != molecules_.end() ) {

        assert(mol == i->second);
        
        nAtoms_ -= mol->getNAtoms();
        nBonds_ -= mol->getNBonds();
        nBends_ -= mol->getNBends();
        nTorsions_ -= mol->getNTorsions();
        nRigidBodies_ -= mol->getNRigidBodies();
        nIntegrableObjects_ -= mol->getNIntegrableObjects();
        nCutoffGroups_ -= mol->getNCutoffGroups();
        nConstraints_ -= mol->getNConstraintPairs();

        removeExcludePairs(mol);
        molecules_.erase(mol->getGlobalIndex());

        delete mol;
        
        return true;
    } else {
        return false;
    }


}    

        
Molecule* SimInfo::beginMolecule(MoleculeIterator& i) {
    i = molecules_.begin();
    return i == molecules_.end() ? NULL : i->second;
}    

Molecule* SimInfo::nextMolecule(MoleculeIterator& i) {
    ++i;
    return i == molecules_.end() ? NULL : i->second;    
}


void SimInfo::calcNdf() {
    int ndf_local;
    MoleculeIterator i;
    std::vector<StuntDouble*>::iterator j;
    Molecule* mol;
    StuntDouble* integrableObject;

    ndf_local = 0;
    
    for (mol = beginMolecule(i); mol != NULL; mol = nextMolecule(i)) {
        for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL; 
               integrableObject = mol->nextIntegrableObject(j)) {

            ndf_local += 3;

            if (integrableObject->isDirectional()) {
                if (integrableObject->isLinear()) {
                    ndf_local += 2;
                } else {
                    ndf_local += 3;
                }
            }
            
        }//end for (integrableObject)
    }// end for (mol)
    
    // n_constraints is local, so subtract them on each processor
    ndf_local -= nConstraints_;

#ifdef IS_MPI
    MPI_Allreduce(&ndf_local,&ndf_,1,MPI_INT,MPI_SUM, MPI_COMM_WORLD);
#else
    ndf_ = ndf_local;
#endif

    // nZconstraints_ is global, as are the 3 COM translations for the 
    // entire system:
    ndf_ = ndf_ - 3 - nZconstraint_;

}

void SimInfo::calcNdfRaw() {
    int ndfRaw_local;

    MoleculeIterator i;
    std::vector<StuntDouble*>::iterator j;
    Molecule* mol;
    StuntDouble* integrableObject;

    // Raw degrees of freedom that we have to set
    ndfRaw_local = 0;
    
    for (mol = beginMolecule(i); mol != NULL; mol = nextMolecule(i)) {
        for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL;
               integrableObject = mol->nextIntegrableObject(j)) {

            ndfRaw_local += 3;

            if (integrableObject->isDirectional()) {
                if (integrableObject->isLinear()) {
                    ndfRaw_local += 2;
                } else {
                    ndfRaw_local += 3;
                }
            }
            
        }
    }
    
#ifdef IS_MPI
    MPI_Allreduce(&ndfRaw_local,&ndfRaw_,1,MPI_INT,MPI_SUM, MPI_COMM_WORLD);
#else
    ndfRaw_ = ndfRaw_local;
#endif
}

void SimInfo::calcNdfTrans() {
    int ndfTrans_local;

    ndfTrans_local = 3 * nIntegrableObjects_ - nConstraints_;


#ifdef IS_MPI
    MPI_Allreduce(&ndfTrans_local,&ndfTrans_,1,MPI_INT,MPI_SUM, MPI_COMM_WORLD);
#else
    ndfTrans_ = ndfTrans_local;
#endif

    ndfTrans_ = ndfTrans_ - 3 - nZconstraint_;
 
}

void SimInfo::addExcludePairs(Molecule* mol) {
    std::vector<Bond*>::iterator bondIter;
    std::vector<Bend*>::iterator bendIter;
    std::vector<Torsion*>::iterator torsionIter;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    int a;
    int b;
    int c;
    int d;
    
    for (bond= mol->beginBond(bondIter); bond != NULL; bond = mol->nextBond(bondIter)) {
        a = bond->getAtomA()->getGlobalIndex();
        b = bond->getAtomB()->getGlobalIndex();        
        exclude_.addPair(a, b);
    }

    for (bend= mol->beginBend(bendIter); bend != NULL; bend = mol->nextBend(bendIter)) {
        a = bend->getAtomA()->getGlobalIndex();
        b = bend->getAtomB()->getGlobalIndex();        
        c = bend->getAtomC()->getGlobalIndex();

        exclude_.addPair(a, b);
        exclude_.addPair(a, c);
        exclude_.addPair(b, c);        
    }

    for (torsion= mol->beginTorsion(torsionIter); torsion != NULL; torsion = mol->nextTorsion(torsionIter)) {
        a = torsion->getAtomA()->getGlobalIndex();
        b = torsion->getAtomB()->getGlobalIndex();        
        c = torsion->getAtomC()->getGlobalIndex();        
        d = torsion->getAtomD()->getGlobalIndex();        

        exclude_.addPair(a, b);
        exclude_.addPair(a, c);
        exclude_.addPair(a, d);
        exclude_.addPair(b, c);
        exclude_.addPair(b, d);
        exclude_.addPair(c, d);        
    }

    
}

void SimInfo::removeExcludePairs(Molecule* mol) {
    std::vector<Bond*>::iterator bondIter;
    std::vector<Bend*>::iterator bendIter;
    std::vector<Torsion*>::iterator torsionIter;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    int a;
    int b;
    int c;
    int d;
    
    for (bond= mol->beginBond(bondIter); bond != NULL; bond = mol->nextBond(bondIter)) {
        a = bond->getAtomA()->getGlobalIndex();
        b = bond->getAtomB()->getGlobalIndex();        
        exclude_.removePair(a, b);
    }

    for (bend= mol->beginBend(bendIter); bend != NULL; bend = mol->nextBend(bendIter)) {
        a = bend->getAtomA()->getGlobalIndex();
        b = bend->getAtomB()->getGlobalIndex();        
        c = bend->getAtomC()->getGlobalIndex();

        exclude_.removePair(a, b);
        exclude_.removePair(a, c);
        exclude_.removePair(b, c);        
    }

    for (torsion= mol->beginTorsion(torsionIter); torsion != NULL; torsion = mol->nextTorsion(torsionIter)) {
        a = torsion->getAtomA()->getGlobalIndex();
        b = torsion->getAtomB()->getGlobalIndex();        
        c = torsion->getAtomC()->getGlobalIndex();        
        d = torsion->getAtomD()->getGlobalIndex();        

        exclude_.removePair(a, b);
        exclude_.removePair(a, c);
        exclude_.removePair(a, d);
        exclude_.removePair(b, c);
        exclude_.removePair(b, d);
        exclude_.removePair(c, d);        
    }

}


void SimInfo::addMoleculeStamp(MoleculeStamp* molStamp, int nmol) {
    int curStampId;

    //index from 0
    curStampId = moleculeStamps_.size();

    moleculeStamps_.push_back(molStamp);
    molStampIds_.insert(molStampIds_.end(), nmol, curStampId);
}

void SimInfo::update() {

    setupSimType();

#ifdef IS_MPI
    setupFortranParallel();
#endif

    setupFortranSim();

    //setup fortran force field
    /** @deprecate */    
    int isError = 0;
    initFortranFF( &fInfo_.SIM_uses_RF , &isError );
    if(isError){
        sprintf( painCave.errMsg,
         "ForceField error: There was an error initializing the forceField in fortran.\n" );
        painCave.isFatal = 1;
        simError();
    }
  
    
    setupCutoff();

    calcNdf();
    calcNdfRaw();
    calcNdfTrans();

    fortranInitialized_ = true;
}

std::set<AtomType*> SimInfo::getUniqueAtomTypes() {
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::AtomIterator ai;
    Atom* atom;
    std::set<AtomType*> atomTypes;

    for(mol = beginMolecule(mi); mol != NULL; mol = nextMolecule(mi)) {

        for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
            atomTypes.insert(atom->getAtomType());
        }
        
    }

    return atomTypes;        
}

void SimInfo::setupSimType() {
    std::set<AtomType*>::iterator i;
    std::set<AtomType*> atomTypes;
    atomTypes = getUniqueAtomTypes();
    
    int useLennardJones = 0;
    int useElectrostatic = 0;
    int useEAM = 0;
    int useCharge = 0;
    int useDirectional = 0;
    int useDipole = 0;
    int useGayBerne = 0;
    int useSticky = 0;
    int useShape = 0; 
    int useFLARB = 0; //it is not in AtomType yet
    int useDirectionalAtom = 0;    
    int useElectrostatics = 0;
    //usePBC and useRF are from simParams
    int usePBC = simParams_->getPBC();
    int useRF = simParams_->getUseRF();

    //loop over all of the atom types
    for (i = atomTypes.begin(); i != atomTypes.end(); ++i) {
        useLennardJones |= (*i)->isLennardJones();
        useElectrostatic |= (*i)->isElectrostatic();
        useEAM |= (*i)->isEAM();
        useCharge |= (*i)->isCharge();
        useDirectional |= (*i)->isDirectional();
        useDipole |= (*i)->isDipole();
        useGayBerne |= (*i)->isGayBerne();
        useSticky |= (*i)->isSticky();
        useShape |= (*i)->isShape(); 
    }

    if (useSticky || useDipole || useGayBerne || useShape) {
        useDirectionalAtom = 1;
    }

    if (useCharge || useDipole) {
        useElectrostatics = 1;
    }

#ifdef IS_MPI    
    int temp;

    temp = usePBC;
    MPI_Allreduce(&temp, &usePBC, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);    

    temp = useDirectionalAtom;
    MPI_Allreduce(&temp, &useDirectionalAtom, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);    

    temp = useLennardJones;
    MPI_Allreduce(&temp, &useLennardJones, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);    

    temp = useElectrostatics;
    MPI_Allreduce(&temp, &useElectrostatics, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);    

    temp = useCharge;
    MPI_Allreduce(&temp, &useCharge, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);    

    temp = useDipole;
    MPI_Allreduce(&temp, &useDipole, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);    

    temp = useSticky;
    MPI_Allreduce(&temp, &useSticky, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);    

    temp = useGayBerne;
    MPI_Allreduce(&temp, &useGayBerne, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);    

    temp = useEAM;
    MPI_Allreduce(&temp, &useEAM, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);    

    temp = useShape;
    MPI_Allreduce(&temp, &useShape, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);   

    temp = useFLARB;
    MPI_Allreduce(&temp, &useFLARB, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);    

    temp = useRF;
    MPI_Allreduce(&temp, &useRF, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);    
    
#endif

    fInfo_.SIM_uses_PBC = usePBC;    
    fInfo_.SIM_uses_DirectionalAtoms = useDirectionalAtom;
    fInfo_.SIM_uses_LennardJones = useLennardJones;
    fInfo_.SIM_uses_Electrostatics = useElectrostatics;    
    fInfo_.SIM_uses_Charges = useCharge;
    fInfo_.SIM_uses_Dipoles = useDipole;
    fInfo_.SIM_uses_Sticky = useSticky;
    fInfo_.SIM_uses_GayBerne = useGayBerne;
    fInfo_.SIM_uses_EAM = useEAM;
    fInfo_.SIM_uses_Shapes = useShape;
    fInfo_.SIM_uses_FLARB = useFLARB;
    fInfo_.SIM_uses_RF = useRF;

    if( fInfo_.SIM_uses_Dipoles && fInfo_.SIM_uses_RF) {

        if (simParams_->haveDielectric()) {
            fInfo_.dielect = simParams_->getDielectric();
        } else {
            sprintf(painCave.errMsg,
                    "SimSetup Error: No Dielectric constant was set.\n"
                    "\tYou are trying to use Reaction Field without"
                    "\tsetting a dielectric constant!\n");
            painCave.isFatal = 1;
            simError();
        }
        
    } else {
        fInfo_.dielect = 0.0;
    }

}

void SimInfo::setupFortranSim() {
    int isError;
    int nExclude;
    std::vector<int> fortranGlobalGroupMembership;
    
    nExclude = exclude_.getSize();
    isError = 0;

    //globalGroupMembership_ is filled by SimCreator    
    for (int i = 0; i < nGlobalAtoms_; i++) {
        fortranGlobalGroupMembership.push_back(globalGroupMembership_[i] + 1);
    }

    //calculate mass ratio of cutoff group
    std::vector<double> mfact;
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::CutoffGroupIterator ci;
    CutoffGroup* cg;
    Molecule::AtomIterator ai;
    Atom* atom;
    double totalMass;

    //to avoid memory reallocation, reserve enough space for mfact
    mfact.reserve(getNCutoffGroups());
    
    for(mol = beginMolecule(mi); mol != NULL; mol = nextMolecule(mi)) {        
        for (cg = mol->beginCutoffGroup(ci); cg != NULL; cg = mol->nextCutoffGroup(ci)) {

            totalMass = cg->getMass();
            for(atom = cg->beginAtom(ai); atom != NULL; atom = cg->nextAtom(ai)) {
                        mfact.push_back(atom->getMass()/totalMass);
            }

        }       
    }

    //fill ident array of local atoms (it is actually ident of AtomType, it is so confusing !!!)
    std::vector<int> identArray;

    //to avoid memory reallocation, reserve enough space identArray
    identArray.reserve(getNAtoms());
    
    for(mol = beginMolecule(mi); mol != NULL; mol = nextMolecule(mi)) {        
        for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
            identArray.push_back(atom->getIdent());
        }
    }    

    //fill molMembershipArray
    //molMembershipArray is filled by SimCreator    
    std::vector<int> molMembershipArray(nGlobalAtoms_);
    for (int i = 0; i < nGlobalAtoms_; i++) {
        molMembershipArray[i] = globalMolMembership_[i] + 1;
    }
    
    //setup fortran simulation
    //gloalExcludes and molMembershipArray should go away (They are never used)
    //why the hell fortran need to know molecule?
    //OOPSE = Object-Obfuscated Parallel Simulation Engine
    int nGlobalExcludes = 0;
    int* globalExcludes = NULL; 
    int* excludeList = exclude_.getExcludeList();
    setFortranSim( &fInfo_, &nGlobalAtoms_, &nAtoms_, &identArray[0], &nExclude, excludeList , 
                  &nGlobalExcludes, globalExcludes, &molMembershipArray[0], 
                  &mfact[0], &nCutoffGroups_, &fortranGlobalGroupMembership[0], &isError); 

    if( isError ){

        sprintf( painCave.errMsg,
                 "There was an error setting the simulation information in fortran.\n" );
        painCave.isFatal = 1;
        painCave.severity = OOPSE_ERROR;
        simError();
    }

#ifdef IS_MPI
    sprintf( checkPointMsg,
       "succesfully sent the simulation information to fortran.\n");
    MPIcheckPoint();
#endif // is_mpi
}


#ifdef IS_MPI
void SimInfo::setupFortranParallel() {
    
    //SimInfo is responsible for creating localToGlobalAtomIndex and localToGlobalGroupIndex
    std::vector<int> localToGlobalAtomIndex(getNAtoms(), 0);
    std::vector<int> localToGlobalCutoffGroupIndex;
    SimInfo::MoleculeIterator mi;
    Molecule::AtomIterator ai;
    Molecule::CutoffGroupIterator ci;
    Molecule* mol;
    Atom* atom;
    CutoffGroup* cg;
    mpiSimData parallelData;
    int isError;

    for (mol = beginMolecule(mi); mol != NULL; mol  = nextMolecule(mi)) {

        //local index(index in DataStorge) of atom is important
        for (atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
            localToGlobalAtomIndex[atom->getLocalIndex()] = atom->getGlobalIndex() + 1;
        }

        //local index of cutoff group is trivial, it only depends on the order of travesing
        for (cg = mol->beginCutoffGroup(ci); cg != NULL; cg = mol->nextCutoffGroup(ci)) {
            localToGlobalCutoffGroupIndex.push_back(cg->getGlobalIndex() + 1);
        }        
        
    }

    //fill up mpiSimData struct
    parallelData.nMolGlobal = getNGlobalMolecules();
    parallelData.nMolLocal = getNMolecules();
    parallelData.nAtomsGlobal = getNGlobalAtoms();
    parallelData.nAtomsLocal = getNAtoms();
    parallelData.nGroupsGlobal = getNGlobalCutoffGroups();
    parallelData.nGroupsLocal = getNCutoffGroups();
    parallelData.myNode = worldRank;
    MPI_Comm_size(MPI_COMM_WORLD, &(parallelData.nProcessors));

    //pass mpiSimData struct and index arrays to fortran
    setFsimParallel(&parallelData, &(parallelData.nAtomsLocal),
                    &localToGlobalAtomIndex[0],  &(parallelData.nGroupsLocal),
                    &localToGlobalCutoffGroupIndex[0], &isError);

    if (isError) {
        sprintf(painCave.errMsg,
                "mpiRefresh errror: fortran didn't like something we gave it.\n");
        painCave.isFatal = 1;
        simError();
    }

    sprintf(checkPointMsg, " mpiRefresh successful.\n");
    MPIcheckPoint();


}

#endif

double SimInfo::calcMaxCutoffRadius() {


    std::set<AtomType*> atomTypes;
    std::set<AtomType*>::iterator i;
    std::vector<double> cutoffRadius;

    //get the unique atom types
    atomTypes = getUniqueAtomTypes();

    //query the max cutoff radius among these atom types
    for (i = atomTypes.begin(); i != atomTypes.end(); ++i) {
        cutoffRadius.push_back(forceField_->getRcutFromAtomType(*i));
    }

    double maxCutoffRadius = *(std::max_element(cutoffRadius.begin(), cutoffRadius.end()));
#ifdef IS_MPI
    //pick the max cutoff radius among the processors
#endif

    return maxCutoffRadius;
}

void SimInfo::getCutoff(double& rcut, double& rsw) {
    
    if (fInfo_.SIM_uses_Charges | fInfo_.SIM_uses_Dipoles | fInfo_.SIM_uses_RF) {
        
        if (!simParams_->haveRcut()){
            sprintf(painCave.errMsg,
                "SimCreator Warning: No value was set for the cutoffRadius.\n"
                "\tOOPSE will use a default value of 15.0 angstroms"
                "\tfor the cutoffRadius.\n");
            painCave.isFatal = 0;
            simError();
            rcut_ = 15.0;
        } else{
            rcut_ = simParams_->getRcut();
        }

        if (!simParams_->haveRsw()){
            sprintf(painCave.errMsg,
                "SimCreator Warning: No value was set for switchingRadius.\n"
                "\tOOPSE will use a default value of\n"
                "\t0.95 * cutoffRadius for the switchingRadius\n");
            painCave.isFatal = 0;
            simError();
            rsw_ = 0.95 * rcut_;
        } else{
            rsw_ = simParams_->getRsw();
        }

    } else {
        // if charge, dipole or reaction field is not used and the cutofff radius is not specified in
        //meta-data file, the maximum cutoff radius calculated from forcefiled will be used
        
        if (simParams_->haveRcut()) {
            rcut_ = simParams_->getRcut();
        } else {
            //set cutoff radius to the maximum cutoff radius based on atom types in the whole system
            rcut_ = calcMaxCutoffRadius();
        }

        if (simParams_->haveRsw()) {
            rsw_  = simParams_->getRsw();
        } else {
            rsw_ = rcut_;
        }
    
    }
}

void SimInfo::setupCutoff() {
    getCutoff(rcut_, rsw_);    
    double rnblist = rcut_ + 1; // skin of neighbor list

    //Pass these cutoff radius etc. to fortran. This function should be called once and only once
    notifyFortranCutoffs(&rcut_, &rsw_, &rnblist);
}

void SimInfo::addProperty(GenericData* genData) {
    properties_.addProperty(genData);  
}

void SimInfo::removeProperty(const std::string& propName) {
    properties_.removeProperty(propName);  
}

void SimInfo::clearProperties() {
    properties_.clearProperties(); 
}

std::vector<std::string> SimInfo::getPropertyNames() {
    return properties_.getPropertyNames();  
}
      
std::vector<GenericData*> SimInfo::getProperties() { 
    return properties_.getProperties(); 
}

GenericData* SimInfo::getPropertyByName(const std::string& propName) {
    return properties_.getPropertyByName(propName); 
}

void SimInfo::setSnapshotManager(SnapshotManager* sman) {
    sman_ = sman;

    Molecule* mol;
    RigidBody* rb;
    Atom* atom;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    Molecule::AtomIterator atomIter;;
 
    for (mol = beginMolecule(mi); mol != NULL; mol = nextMolecule(mi)) {
        
        for (atom = mol->beginAtom(atomIter); atom != NULL; atom = mol->nextAtom(atomIter)) {
            atom->setSnapshotManager(sman_);
        }
        
        for (rb = mol->beginRigidBody(rbIter); rb != NULL; rb = mol->nextRigidBody(rbIter)) {
            rb->setSnapshotManager(sman_);
        }
    }    
    
}

Vector3d SimInfo::getComVel(){ 
    SimInfo::MoleculeIterator i;
    Molecule* mol;

    Vector3d comVel(0.0);
    double totalMass = 0.0;
    
 
    for (mol = beginMolecule(i); mol != NULL; mol = nextMolecule(i)) {
        double mass = mol->getMass();
        totalMass += mass;
        comVel += mass * mol->getComVel();
    }  

#ifdef IS_MPI
    double tmpMass = totalMass;
    Vector3d tmpComVel(comVel);    
    MPI_Allreduce(&tmpMass,&totalMass,1,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(tmpComVel.getArrayPointer(), comVel.getArrayPointer(),3,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
#endif

    comVel /= totalMass;

    return comVel;
}

Vector3d SimInfo::getCom(){ 
    SimInfo::MoleculeIterator i;
    Molecule* mol;

    Vector3d com(0.0);
    double totalMass = 0.0;
     
    for (mol = beginMolecule(i); mol != NULL; mol = nextMolecule(i)) {
        double mass = mol->getMass();
        totalMass += mass;
        com += mass * mol->getCom();
    }  

#ifdef IS_MPI
    double tmpMass = totalMass;
    Vector3d tmpCom(com);    
    MPI_Allreduce(&tmpMass,&totalMass,1,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(tmpCom.getArrayPointer(), com.getArrayPointer(),3,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
#endif

    com /= totalMass;

    return com;

}        

std::ostream& operator <<(std::ostream& o, SimInfo& info) {

    return o;
}

}//end namespace oopse

