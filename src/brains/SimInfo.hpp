/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
/**
 * @file SimInfo.hpp
 * @author    tlin
 * @date  11/02/2004
 * @version 1.0
 */ 

#ifndef BRAINS_SIMMODEL_HPP
#define BRAINS_SIMMODEL_HPP

#include <iostream>
#include <set>
#include <utility>
#include <vector>

#include "brains/PairList.hpp"
#include "io/Globals.hpp"
#include "math/Vector3.hpp"
#include "math/SquareMatrix3.hpp"
#include "types/MoleculeStamp.hpp"
#include "brains/ForceField.hpp"
#include "utils/PropertyMap.hpp"
#include "utils/LocalIndexManager.hpp"
#include "nonbonded/SwitchingFunction.hpp"

using namespace std;
namespace OpenMD{
  //forward declaration 
  class SnapshotManager;
  class Molecule;
  class SelectionManager;
  class StuntDouble;

  /**
   * @class SimInfo SimInfo.hpp "brains/SimInfo.hpp"
   *
   * @brief One of the heavy-weight classes of OpenMD, SimInfo
   * maintains objects and variables relating to the current
   * simulation.  This includes the master list of Molecules.  The
   * Molecule class maintains all of the concrete objects (Atoms,
   * Bond, Bend, Torsions, Inversions, RigidBodies, CutoffGroups,
   * Constraints). In both the single and parallel versions, Atoms and
   * RigidBodies have both global and local indices.
   */
  class SimInfo {
  public:
    typedef map<int, Molecule*>::iterator  MoleculeIterator;
    
    /**
     * Constructor of SimInfo
     *
     * @param molStampPairs MoleculeStamp Array. The first element of
     * the pair is molecule stamp, the second element is the total
     * number of molecules with the same molecule stamp in the system
     *
     * @param ff pointer of a concrete ForceField instance
     *
     * @param simParams 
     */
    SimInfo(ForceField* ff, Globals* simParams);
    virtual ~SimInfo();

    /**
     * Adds a molecule
     *
     * @return return true if adding successfully, return false if the
     * molecule is already in SimInfo
     *
     * @param mol molecule to be added
     */
    bool addMolecule(Molecule* mol);

    /**
     * Removes a molecule from SimInfo
     *
     * @return true if removing successfully, return false if molecule
     * is not in this SimInfo
     */
    bool removeMolecule(Molecule* mol);

    /** Returns the total number of molecules in the system. */
    int getNGlobalMolecules() {
      return nGlobalMols_;
    }

    /** Returns the total number of atoms in the system. */
    int getNGlobalAtoms() {
      return nGlobalAtoms_;
    }

    /** Returns the total number of cutoff groups in the system. */
    int getNGlobalCutoffGroups() {
      return nGlobalCutoffGroups_;
    }

    /**
     * Returns the total number of integrable objects (total number of
     * rigid bodies plus the total number of atoms which do not belong
     * to the rigid bodies) in the system
     */
    int getNGlobalIntegrableObjects() {
      return nGlobalIntegrableObjects_;
    }

    /**
     * Returns the total number of integrable objects (total number of
     * rigid bodies plus the total number of atoms which do not belong
     * to the rigid bodies) in the system
     */
    int getNGlobalRigidBodies() {
      return nGlobalRigidBodies_;
    }

    int getNGlobalConstraints();
    /** 
     * Returns the number of local molecules.
     * @return the number of local molecules 
     */
    int getNMolecules() {
      return molecules_.size();
    }

    /** Returns the number of local atoms */
    unsigned int getNAtoms() {
      return nAtoms_;
    }

    /** Returns the number of effective cutoff groups on local processor */
    unsigned int getNLocalCutoffGroups();

    /** Returns the number of local bonds */        
    unsigned int getNBonds(){
      return nBonds_;
    }

    /** Returns the number of local bends */        
    unsigned int getNBends() {
      return nBends_;
    }

    /** Returns the number of local torsions */        
    unsigned int getNTorsions() {
      return nTorsions_;
    }

    /** Returns the number of local torsions */        
    unsigned int getNInversions() {
      return nInversions_;
    }
    /** Returns the number of local rigid bodies */        
    unsigned int getNRigidBodies() {
      return nRigidBodies_;
    }

    /** Returns the number of local integrable objects */
    unsigned int getNIntegrableObjects() {
      return nIntegrableObjects_;
    }

    /** Returns the number of local cutoff groups */
    unsigned int getNCutoffGroups() {
      return nCutoffGroups_;
    }

    /** Returns the total number of constraints in this SimInfo */
    unsigned int getNConstraints() {
      return nConstraints_;
    }
        
    /**
     * Returns the first molecule in this SimInfo and intialize the iterator.
     * @return the first molecule, return NULL if there is not molecule in this SimInfo
     * @param i the iterator of molecule array (user shouldn't change it)
     */
    Molecule* beginMolecule(MoleculeIterator& i);

    /** 
     * Returns the next avaliable Molecule based on the iterator.
     * @return the next avaliable molecule, return NULL if reaching the end of the array 
     * @param i the iterator of molecule array
     */
    Molecule* nextMolecule(MoleculeIterator& i);

    /** Returns the total number of fluctuating charges that are present */
    int getNFluctuatingCharges() {
      return nGlobalFluctuatingCharges_;
    }

    /** Returns the number of degrees of freedom */
    int getNdf() {
      return ndf_ - getFdf();
    }

    /** Returns the number of raw degrees of freedom */
    int getNdfRaw() {
      return ndfRaw_;
    }

    /** Returns the number of translational degrees of freedom */
    int getNdfTrans() {
      return ndfTrans_;
    }

    /** sets the current number of frozen degrees of freedom */
    void setFdf(int fdf) {
      fdf_local = fdf;
    }

    int getFdf(); 
    
    //getNZconstraint and setNZconstraint ruin the coherence of
    //SimInfo class, need refactoring
        
    /** Returns the total number of z-constraint molecules in the system */
    int getNZconstraint() {
      return nZconstraint_;
    }

    /** 
     * Sets the number of z-constraint molecules in the system.
     */
    void setNZconstraint(int nZconstraint) {
      nZconstraint_ = nZconstraint;
    }
        
    /** Returns the snapshot manager. */
    SnapshotManager* getSnapshotManager() {
      return sman_;
    }

    /** Sets the snapshot manager. */
    void setSnapshotManager(SnapshotManager* sman);
        
    /** Returns the force field */
    ForceField* getForceField() {
      return forceField_;
    }

    Globals* getSimParams() {
      return simParams_;
    }

    /** Returns the velocity of center of mass of the whole system.*/
    Vector3d getComVel();

    /** Returns the center of the mass of the whole system.*/
    Vector3d getCom();
    /** Returns the center of the mass and Center of Mass velocity of
        the whole system.*/ 
    void getComAll(Vector3d& com,Vector3d& comVel);

    /** Returns intertia tensor for the entire system and system
        Angular Momentum.*/
    void getInertiaTensor(Mat3x3d &intertiaTensor,Vector3d &angularMomentum);
    
    /** Returns system angular momentum */
    Vector3d getAngularMomentum();

    /** Returns volume of system as estimated by an ellipsoid defined
        by the radii of gyration*/
    void getGyrationalVolume(RealType &vol);
    /** Overloaded version of gyrational volume that also returns
        det(I) so dV/dr can be calculated*/
    void getGyrationalVolume(RealType &vol, RealType &detI);

    void update();
    /**
     * Do final bookkeeping before Force managers need their data.
     */
    void prepareTopology();


    /** Returns the local index manager */
    LocalIndexManager* getLocalIndexManager() {
      return &localIndexMan_;
    }

    int getMoleculeStampId(int globalIndex) {
      //assert(globalIndex < molStampIds_.size())
      return molStampIds_[globalIndex];
    }

    /** Returns the molecule stamp */
    MoleculeStamp* getMoleculeStamp(int id) {
      return moleculeStamps_[id];
    }

    /** Return the total number of the molecule stamps */
    int getNMoleculeStamp() {
      return moleculeStamps_.size();
    }
    /**
     * Finds a molecule with a specified global index
     * @return a pointer point to found molecule
     * @param index
     */
    Molecule* getMoleculeByGlobalIndex(int index) {
      MoleculeIterator i;
      i = molecules_.find(index);

      return i != molecules_.end() ? i->second : NULL;
    }

    int getGlobalMolMembership(int id){
      return globalMolMembership_[id];
    }

    /** 
     * returns a vector which maps the local atom index on this
     * processor to the global atom index.  With only one processor,
     * these should be identical.
     */
    vector<int> getGlobalAtomIndices();

    /** 
     * returns a vector which maps the local cutoff group index on
     * this processor to the global cutoff group index.  With only one
     * processor, these should be identical.
     */
    vector<int> getGlobalGroupIndices();

        
    string getFinalConfigFileName() {
      return finalConfigFileName_;
    }

    void setFinalConfigFileName(const string& fileName) {
      finalConfigFileName_ = fileName;
    }

    string getRawMetaData() {
      return rawMetaData_;
    }
    void setRawMetaData(const string& rawMetaData) {
      rawMetaData_ = rawMetaData;
    }
        
    string getDumpFileName() {
      return dumpFileName_;
    }
        
    void setDumpFileName(const string& fileName) {
      dumpFileName_ = fileName;
    }

    string getStatFileName() {
      return statFileName_;
    }
        
    void setStatFileName(const string& fileName) {
      statFileName_ = fileName;
    }
        
    string getRestFileName() {
      return restFileName_;
    }
        
    void setRestFileName(const string& fileName) {
      restFileName_ = fileName;
    }

    /** 
     * Sets GlobalGroupMembership
     * @see #SimCreator::setGlobalIndex
     */  
    void setGlobalGroupMembership(const vector<int>& globalGroupMembership) {
      assert(globalGroupMembership.size() == static_cast<size_t>(nGlobalAtoms_));
      globalGroupMembership_ = globalGroupMembership;
    }

    /** 
     * Sets GlobalMolMembership
     * @see #SimCreator::setGlobalIndex
     */        
    void setGlobalMolMembership(const vector<int>& globalMolMembership) {
      assert(globalMolMembership.size() == static_cast<size_t>(nGlobalAtoms_));
      globalMolMembership_ = globalMolMembership;
    }


    bool isTopologyDone() {
      return topologyDone_;
    }
        
    bool getCalcBoxDipole() {
      return calcBoxDipole_;
    }

    bool getUseAtomicVirial() {
      return useAtomicVirial_;
    }

    /**
     * Adds property into property map
     * @param genData GenericData to be added into PropertyMap
     */
    void addProperty(GenericData* genData);

    /**
     * Removes property from PropertyMap by name
     * @param propName the name of property to be removed
     */
    void removeProperty(const string& propName);

    /**
     * clear all of the properties
     */
    void clearProperties();

    /**
     * Returns all names of properties
     * @return all names of properties
     */
    vector<string> getPropertyNames();

    /**
     * Returns all of the properties in PropertyMap
     * @return all of the properties in PropertyMap
     */      
    vector<GenericData*> getProperties();

    /**
     * Returns property 
     * @param propName name of property
     * @return a pointer point to property with propName. If no property named propName
     * exists, return NULL
     */      
    GenericData* getPropertyByName(const string& propName);

    /**
     * add all special interaction pairs (including excluded
     * interactions) in a molecule into the appropriate lists.
     */
    void addInteractionPairs(Molecule* mol);

    /**
     * remove all special interaction pairs which belong to a molecule
     * from the appropriate lists.
     */
    void removeInteractionPairs(Molecule* mol);

    /** Returns the set of atom types present in this simulation */
    set<AtomType*> getSimulatedAtomTypes();
        
    friend ostream& operator <<(ostream& o, SimInfo& info);

    void getCutoff(RealType& rcut, RealType& rsw);
        
  private:

    /** fill up the simtype struct and other simulation-related variables */
    void setupSimVariables();


    /** Determine if we need to accumulate the simulation box dipole */
    void setupAccumulateBoxDipole();

    /** Calculates the number of degress of freedom in the whole system */
    void calcNdf();
    void calcNdfRaw();
    void calcNdfTrans();

    /**
     * Adds molecule stamp and the total number of the molecule with
     * same molecule stamp in the whole system.
     */
    void addMoleculeStamp(MoleculeStamp* molStamp, int nmol);

    // Other classes holdingn important information
    ForceField* forceField_; /**< provides access to defined atom types, bond types, etc. */
    Globals* simParams_;     /**< provides access to simulation parameters set by user */

    ///  Counts of local objects
    int nAtoms_;              /**< number of atoms in local processor */
    int nBonds_;              /**< number of bonds in local processor */
    int nBends_;              /**< number of bends in local processor */
    int nTorsions_;           /**< number of torsions in local processor */
    int nInversions_;         /**< number of inversions in local processor */
    int nRigidBodies_;        /**< number of rigid bodies in local processor */
    int nIntegrableObjects_;  /**< number of integrable objects in local processor */
    int nCutoffGroups_;       /**< number of cutoff groups in local processor */
    int nConstraints_;        /**< number of constraints in local processors */
    int nFluctuatingCharges_; /**< number of fluctuating charges in local processor */
        
    /// Counts of global objects
    int nGlobalMols_;              /**< number of molecules in the system (GLOBAL) */
    int nGlobalAtoms_;             /**< number of atoms in the system (GLOBAL) */
    int nGlobalCutoffGroups_;      /**< number of cutoff groups in this system (GLOBAL) */
    int nGlobalIntegrableObjects_; /**< number of integrable objects in this system */
    int nGlobalRigidBodies_;       /**< number of rigid bodies in this system (GLOBAL) */
    int nGlobalFluctuatingCharges_;/**< number of fluctuating charges in this system (GLOBAL) */
    
       
    /// Degress of freedom
    int ndf_;          /**< number of degress of freedom (excludes constraints) (LOCAL) */
    int fdf_local;     /**< number of frozen degrees of freedom (LOCAL) */
    int fdf_;          /**< number of frozen degrees of freedom (GLOBAL) */
    int ndfRaw_;       /**< number of degress of freedom (includes constraints),  (LOCAL) */
    int ndfTrans_;     /**< number of translation degress of freedom, (LOCAL) */
    int nZconstraint_; /**< number of  z-constraint molecules (GLOBAL) */

    /// logicals 
    bool usesPeriodicBoundaries_; /**< use periodic boundary conditions? */
    bool usesDirectionalAtoms_;   /**< are there atoms with position AND orientation? */
    bool usesMetallicAtoms_;      /**< are there transition metal atoms? */
    bool usesElectrostaticAtoms_; /**< are there electrostatic atoms? */
    bool usesFluctuatingCharges_; /**< are there fluctuating charges? */
    bool usesAtomicVirial_;       /**< are we computing atomic virials? */
    bool requiresPrepair_;        /**< does this simulation require a pre-pair loop? */
    bool requiresSkipCorrection_; /**< does this simulation require a skip-correction? */
    bool requiresSelfCorrection_; /**< does this simulation require a self-correction? */

  public:
    bool usesElectrostaticAtoms() { return usesElectrostaticAtoms_; }
    bool usesDirectionalAtoms() { return usesDirectionalAtoms_; }
    bool usesFluctuatingCharges() { return usesFluctuatingCharges_; }
    bool usesAtomicVirial() { return usesAtomicVirial_; }
    bool requiresPrepair() { return requiresPrepair_; }
    bool requiresSkipCorrection() { return requiresSkipCorrection_;}
    bool requiresSelfCorrection() { return requiresSelfCorrection_;}

  private:
    /// Data structures holding primary simulation objects
    map<int, Molecule*>  molecules_;  /**< map holding pointers to LOCAL molecules */

    /// Stamps are templates for objects that are then used to create
    /// groups of objects.  For example, a molecule stamp contains
    /// information on how to build that molecule (i.e. the topology,
    /// the atoms, the bonds, etc.)  Once the system is built, the
    /// stamps are no longer useful.
    vector<int> molStampIds_;                /**< stamp id for molecules in the system */
    vector<MoleculeStamp*> moleculeStamps_;  /**< molecule stamps array */        

    /**
     * A vector that maps between the global index of an atom, and the
     * global index of cutoff group the atom belong to.  It is filled
     * by SimCreator once and only once, since it never changed during
     * the simulation.  It should be nGlobalAtoms_ in size.
     */
    vector<int> globalGroupMembership_; 
  public:
    vector<int> getGlobalGroupMembership() { return globalGroupMembership_; }
  private:

    /**
     * A vector that maps between the global index of an atom and the
     * global index of the molecule the atom belongs to.  It is filled
     * by SimCreator once and only once, since it is never changed
     * during the simulation. It shoudl be nGlobalAtoms_ in size.
     */
    vector<int> globalMolMembership_; 

    /** 
     * A vector that maps between the local index of an atom and the
     * index of the AtomType.
     */
    vector<int> identArray_;
  public:
    vector<int> getIdentArray() { return identArray_; }
  private:
    
    /** 
     * A vector which contains the fractional contribution of an
     * atom's mass to the total mass of the cutoffGroup that atom
     * belongs to.  In the case of single atom cutoff groups, the mass
     * factor for that atom is 1.  For massless atoms, the factor is
     * also 1.
     */
    vector<RealType> massFactors_;
  public:
    vector<RealType> getMassFactors() { return massFactors_; }

    PairList* getExcludedInteractions() { return &excludedInteractions_; }
    PairList* getOneTwoInteractions() { return &oneTwoInteractions_; }
    PairList* getOneThreeInteractions() { return &oneThreeInteractions_; }
    PairList* getOneFourInteractions() { return &oneFourInteractions_; }

  private:
               
    /// lists to handle atoms needing special treatment in the non-bonded interactions
    PairList excludedInteractions_;  /**< atoms excluded from interacting with each other */
    PairList oneTwoInteractions_;    /**< atoms that are directly Bonded */ 
    PairList oneThreeInteractions_;  /**< atoms sharing a Bend */    
    PairList oneFourInteractions_;   /**< atoms sharing a Torsion */

    PropertyMap properties_;       /**< Generic Properties can be added */
    SnapshotManager* sman_;        /**< SnapshotManager (handles particle positions, etc.) */

    /** 
     * The reason to have a local index manager is that when molecule
     * is migrating to other processors, the atoms and the
     * rigid-bodies will release their local indices to
     * LocalIndexManager. Combining the information of molecule
     * migrating to current processor, Migrator class can query the
     * LocalIndexManager to make a efficient data moving plan.
     */        
    LocalIndexManager localIndexMan_;

    // unparsed MetaData block for storing in Dump and EOR files:
    string rawMetaData_;

    // file names
    string finalConfigFileName_;
    string dumpFileName_;
    string statFileName_;
    string restFileName_;
        

    bool topologyDone_;  /** flag to indicate whether the topology has
                             been scanned and all the relevant
                             bookkeeping has been done*/
    
    bool calcBoxDipole_; /**< flag to indicate whether or not we calculate 
                            the simulation box dipole moment */
    
    bool useAtomicVirial_; /**< flag to indicate whether or not we use 
                              Atomic Virials to calculate the pressure */
    
  public:
    /**
     * return an integral objects by its global index. In MPI
     * version, if the StuntDouble with specified global index does
      * not belong to local processor, a NULL will be return.
      */
    StuntDouble* getIOIndexToIntegrableObject(int index);
    void setIOIndexToIntegrableObject(const vector<StuntDouble*>& v);
    
  private:
    vector<StuntDouble*> IOIndexToIntegrableObject;
    
  public:
                
    /**
     * Finds the processor where a molecule resides
     * @return the id of the processor which contains the molecule
     * @param globalIndex global Index of the molecule
     */
    int getMolToProc(int globalIndex) {
      //assert(globalIndex < molToProcMap_.size());
      return molToProcMap_[globalIndex];
    }
    
    /** 
     * Set MolToProcMap array
     * @see #SimCreator::divideMolecules
     */
    void setMolToProcMap(const vector<int>& molToProcMap) {
      molToProcMap_ = molToProcMap;
    }
        
  private:
        
    /** 
     * The size of molToProcMap_ is equal to total number of molecules
     * in the system.  It maps a molecule to the processor on which it
     * resides. it is filled by SimCreator once and only once.
     */        
    vector<int> molToProcMap_; 

  };

} //namespace OpenMD
#endif //BRAINS_SIMMODEL_HPP

