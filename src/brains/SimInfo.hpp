#ifndef __SIMINFO_H__
#define __SIMINFO_H__

#include <map>
#include <string>
#include <vector>

#include "Atom.hpp"
#include "RigidBody.hpp"
#include "Molecule.hpp"
#include "Exclude.hpp"
#include "SkipList.hpp"
#include "AbstractClasses.hpp"
#include "MakeStamps.hpp"
#include "SimState.hpp"
#include "Restraints.hpp"

#define __C
#include "fSimulation.h"
#include "fortranWrapDefines.hpp"
#include "GenericData.hpp"


//#include "Minimizer.hpp"
//#include "OOPSEMinimizer.hpp"


double roundMe( double x );
class OOPSEMinimizer;
class SimInfo{

public:

  SimInfo();
  ~SimInfo();

  int n_atoms; // the number of atoms
  Atom **atoms; // the array of atom objects

  vector<RigidBody*> rigidBodies;  // A vector of rigid bodies
  vector<StuntDouble*> integrableObjects;
  
  double tau[9]; // the stress tensor

  int n_bonds;    // number of bends
  int n_bends;    // number of bends
  int n_torsions; // number of torsions
  int n_oriented; // number of of atoms with orientation
  int ndf;        // number of actual degrees of freedom
  int ndfRaw;     // number of settable degrees of freedom
  int ndfTrans;   // number of translational degrees of freedom
  int nZconstraints; // the number of zConstraints

  int setTemp;   // boolean to set the temperature at each sampleTime
  int resetIntegrator; // boolean to reset the integrator

  int n_dipoles; // number of dipoles

  int n_exclude;
  Exclude* excludes;  // the exclude list for ignoring pairs in fortran
  int nGlobalExcludes;
  int* globalExcludes; // same as above, but these guys participate in
		       // no long range forces.

  int* identArray;     // array of unique identifiers for the atoms
  int* molMembershipArray;  // map of atom numbers onto molecule numbers

  int n_constraints; // the number of constraints on the system

  int n_SRI;   // the number of short range interactions

  double lrPot; // the potential energy from the long range calculations.

  double Hmat[3][3];  // the periodic boundry conditions. The Hmat is the
                      // column vectors of the x, y, and z box vectors.
                      //   h1  h2  h3
                      // [ Xx  Yx  Zx ]
                      // [ Xy  Yy  Zy ]
                      // [ Xz  Yz  Zz ]
                      //   
  double HmatInv[3][3];

  double boxL[3]; // The Lengths of the 3 column vectors of Hmat
  double boxVol;
  int orthoRhombic;
  

  double dielectric;      // the dielectric of the medium for reaction field

  
  int usePBC; // whether we use periodic boundry conditions.
  int useLJ; 
  int useSticky;
  int useCharges;
  int useDipoles;
  int useReactionField;
  int useGB;
  int useEAM;
  bool haveCutoffGroups;
  bool useInitXSstate;
  double orthoTolerance;

  double dt, run_time;           // the time step and total time
  double sampleTime, statusTime; // the position and energy dump frequencies
  double target_temp;            // the target temperature of the system
  double thermalTime;            // the temp kick interval
  double currentTime;            // Used primarily for correlation Functions
  double resetTime;              // Use to reset the integrator periodically
  short int have_target_temp;

  int n_mol;           // n_molecules;
  Molecule* molecules; // the array of molecules
  
  int nComponents;           // the number of components in the system
  int* componentsNmol;       // the number of molecules of each component
  MoleculeStamp** compStamps;// the stamps matching the components
  LinkedMolStamp* headStamp; // list of stamps used in the simulation
  
  
  char ensemble[100]; // the enesemble of the simulation (NVT, NVE, etc. )
  char mixingRule[100]; // the mixing rules for Lennard jones/van der walls 
  BaseIntegrator *the_integrator; // the integrator of the simulation

  OOPSEMinimizer* the_minimizer; // the energy minimizer
  Restraints* restraint;
  bool has_minimizer;

  string finalName;  // the name of the eor file to be written
  string sampleName; // the name of the dump file to be written
  string statusName; // the name of the stat file to be written

  int seed;                    //seed for random number generator

  int useSolidThermInt;  // is solid-state thermodynamic integration being used
  int useLiquidThermInt; // is liquid thermodynamic integration being used
  double thermIntLambda; // lambda for TI
  double thermIntK;      // power of lambda for TI
  double vRaw;           // unperturbed potential for TI
  double vHarm;          // harmonic potential for TI
  int i;                 // just an int

  vector<double> mfact;
  vector<int> FglobalGroupMembership;
  int ngroup;
  int* globalGroupMembership;

  // refreshes the sim if things get changed (load balanceing, volume
  // adjustment, etc.)

  void refreshSim( void );
  

  // sets the internal function pointer to fortran.

  void setInternal( setFortranSim_TD fSetup,
		    setFortranBox_TD fBox,
		    notifyFortranCutOff_TD fCut){
    setFsimulation = fSetup;
    setFortranBoxSize = fBox;
    notifyFortranCutOffs = fCut;
  }

  int getNDF();
  int getNDFraw();
  int getNDFtranslational();
  int getTotIntegrableObjects();
  void setBox( double newBox[3] );
  void setBoxM( double newBox[3][3] );
  void getBoxM( double theBox[3][3] );
  void scaleBox( double scale );
  
  void setDefaultRcut( double theRcut );
  void setDefaultRcut( double theRcut, double theRsw );
  void checkCutOffs( void );

  double getRcut( void )  { return rCut; }
  double getRlist( void ) { return rList; }
  double getRsw( void )   { return rSw; }
  double getMaxCutoff( void ) { return maxCutoff; }
  
  void setTime( double theTime ) { currentTime = theTime; }
  void incrTime( double the_dt ) { currentTime += the_dt; }
  void decrTime( double the_dt ) { currentTime -= the_dt; }
  double getTime( void ) { return currentTime; }

  void wrapVector( double thePos[3] );

  SimState* getConfiguration( void ) { return myConfiguration; }
  
  void addProperty(GenericData* prop); 
  GenericData* getProperty(const string& propName);
  //vector<GenericData*>& getProperties()  {return properties;}     

  int getSeed(void) {  return seed; }
  void setSeed(int theSeed) {  seed = theSeed;}

private:

  SimState* myConfiguration;

  int boxIsInit, haveRcut, haveRsw;

  double rList, rCut; // variables for the neighborlist
  double rSw;         // the switching radius

  double maxCutoff;

  double distXY;
  double distYZ;
  double distZX;
  
  void calcHmatInv( void );
  void calcBoxL();
  double calcMaxCutOff();

  // private function to initialize the fortran side of the simulation
  setFortranSim_TD setFsimulation;

  setFortranBox_TD setFortranBoxSize;
  
  notifyFortranCutOff_TD notifyFortranCutOffs;
  
  //Addtional Properties of SimInfo
  map<string, GenericData*> properties;
  void getFortranGroupArrays(SimInfo* info, 
                             vector<int>& FglobalGroupMembership,
                             vector<double>& mfact);


};


#endif
