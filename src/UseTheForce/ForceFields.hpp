#ifndef __FORCEFIELDS_H__
#define __FORCEFIELDS_H__

#define MK_STR(s) # s
#define STR_DEFINE(t, s) t = MK_STR(s)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <map>

#include "primitives/Atom.hpp"
#include "brains/SimInfo.hpp"
#include "primitives/StuntDouble.hpp"
#include "types/ShapeAtomType.hpp"
#include "io/basic_ifstrstream.hpp"

#ifdef IS_MPI
#include "UseTheForce/mpiForceField.h"
#endif

using namespace std;
using namespace oopse;

class bond_pair{
public:
  bond_pair(){}
  ~bond_pair(){}

  int a;
  int b;
};

class bend_set{
public:
  bend_set(){ isGhost = 0; }
  ~bend_set(){}

  int ghost;
  int isGhost;

  int a;
  int b;
  int c;
};

class torsion_set{
public:
  torsion_set(){}
  ~torsion_set(){}

  int a;
  int b;
  int c;
  int d;
};



class ForceFields{

public:
  ForceFields(){ 
    frcFile = NULL; 
    entry_plug = NULL; 
    has_variant=0;
    ffPath = getenv("FORCE_PARAM_PATH");
    if( ffPath.empty() ) {
      STR_DEFINE(ffPath, FRC_PATH );
    }   
  }
  ForceFields(string theVariant){ 
    frcFile = NULL; 
    entry_plug = NULL; 
    has_variant=1; 
    variant = theVariant;
    ffPath = getenv("FORCE_PARAM_PATH");
    if( ffPath.empty() ) {
      STR_DEFINE(ffPath, FRC_PATH );
    }   
  }
  virtual ~ForceFields(){}
  
  void setSimInfo( SimInfo* the_entry_plug ) { entry_plug = the_entry_plug; }
  
  virtual void readParams( void ) = 0;
  
  virtual void cleanMe( void ) = 0;


  virtual void initializeAtoms( int nAtoms, Atom** atomArray ) = 0;
  virtual void initializeBonds( int nBonds, Bond** bondArray,
				bond_pair* the_bonds ) = 0;
  virtual void initializeBends( int nBends, Bend** bendArray,
				bend_set* the_bends ) = 0;
  virtual void initializeTorsions( int nTorsions, Torsion** torsionArray,
				   torsion_set* the_torsions ) = 0;
  virtual void initForceField() = 0;
  virtual void initRestraints();
  virtual void dumpzAngle();

  virtual void calcRcut( void );
  virtual void setRcut( double LJrcut );
  virtual void doForces( int calcPot, int calcStress );

protected:
  
  void initFortran( int useReactionField );
  
  FILE *frcFile;
  SimInfo* entry_plug;
  
  int lineNum;
  char readLine[500];
  char* eof_test;
  short int has_variant;
  double bigSigma;

  string ffPath;
  ifstrstream forceFile;
  bool hasVariant;
  string variant;
  map<string, AtomType*> atomTypeMap;
  // map<pair<string,string>, BondType*> bondTypeMap;
  // map<tuple3<string,string,string>, BendType*> bendTypeMap;
  // map<tuple4<string,string,string,string>, TorsionType*> torsionTypeMap;

};


class DUFF : public ForceFields{

public:
  DUFF();
  virtual ~DUFF();

  void readParams();
  void cleanMe( void );

  void initializeAtoms( int nAtoms, Atom** atomArray );
  void initializeBonds( int nBonds, Bond** bondArray,
			bond_pair* the_bonds );
  void initializeBends( int nBends, Bend** bendArray,
			bend_set* the_bends );
  void initializeTorsions( int nTorsions, Torsion** torsionArray,
			   torsion_set* the_torsions );

  void initForceField();

private:
  
  void fastForward( char* stopText, char* searchOwner );
};

class LJFF : public ForceFields{

public:
  LJFF();
  virtual ~LJFF();
  

  void readParams();
  void cleanMe( void );

  void initializeAtoms( int nAtoms, Atom** atomArray );
  void initializeBonds( int nBonds, Bond** bondArray,
			bond_pair* the_bonds );
  void initializeBends( int nBends, Bend** bendArray,
			bend_set* the_bends );
  void initializeTorsions( int nTorsions, Torsion** torsionArray,
			   torsion_set* the_torsions );

  void initForceField( );


private:

  void fastForward( char* stopText, char* searchOwner );

};

class EAM_FF : public ForceFields{

public:
  EAM_FF();
  EAM_FF(string the_variant);
  virtual ~EAM_FF();
  

  void readParams();
  void cleanMe( void );

  void initializeAtoms( int nAtoms, Atom** atomArray );
  void initializeBonds( int nBonds, Bond** bondArray,
			bond_pair* the_bonds );
  void initializeBends( int nBends, Bend** bendArray,
			bend_set* the_bends );
  void initializeTorsions( int nTorsions, Torsion** torsionArray,
			   torsion_set* the_torsions );

  void initForceField();

  void calcRcut( void );

private:

  void fastForward( char* stopText, char* searchOwner );
  
  double eamRcut;
};

class WATER : public ForceFields{

public:
  WATER();
  virtual ~WATER();

  void readParams();
  void cleanMe( void );
  void initializeAtoms( int nAtoms, Atom** atomArray );
  void initializeBonds( int nBonds, Bond** bondArray,
			bond_pair* the_bonds );
  void initializeBends( int nBends, Bend** bendArray,
			bend_set* the_bends );
  void initializeTorsions( int nTorsions, Torsion** torsionArray,
			   torsion_set* the_torsions );
  void initForceField();


private:
  
  void fastForward( char* stopText, char* searchOwner );
  void sectionSearch( char* secHead, char* stopText, char* searchOwner );

};

class Shapes_FF : public ForceFields{

public:
  Shapes_FF() : ForceFields() {};
  Shapes_FF(char* the_variant) : ForceFields(the_variant) {};
  virtual ~Shapes_FF();
  
  void readParams();
  void cleanMe( void );

  void initializeAtoms( int nAtoms, Atom** atomArray );
  void initializeBonds( int nBonds, Bond** bondArray,
			bond_pair* the_bonds );
  void initializeBends( int nBends, Bend** bendArray,
			bend_set* the_bends );
  void initializeTorsions( int nTorsions, Torsion** torsionArray,
			   torsion_set* the_torsions );

  void initForceField();

  void parseShapeFile(string shapeFileName, ShapeAtomType* st);

  void calcRcut( void );

private:

  double findLargestContactDistance(ShapeAtomType* st);
  double shapesRcut;

};


#endif

