#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
using namespace std;

#include "UseTheForce/ForceFields.hpp"
#include "primitives/SRI.hpp"
#include "utils/simError.h"
#include "UseTheForce/DarkSide/sticky_interface.h"
#include "UseTheForce/DarkSide/atype_interface.h"

//#include "UseTheForce/fortranWrappers.hpp"


#ifdef IS_MPI
#include "UseTheForce/mpiForceField.h"
#endif // is_mpi


// define some bond Types

#define FIXED_BOND    0
#define HARMONIC_BOND 1


namespace DUFF_NS {  // restrict the access of the folowing to this file only.


  // Declare the structures that will be passed by MPI
  
  typedef struct{
    char name[15];
    double mass;
    double epslon;
    double sigma;
    double charge;
    double dipole;
    double w0;
    double v0;
    double v0p;
    double rl;
    double ru;
    double rlp;
    double rup;
    int isSSD;
    int isCharge;
    int isDipole;
    int ident;
    int last;      //  0  -> default
                   //  1  -> tells nodes to stop listening
  } atomStruct;
  
  
  typedef struct{
    char nameA[15];
    char nameB[15];
    double d0;
    double k0;
    int last;      //  0  -> default
                   //  1  -> tells nodes to stop listening
    int type;
  } bondStruct;
  
  
  typedef struct{
    char nameA[15];
    char nameB[15];
    char nameC[15];
    char type[30];
    double k1, k2, k3, t0;
    int last;      //  0  -> default
                   //  1  -> tells nodes to stop listening
  } bendStruct;
  

  typedef struct{
    char nameA[15];
    char nameB[15];
    char nameC[15];
    char nameD[15];
    char type[30];
    double k1, k2, k3, k4;
    int last;      //  0  -> default
                   //  1  -> tells nodes to stop listening
  } torsionStruct;
  
  
  int parseAtom(    char *lineBuffer, int lineNum, atomStruct     &info );
  int parseBond(    char *lineBuffer, int lineNum, bondStruct     &info );
  int parseBend(    char *lineBuffer, int lineNum, bendStruct     &info );
  int parseTorsion( char *lineBuffer, int lineNum, torsionStruct  &info );
  
  
#ifdef IS_MPI

  MPI_Datatype mpiAtomStructType;
  MPI_Datatype mpiBondStructType;
  MPI_Datatype mpiBendStructType;
  MPI_Datatype mpiTorsionStructType;

#endif

  class LinkedAtomType {
  public:
    LinkedAtomType(){ 
      next = NULL;
      name[0] = '\0';
    }
    ~LinkedAtomType(){ if( next != NULL ) delete next; }

    LinkedAtomType* find(char* key){ 
      if( !strcmp(name, key) ) return this;
      if( next != NULL ) return next->find(key);
      return NULL;
    }
    
    void printMe( void ){
      
      std::cerr << "LinkedAtype " << name << ": ident = " << ident << "\n";
      //      if( next != NULL ) next->printMe();

    }

    void add( atomStruct &info ){

      // check for duplicates
      
      if( !strcmp( info.name, name ) ){
	sprintf( painCave.errMsg,
		 "Duplicate DUFF atom type \"%s\" found in "
		 "the DUFF param file./n",
		 name );
	painCave.isFatal = 1;
	simError();
      }

      if( next != NULL ) next->add(info);
      else{
	next = new LinkedAtomType();
	strcpy(next->name, info.name);
	next->isDipole = info.isDipole;
	next->isSSD    = info.isSSD;
	next->mass     = info.mass;
	next->epslon   = info.epslon;
	next->sigma    = info.sigma;
	next->dipole   = info.dipole;
	next->w0       = info.w0;
	next->v0       = info.v0;
	next->v0p      = info.v0p;
	next->rl       = info.rl;
	next->ru       = info.ru;
	next->rlp      = info.rlp;
	next->rup      = info.rup;
	next->ident    = info.ident;
      }
    }

#ifdef IS_MPI
    
    void duplicate( atomStruct &info ){
      strcpy(info.name, name);
      info.isDipole = isDipole;
      info.isSSD    = isSSD;
      info.mass     = mass;
      info.epslon   = epslon;
      info.sigma    = sigma;
      info.dipole   = dipole;
      info.w0       = w0;
      info.v0       = v0;
      info.v0p      = v0p;
      info.rl       = rl;
      info.ru       = ru;
      info.rlp      = rlp;
      info.rup      = rup;
      info.ident    = ident;
      info.last     = 0;
    }


#endif

    char name[15];
    int isDipole;
    int isSSD;
    double mass;
    double epslon;
    double sigma;
    double dipole;
    double w0;
    double v0;
    double v0p;
    double rl;
    double ru;
    double rlp;
    double rup;
    int ident;
    LinkedAtomType* next;
  };

  class LinkedBondType {
  public:
    LinkedBondType(){ 
      next = NULL;
      nameA[0] = '\0';
      nameB[0] = '\0';
    }
    ~LinkedBondType(){ if( next != NULL ) delete next; }

    LinkedBondType* find(char* key1, char* key2){ 
      if( !strcmp(nameA, key1 ) && !strcmp( nameB, key2 ) ) return this;
      if( !strcmp(nameA, key2 ) && !strcmp( nameB, key1 ) ) return this;
      if( next != NULL ) return next->find(key1, key2);
      return NULL;
    }
    

    void add( bondStruct &info ){
      
      // check for duplicates
      int dup = 0;

      if( !strcmp(nameA, info.nameA ) && !strcmp( nameB, info.nameB ) ) dup = 1;
      if( !strcmp(nameA, info.nameB ) && !strcmp( nameB, info.nameA ) ) dup = 1;
      
      if(dup){
	sprintf( painCave.errMsg,
		 "Duplicate DUFF bond type \"%s - %s\" found in "
		 "the DUFF param file./n",
		 nameA, nameB );
	painCave.isFatal = 1;
	simError();
      }

	
      if( next != NULL ) next->add(info);
      else{
	next = new LinkedBondType();
	strcpy(next->nameA, info.nameA);
	strcpy(next->nameB, info.nameB);
	next->type = info.type;
	next->d0 = info.d0;
	next->k0 = info.k0;
      }
    }
    
#ifdef IS_MPI
    void duplicate( bondStruct &info ){
      strcpy(info.nameA, nameA);
      strcpy(info.nameB, nameB);
      info.type = type;
      info.d0   = d0;
      info.k0   = k0;
      info.last = 0;
    }


#endif

    char nameA[15];
    char nameB[15];
    int type;
    double d0;
    double k0;

    LinkedBondType* next;
  };

  class LinkedBendType {
  public:
    LinkedBendType(){ 
      next = NULL;
      nameA[0] = '\0';
      nameB[0] = '\0';
      nameC[0] = '\0';
      type[0] = '\0';
    }
    ~LinkedBendType(){ if( next != NULL ) delete next; }

    LinkedBendType* find( char* key1, char* key2, char* key3 ){ 
      if( !strcmp( nameA, key1 ) && !strcmp( nameB, key2 ) 
	  && !strcmp( nameC, key3 ) ) return this;
      if( !strcmp( nameA, key3 ) && !strcmp( nameB, key2 ) 
	  && !strcmp( nameC, key1 ) ) return this;
      if( next != NULL ) return next->find(key1, key2, key3);
      return NULL;
    }
    
    void add( bendStruct &info ){

      // check for duplicates
      int dup = 0;
      
      if( !strcmp( nameA, info.nameA ) && !strcmp( nameB, info.nameB ) 
	  && !strcmp( nameC, info.nameC ) ) dup = 1;
      if( !strcmp( nameA, info.nameC ) && !strcmp( nameB, info.nameB ) 
	  && !strcmp( nameC, info.nameA ) ) dup = 1;

      if(dup){
	sprintf( painCave.errMsg,
		 "Duplicate DUFF bend type \"%s - %s - %s\" found in "
		 "the DUFF param file./n",
		 nameA, nameB, nameC );
	painCave.isFatal = 1;
	simError();
      }

      if( next != NULL ) next->add(info);
      else{
	next = new LinkedBendType();
	strcpy(next->nameA, info.nameA);
	strcpy(next->nameB, info.nameB);
	strcpy(next->nameC, info.nameC);
	strcpy(next->type,  info.type);
	next->k1 = info.k1;
	next->k2 = info.k2;
	next->k3 = info.k3;
	next->t0 = info.t0;
      }
    }

#ifdef IS_MPI    

    void duplicate( bendStruct &info ){
      strcpy(info.nameA, nameA);
      strcpy(info.nameB, nameB);
      strcpy(info.nameC, nameC);
      strcpy(info.type,  type);
      info.k1   = k1;
      info.k2   = k2;
      info.k3   = k3;
      info.t0   = t0;
      info.last = 0;
    }

#endif // is_mpi

    char nameA[15];
    char nameB[15];
    char nameC[15];
    char type[30];
    double k1, k2, k3, t0;

    LinkedBendType* next;
  };

  class LinkedTorsionType {
  public:
    LinkedTorsionType(){ 
      next = NULL;
      nameA[0] = '\0';
      nameB[0] = '\0';
      nameC[0] = '\0';
      type[0] = '\0';
    }
    ~LinkedTorsionType(){ if( next != NULL ) delete next; }

    LinkedTorsionType* find( char* key1, char* key2, char* key3, char* key4 ){ 
      
 


      if( !strcmp( nameA, key1 ) && !strcmp( nameB, key2 ) &&
	  !strcmp( nameC, key3 ) && !strcmp( nameD, key4 ) ) return this;

      if( !strcmp( nameA, key4 ) && !strcmp( nameB, key3 ) &&
	  !strcmp( nameC, key2 ) && !strcmp( nameD, key1 ) ) return this;

      if( next != NULL ) return next->find(key1, key2, key3, key4);
      return NULL;
    }

    void add( torsionStruct &info ){

      // check for duplicates
      int dup = 0;

      if( !strcmp( nameA, info.nameA ) && !strcmp( nameB, info.nameB ) &&
	  !strcmp( nameC, info.nameC ) && !strcmp( nameD, info.nameD ) ) dup = 1;
      
      if( !strcmp( nameA, info.nameD ) && !strcmp( nameB, info.nameC ) &&
	  !strcmp( nameC, info.nameB ) && !strcmp( nameD, info.nameA ) ) dup = 1;
      
      if(dup){
	sprintf( painCave.errMsg,
		 "Duplicate DUFF torsion type \"%s - %s - %s - %s\" found in "
		 "the DUFF param file./n", nameA, nameB, nameC, nameD );
	painCave.isFatal = 1;
	simError();
      }

      if( next != NULL ) next->add(info);
      else{
	next = new LinkedTorsionType();
	strcpy(next->nameA, info.nameA);
	strcpy(next->nameB, info.nameB);
	strcpy(next->nameC, info.nameC);
	strcpy(next->nameD, info.nameD);
	strcpy(next->type,  info.type);
	next->k1 = info.k1;
	next->k2 = info.k2;
	next->k3 = info.k3;
	next->k4 = info.k4;

      }
    }

#ifdef IS_MPI
    
    void duplicate( torsionStruct &info ){
      strcpy(info.nameA, nameA);
      strcpy(info.nameB, nameB);
      strcpy(info.nameC, nameC);
      strcpy(info.nameD, nameD);
      strcpy(info.type,  type);
      info.k1   = k1;
      info.k2   = k2;
      info.k3   = k3;
      info.k4   = k4;
      info.last = 0;
    }

#endif

    char nameA[15];
    char nameB[15];
    char nameC[15];
    char nameD[15];
    char type[30];
    double k1, k2, k3, k4;

    LinkedTorsionType* next;
  };


  LinkedAtomType* headAtomType; 
  LinkedAtomType* currentAtomType;
  LinkedBondType* headBondType; 
  LinkedBondType* currentBondType;
  LinkedBendType* headBendType; 
  LinkedBendType* currentBendType;
  LinkedTorsionType* headTorsionType; 
  LinkedTorsionType* currentTorsionType;

} // namespace

using namespace DUFF_NS;


//****************************************************************
// begins the actual forcefield stuff.	
//****************************************************************


DUFF::DUFF(){

  char fileName[200];
  char* ffPath_env = "FORCE_PARAM_PATH";
  char* ffPath;
  char temp[200];

  headAtomType       = NULL; 
  currentAtomType    = NULL;
  headBondType       = NULL; 
  currentBondType    = NULL;
  headBendType       = NULL;
  currentBendType    = NULL;
  headTorsionType    = NULL; 
  currentTorsionType = NULL;


#ifdef IS_MPI
  int i;
  
  // **********************************************************************
  // Init the atomStruct mpi type

  atomStruct atomProto; // mpiPrototype
  int atomBC[3] = {15,12,5};  // block counts
  MPI_Aint atomDspls[3];           // displacements
  MPI_Datatype atomMbrTypes[3];    // member mpi types

  MPI_Address(&atomProto.name, &atomDspls[0]);
  MPI_Address(&atomProto.mass, &atomDspls[1]);
  MPI_Address(&atomProto.isSSD, &atomDspls[2]);
  
  atomMbrTypes[0] = MPI_CHAR;
  atomMbrTypes[1] = MPI_DOUBLE;
  atomMbrTypes[2] = MPI_INT;
  
  for (i=2; i >= 0; i--) atomDspls[i] -= atomDspls[0];
  
  MPI_Type_struct(3, atomBC, atomDspls, atomMbrTypes, &mpiAtomStructType);
  MPI_Type_commit(&mpiAtomStructType);


  // **********************************************************************
  // Init the bondStruct mpi type
  
  bondStruct bondProto; // mpiPrototype
  int bondBC[3] = {30,2,2};  // block counts
  MPI_Aint bondDspls[3];           // displacements
  MPI_Datatype bondMbrTypes[3];    // member mpi types
  
  MPI_Address(&bondProto.nameA, &bondDspls[0]);
  MPI_Address(&bondProto.d0,    &bondDspls[1]);
  MPI_Address(&bondProto.last,  &bondDspls[2]);
  
  bondMbrTypes[0] = MPI_CHAR;
  bondMbrTypes[1] = MPI_DOUBLE;
  bondMbrTypes[2] = MPI_INT;
  
  for (i=2; i >= 0; i--) bondDspls[i] -= bondDspls[0];
  
  MPI_Type_struct(3, bondBC, bondDspls, bondMbrTypes, &mpiBondStructType);
  MPI_Type_commit(&mpiBondStructType);


  // **********************************************************************
  // Init the bendStruct mpi type
  
  bendStruct bendProto; // mpiPrototype
  int bendBC[3] = {75,4,1};  // block counts
  MPI_Aint bendDspls[3];           // displacements
  MPI_Datatype bendMbrTypes[3];    // member mpi types
  
  MPI_Address(&bendProto.nameA, &bendDspls[0]);
  MPI_Address(&bendProto.k1,    &bendDspls[1]);
  MPI_Address(&bendProto.last,  &bendDspls[2]);
  
  bendMbrTypes[0] = MPI_CHAR;
  bendMbrTypes[1] = MPI_DOUBLE;
  bendMbrTypes[2] = MPI_INT;
  
  for (i=2; i >= 0; i--) bendDspls[i] -= bendDspls[0];
  
  MPI_Type_struct(3, bendBC, bendDspls, bendMbrTypes, &mpiBendStructType);
  MPI_Type_commit(&mpiBendStructType);


  // **********************************************************************
  // Init the torsionStruct mpi type
  
  torsionStruct torsionProto; // mpiPrototype
  int torsionBC[3] = {90,4,1};  // block counts
  MPI_Aint torsionDspls[3];           // displacements
  MPI_Datatype torsionMbrTypes[3];    // member mpi types
  
  MPI_Address(&torsionProto.nameA, &torsionDspls[0]);
  MPI_Address(&torsionProto.k1,    &torsionDspls[1]);
  MPI_Address(&torsionProto.last,  &torsionDspls[2]);
  
  torsionMbrTypes[0] = MPI_CHAR;
  torsionMbrTypes[1] = MPI_DOUBLE;
  torsionMbrTypes[2] = MPI_INT;
  
  for (i=2; i >= 0; i--) torsionDspls[i] -= torsionDspls[0];
  
  MPI_Type_struct(3, torsionBC, torsionDspls, torsionMbrTypes, 
		  &mpiTorsionStructType);
  MPI_Type_commit(&mpiTorsionStructType);

  // ***********************************************************************
  
  if( worldRank == 0 ){
#endif
    
    // generate the force file name
    
    strcpy( fileName, "DUFF.frc" ); 
    //    fprintf( stderr,"Trying to open %s\n", fileName );
    
    // attempt to open the file in the current directory first.
    
    frcFile = fopen( fileName, "r" );
    
    if( frcFile == NULL ){
      
      // next see if the force path enviorment variable is set
      
      ffPath = getenv( ffPath_env );
      if( ffPath == NULL ) {
	STR_DEFINE(ffPath, FRC_PATH );
      }
      
      
      strcpy( temp, ffPath );
      strcat( temp, "/" );
      strcat( temp, fileName );
      strcpy( fileName, temp );
      
      frcFile = fopen( fileName, "r" );
      
      if( frcFile == NULL ){
	
	sprintf( painCave.errMsg,
		 "Error opening the force field parameter file:\n"
                 "\t%s\n"
		 "\tHave you tried setting the FORCE_PARAM_PATH environment "
		 "variable?\n",
		 fileName );
        painCave.severity = OOPSE_ERROR;
	painCave.isFatal = 1;
	simError();
      }
    }
    
#ifdef IS_MPI
  }
  
  sprintf( checkPointMsg, "DUFF file opened sucessfully." );
  MPIcheckPoint();
  
#endif // is_mpi
}


DUFF::~DUFF(){

  if( headAtomType != NULL ) delete headAtomType;
  if( headBondType != NULL ) delete headBondType;
  if( headBendType != NULL ) delete headBendType;
  if( headTorsionType != NULL ) delete headTorsionType;

#ifdef IS_MPI
  if( worldRank == 0 ){
#endif // is_mpi
    
    fclose( frcFile );
    
#ifdef IS_MPI
  }
#endif // is_mpi
}

void DUFF::cleanMe( void ){

#ifdef IS_MPI
  
  // keep the linked lists in the mpi version

#else // is_mpi

  // delete the linked lists in the single processor version

  if( headAtomType != NULL ) delete headAtomType;
  if( headBondType != NULL ) delete headBondType;
  if( headBendType != NULL ) delete headBendType;
  if( headTorsionType != NULL ) delete headTorsionType;

#endif // is_mpi
}


void DUFF::initForceField( int ljMixRule ){
  
  initFortran( ljMixRule, entry_plug->useReactionField );
}

double DUFF::getAtomTypeMass (char* atomType) {

  currentAtomType = headAtomType->find( atomType );
  if( currentAtomType == NULL ){
    sprintf( painCave.errMsg,
            "AtomType error, %s not found in force file.\n",
             atomType );
    painCave.isFatal = 1;
    simError();
  }

  return currentAtomType->mass;
}

void DUFF::readParams( void ){

  int identNum;
  
  atomStruct atomInfo;
  bondStruct bondInfo;
  bendStruct bendInfo;
  torsionStruct torsionInfo;
  
  bigSigma = 0.0;

  atomInfo.last = 1; 
  bondInfo.last = 1;
  bendInfo.last = 1;
  torsionInfo.last = 1;

  // read in the atom info
  
#ifdef IS_MPI
  if( worldRank == 0 ){
#endif
    
    // read in the atom types.
    
    headAtomType = new LinkedAtomType;
    
    fastForward( "AtomTypes", "initializeAtoms" );

    // we are now at the AtomTypes section.
    
    eof_test = fgets( readLine, sizeof(readLine), frcFile );
    lineNum++;
    
    
    // read a line, and start parseing out the atom types 

    if( eof_test == NULL ){
      sprintf( painCave.errMsg, 
	       "Error in reading Atoms from force file at line %d.\n",
	       lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    identNum = 1;
    // stop reading at end of file, or at next section
    while( readLine[0] != '#' && eof_test != NULL ){

      // toss comment lines
      if( readLine[0] != '!' ){
	
	// the parser returns 0 if the line was blank
	if( parseAtom( readLine, lineNum, atomInfo ) ){
	  atomInfo.ident = identNum;
	  headAtomType->add( atomInfo );;
	  identNum++;
	}
      }
      eof_test = fgets( readLine, sizeof(readLine), frcFile );
      lineNum++;
    }

#ifdef IS_MPI
    
    // send out the linked list to all the other processes

    sprintf( checkPointMsg,
	     "DUFF atom structures read successfully." );
    MPIcheckPoint();

    currentAtomType = headAtomType->next; //skip the first element who is a place holder.
    while( currentAtomType != NULL ){
      currentAtomType->duplicate( atomInfo );

      sendFrcStruct( &atomInfo, mpiAtomStructType );

      sprintf( checkPointMsg,
	       "successfully sent DUFF force type: \"%s\"\n",
	       atomInfo.name );
      MPIcheckPoint();

      currentAtomType = currentAtomType->next;
    }
    atomInfo.last = 1;
    sendFrcStruct( &atomInfo, mpiAtomStructType );
    
  }

  else{
    
    // listen for node 0 to send out the force params

    MPIcheckPoint();

    headAtomType = new LinkedAtomType;
    receiveFrcStruct( &atomInfo, mpiAtomStructType );
    
    while( !atomInfo.last ){

      headAtomType->add( atomInfo );
      
      MPIcheckPoint();

      receiveFrcStruct( &atomInfo, mpiAtomStructType );
    }
  }

#endif // is_mpi



  // call new A_types in fortran
  
  int isError;
  
  // dummy variables
  
  int isGB = 0;
  int isLJ = 1;
  int isEAM =0;
  int isCharge = 0;
  double charge=0.0;
    
  currentAtomType = headAtomType->next;;
  while( currentAtomType != NULL ){
    
    if(currentAtomType->isDipole) entry_plug->useDipoles = 1;
    if(currentAtomType->isSSD) {
      entry_plug->useSticky = 1;
      set_sticky_params( &(currentAtomType->w0), &(currentAtomType->v0), 
                         &(currentAtomType->v0p), 
                         &(currentAtomType->rl), &(currentAtomType->ru), 
                         &(currentAtomType->rlp), &(currentAtomType->rup));
    }

    if( currentAtomType->name[0] != '\0' ){
      isError = 0;
      makeAtype( &(currentAtomType->ident),
		 &isLJ,
		 &(currentAtomType->isSSD),
		 &(currentAtomType->isDipole),
		 &isGB,
		 &isEAM,
		 &isCharge,
		 &(currentAtomType->epslon),
		 &(currentAtomType->sigma),
		 &charge,
		 &(currentAtomType->dipole),
		 &isError );
      if( isError ){
	sprintf( painCave.errMsg,
		 "Error initializing the \"%s\" atom type in fortran\n",
		 currentAtomType->name );
	painCave.isFatal = 1;
	simError();
      }
    }
    currentAtomType = currentAtomType->next;
  }
      
#ifdef IS_MPI
  sprintf( checkPointMsg,
	   "DUFF atom structures successfully sent to fortran\n" );
  MPIcheckPoint();
#endif // is_mpi

  

  // read in the bonds
  
#ifdef IS_MPI
  if( worldRank == 0 ){
#endif
    
    // read in the bond types.
    
    headBondType = new LinkedBondType;
    
    fastForward( "BondTypes", "initializeBonds" );

    // we are now at the bondTypes section

    eof_test =  fgets( readLine, sizeof(readLine), frcFile );
    lineNum++;
    
    
    // read a line, and start parseing out the atom types 

    if( eof_test == NULL ){
      sprintf( painCave.errMsg, 
	       "Error in reading bonds from force file at line %d.\n",
	       lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    // stop reading at end of file, or at next section
    while( readLine[0] != '#' && eof_test != NULL ){

      // toss comment lines
      if( readLine[0] != '!' ){
	
	// the parser returns 0 if the line was blank
	if( parseBond( readLine, lineNum, bondInfo ) ){
	  headBondType->add( bondInfo );
	}
      }
      eof_test = fgets( readLine, sizeof(readLine), frcFile );
      lineNum++;
    }
        
#ifdef IS_MPI
    
    // send out the linked list to all the other processes
    
    sprintf( checkPointMsg,
	     "DUFF bond structures read successfully." );
    MPIcheckPoint();
    
    currentBondType = headBondType->next;
    while( currentBondType != NULL ){
      currentBondType->duplicate( bondInfo );
      sendFrcStruct( &bondInfo, mpiBondStructType );
      currentBondType = currentBondType->next;
    }
    bondInfo.last = 1;
    sendFrcStruct( &bondInfo, mpiBondStructType );
    
  }

  else{
    
    // listen for node 0 to send out the force params
    
    MPIcheckPoint(); 

    headBondType = new LinkedBondType;
    receiveFrcStruct( &bondInfo, mpiBondStructType );
    while( !bondInfo.last ){

      headBondType->add( bondInfo );
      receiveFrcStruct( &bondInfo, mpiBondStructType );
    }
  }

  sprintf( checkPointMsg,
	   "DUFF bond structures broadcast successfully." );
  MPIcheckPoint();

#endif // is_mpi
  

  // read in the bends

#ifdef IS_MPI
  if( worldRank == 0 ){
#endif

    // read in the bend types.

    headBendType = new LinkedBendType;
    
    fastForward( "BendTypes", "initializeBends" );

    // we are now at the bendTypes section
    
    eof_test =  fgets( readLine, sizeof(readLine), frcFile );
    lineNum++;
        
    // read a line, and start parseing out the bend types 

    if( eof_test == NULL ){
      sprintf( painCave.errMsg, 
	       "Error in reading bends from force file at line %d.\n",
	       lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    // stop reading at end of file, or at next section
    while( readLine[0] != '#' && eof_test != NULL ){
      
      // toss comment lines
      if( readLine[0] != '!' ){
	
	// the parser returns 0 if the line was blank
	if( parseBend( readLine, lineNum, bendInfo ) ){
	  headBendType->add( bendInfo );
	}
      }
      eof_test = fgets( readLine, sizeof(readLine), frcFile );
      lineNum++;
    }
    
#ifdef IS_MPI
    
    // send out the linked list to all the other processes

    sprintf( checkPointMsg,
	     "DUFF bend structures read successfully." );
    MPIcheckPoint();

    currentBendType = headBendType->next;
    while( currentBendType != NULL ){
      currentBendType->duplicate( bendInfo );
      sendFrcStruct( &bendInfo, mpiBendStructType );
      currentBendType = currentBendType->next;
    }
    bendInfo.last = 1;
    sendFrcStruct( &bendInfo, mpiBendStructType );
    
  }

  else{
    
    // listen for node 0 to send out the force params
    
    MPIcheckPoint();

    headBendType = new LinkedBendType;
    receiveFrcStruct( &bendInfo, mpiBendStructType );
    while( !bendInfo.last ){

      headBendType->add( bendInfo );
      receiveFrcStruct( &bendInfo, mpiBendStructType );
    }
  }

  sprintf( checkPointMsg,
	   "DUFF bend structures broadcast successfully." );
  MPIcheckPoint();

#endif // is_mpi


  // read in the torsions

#ifdef IS_MPI
  if( worldRank == 0 ){
#endif

    // read in the torsion types.

    headTorsionType = new LinkedTorsionType;
    
    fastForward( "TorsionTypes", "initializeTorsions" );

    // we are now at the torsionTypes section

    eof_test =  fgets( readLine, sizeof(readLine), frcFile );
    lineNum++;
    
    
    // read a line, and start parseing out the atom types 

    if( eof_test == NULL ){
      sprintf( painCave.errMsg, 
	       "Error in reading torsions from force file at line %d.\n",
	       lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    // stop reading at end of file, or at next section
    while( readLine[0] != '#' && eof_test != NULL ){

      // toss comment lines
      if( readLine[0] != '!' ){
	
	// the parser returns 0 if the line was blank
	if( parseTorsion( readLine, lineNum, torsionInfo ) ){
	  headTorsionType->add( torsionInfo );

	}
      }
      eof_test = fgets( readLine, sizeof(readLine), frcFile );
      lineNum++;
    }

#ifdef IS_MPI
    
    // send out the linked list to all the other processes
    
    sprintf( checkPointMsg,
	     "DUFF torsion structures read successfully." );
    MPIcheckPoint();
    
    currentTorsionType = headTorsionType->next;
    while( currentTorsionType != NULL ){
      currentTorsionType->duplicate( torsionInfo );
      sendFrcStruct( &torsionInfo, mpiTorsionStructType );
      currentTorsionType = currentTorsionType->next;
    }
    torsionInfo.last = 1;
    sendFrcStruct( &torsionInfo, mpiTorsionStructType );
    
  }

  else{
    
    // listen for node 0 to send out the force params
    
    MPIcheckPoint();

    headTorsionType = new LinkedTorsionType;
    receiveFrcStruct( &torsionInfo, mpiTorsionStructType );
    while( !torsionInfo.last ){

      headTorsionType->add( torsionInfo );
      receiveFrcStruct( &torsionInfo, mpiTorsionStructType );
    }
  }

  sprintf( checkPointMsg,
	   "DUFF torsion structures broadcast successfully." );
  MPIcheckPoint();

#endif // is_mpi

  entry_plug->useLJ = 1;
}



void DUFF::initializeAtoms( int nAtoms, Atom** the_atoms ){
  
  
  //////////////////////////////////////////////////
  // a quick water fix 

  double waterI[3][3]; 
  waterI[0][0] = 1.76958347772500;
  waterI[0][1] = 0.0;
  waterI[0][2] = 0.0;

  waterI[1][0] = 0.0;
  waterI[1][1] = 0.614537057924513;
  waterI[1][2] = 0.0;

  waterI[2][0] = 0.0;
  waterI[2][1] = 0.0;
  waterI[2][2] = 1.15504641980049;


  double headI[3][3]; 
  headI[0][0] = 1125;
  headI[0][1] = 0.0;
  headI[0][2] = 0.0;

  headI[1][0] = 0.0;
  headI[1][1] = 1125;
  headI[1][2] = 0.0;

  headI[2][0] = 0.0;
  headI[2][1] = 0.0;
  headI[2][2] = 250;

  //////////////////////////////////////////////////

  
  // initialize the atoms
  
  DirectionalAtom* dAtom;
  double ji[3];

  for(int i=0; i<nAtoms; i++ ){

    currentAtomType = headAtomType->find( the_atoms[i]->getType() );
    if( currentAtomType == NULL ){
      sprintf( painCave.errMsg, 
	       "AtomType error, %s not found in force file.\n",
	       the_atoms[i]->getType() );
      painCave.isFatal = 1;
      simError();
    }
    
    the_atoms[i]->setMass( currentAtomType->mass );
    the_atoms[i]->setIdent( currentAtomType->ident );

    if( bigSigma < currentAtomType->sigma ) bigSigma = currentAtomType->sigma;

    if( currentAtomType->isDipole ){
      if( the_atoms[i]->isDirectional() ){
	
	dAtom = (DirectionalAtom *) the_atoms[i];
	dAtom->setHasDipole( 1 );

        ji[0] = 0.0;
        ji[1] = 0.0;
        ji[2] = 0.0;

        dAtom->setJ( ji );
	
	if(!strcmp("SSD",the_atoms[i]->getType())){
	  dAtom->setI( waterI );
	}
	else if(!strcmp("HEAD",the_atoms[i]->getType())){
	  dAtom->setI( headI );
	}
	else{
	  sprintf(painCave.errMsg,
		  "AtmType error, %s does not have a moment of inertia set.\n",
		  the_atoms[i]->getType() );
	  painCave.isFatal = 1;
	  simError();
	}
	entry_plug->n_dipoles++;
      }
      else{
	
	sprintf( painCave.errMsg,
		"DUFF error: Atom \"%s\" is a dipole, yet no standard"
		 " orientation was specifed in the meta-data file.\n",
		 currentAtomType->name );
	painCave.isFatal = 1;
	simError();
      }
    }
    else{
      if( the_atoms[i]->isDirectional() ){
	sprintf( painCave.errMsg,
		 "DUFF error: Atom \"%s\" was given a standard "
		 "orientation in the meta-data file, yet it is not a dipole.\n",
		 currentAtomType->name);
	painCave.isFatal = 1;
	simError();
      }
    }
  }
}

void DUFF::initializeBonds( int nBonds, Bond** bondArray,
				   bond_pair* the_bonds ){
  int i,a,b;
  char* atomA;
  char* atomB;
  
  Atom** the_atoms;
  the_atoms = entry_plug->atoms;
  

  // initialize the Bonds
  
  for( i=0; i<nBonds; i++ ){
    
    a = the_bonds[i].a;
    b = the_bonds[i].b;

    atomA = the_atoms[a]->getType();
    atomB = the_atoms[b]->getType();
    currentBondType = headBondType->find( atomA, atomB );
    if( currentBondType == NULL ){
      sprintf( painCave.errMsg,
	       "BondType error, %s - %s not found in force file.\n",
	       atomA, atomB );
      painCave.isFatal = 1;
      simError();
    }
    
    switch( currentBondType->type ){

    case FIXED_BOND:
            
      bondArray[i] = new ConstrainedBond( *the_atoms[a], 
					  *the_atoms[b],
					  currentBondType->d0 );
      entry_plug->n_constraints++;
      break;

    case HARMONIC_BOND:
      
      bondArray[i] = new HarmonicBond( *the_atoms[a],
				       *the_atoms[b],
				       currentBondType->d0,
				       currentBondType->k0 );
      break;
      
    default:

      break;
      // do nothing
    }
  }
}

void DUFF::initializeBends( int nBends, Bend** bendArray,
				   bend_set* the_bends ){
  
  QuadraticBend* qBend;
  GhostBend* gBend;
  Atom** the_atoms;
  the_atoms = entry_plug->atoms;
  
  int i, a, b, c;
  char* atomA;
  char* atomB;
  char* atomC;
  
  // initialize the Bends

  for( i=0; i<nBends; i++ ){
    
    a = the_bends[i].a;
    b = the_bends[i].b;
    c = the_bends[i].c;

    atomA = the_atoms[a]->getType();
    atomB = the_atoms[b]->getType();

    if( the_bends[i].isGhost ) atomC = "GHOST";
    else atomC = the_atoms[c]->getType();

    currentBendType = headBendType->find( atomA, atomB, atomC );
    if( currentBendType == NULL ){
      sprintf( painCave.errMsg, "BendType error, %s - %s - %s not found"
	       " in force file.\n",
	       atomA, atomB, atomC );
      painCave.isFatal = 1;
      simError();
    }
    
    if( !strcmp( currentBendType->type, "quadratic" ) ){
      
      if( the_bends[i].isGhost){
	
	if( the_bends[i].ghost == b ){
	  // do nothing
	}
	else if( the_bends[i].ghost == a ){
	  c = a;
	  a = b;
	  b = c;
	}
	else{ 
	  sprintf( painCave.errMsg, 
		   "BendType error, %s - %s - %s,\n"
		   "  --> central atom is not "
		   "correctly identified with the "
		   "\"ghostVectorSource = \" tag.\n",
		   atomA, atomB, atomC );
	  painCave.isFatal = 1;
	  simError();
	}
	
	gBend = new GhostBend( *the_atoms[a], 
			       *the_atoms[b]);
			       			      		       
	gBend->setConstants( currentBendType->k1,
			     currentBendType->k2,
			     currentBendType->k3,
			     currentBendType->t0 );
	bendArray[i] = gBend;
      }
      else{
	qBend = new QuadraticBend( *the_atoms[a], 
				   *the_atoms[b],
				   *the_atoms[c] );
	qBend->setConstants( currentBendType->k1,
			     currentBendType->k2,
			     currentBendType->k3,
			     currentBendType->t0 );
	bendArray[i] = qBend;
      }      
    }
  }
}

void DUFF::initializeTorsions( int nTorsions, Torsion** torsionArray,
				      torsion_set* the_torsions ){

  int i, a, b, c, d;
  char* atomA;
  char* atomB;
  char* atomC;
  char* atomD;

  CubicTorsion* cTors;
  Atom** the_atoms;
  the_atoms = entry_plug->atoms;

  // initialize the Torsions

  for( i=0; i<nTorsions; i++ ){
    
    a = the_torsions[i].a;
    b = the_torsions[i].b;
    c = the_torsions[i].c;
    d = the_torsions[i].d;

    atomA = the_atoms[a]->getType();
    atomB = the_atoms[b]->getType();
    atomC = the_atoms[c]->getType();
    atomD = the_atoms[d]->getType();
    currentTorsionType = headTorsionType->find( atomA, atomB, atomC, atomD );
    if( currentTorsionType == NULL ){
      sprintf( painCave.errMsg,
	       "TorsionType error, %s - %s - %s - %s not found"
	       " in force file.\n",
	       atomA, atomB, atomC, atomD );
      painCave.isFatal = 1;
      simError();
    }
    
    if( !strcmp( currentTorsionType->type, "cubic" ) ){
      
      cTors = new CubicTorsion( *the_atoms[a], *the_atoms[b],
				*the_atoms[c], *the_atoms[d] );
      cTors->setConstants( currentTorsionType->k1, currentTorsionType->k2,
			   currentTorsionType->k3, currentTorsionType->k4 );
      torsionArray[i] = cTors;
    }
  }
}

void DUFF::fastForward( char* stopText, char* searchOwner ){

  int foundText = 0;
  char* the_token;

  rewind( frcFile );
  lineNum = 0;

  eof_test = fgets( readLine, sizeof(readLine), frcFile );
  lineNum++;
  if( eof_test == NULL ){
    sprintf( painCave.errMsg, "Error fast forwarding force file for %s: "
	     " file is empty.\n",
	     searchOwner );
    painCave.isFatal = 1;
    simError();
  }
  
  
  while( !foundText ){
    while( eof_test != NULL && readLine[0] != '#' ){
      eof_test = fgets( readLine, sizeof(readLine), frcFile );
      lineNum++;
    }
    if( eof_test == NULL ){
      sprintf( painCave.errMsg,
	       "Error fast forwarding force file for %s at "
	       "line %d: file ended unexpectedly.\n",
	       searchOwner,
	       lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    the_token = strtok( readLine, " ,;\t#\n" );
    foundText = !strcmp( stopText, the_token );
    
    if( !foundText ){
      eof_test = fgets( readLine, sizeof(readLine), frcFile );
      lineNum++;
      
      if( eof_test == NULL ){
	sprintf( painCave.errMsg, 
		 "Error fast forwarding force file for %s at "
		 "line %d: file ended unexpectedly.\n",
		 searchOwner,
		 lineNum );
	painCave.isFatal = 1;
	simError();
      } 
    }
  }  
}


int DUFF_NS::parseAtom( char *lineBuffer, int lineNum, atomStruct &info ){
 
  char* the_token;
  
  the_token = strtok( lineBuffer, " \n\t,;" );
  if( the_token != NULL ){
    
    strcpy( info.name, the_token );

    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing AtomTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    info.isDipole = atoi( the_token );

    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing AtomTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }

    info.isSSD = atoi( the_token );

    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing AtomTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    info.mass = atof( the_token );
    
    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing AtomTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }
	
    info.epslon = atof( the_token );
	  
    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing AtomTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }
	
    info.sigma = atof( the_token );

    if( info.isDipole ){
      
      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg, 
		 "Error parseing AtomTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.dipole = atof( the_token );
    }
    else info.dipole = 0.0;

    if( info.isSSD ){

      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg, 
		 "Error parseing AtomTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.w0 = atof( the_token );

      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg, 
		 "Error parseing AtomTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.v0 = atof( the_token );
      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg, 
		 "Error parseing AtomTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.v0p = atof( the_token );

      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg, 
		 "Error parseing AtomTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.rl = atof( the_token );

      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg, 
		 "Error parseing AtomTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.ru = atof( the_token );

      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg, 
		 "Error parseing AtomTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.rlp = atof( the_token );

      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg, 
		 "Error parseing AtomTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.rup = atof( the_token );
    }
    else info.v0 = info.w0 = info.v0p = info.rl = info.ru = info.rlp = info.rup = 0.0;

    return 1;
  }
  else return 0;
}

int DUFF_NS::parseBond( char *lineBuffer, int lineNum, bondStruct &info ){
 
  char* the_token;
  char bondType[30];
  
  the_token = strtok( lineBuffer, " \n\t,;" );
  if( the_token != NULL ){
    
    strcpy( info.nameA, the_token );

    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing BondTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    strcpy( info.nameB, the_token );

    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing BondTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    strcpy( bondType, the_token );
    
    if( !strcmp( bondType, "fixed" ) ){
      info.type = FIXED_BOND;
      
      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg,
		 "Error parseing BondTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.d0 = atof( the_token );
      
      info.k0=0.0;
    }
    else if( !strcmp( bondType, "harmonic" ) ){
      info.type = HARMONIC_BOND;
      
      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg,
		 "Error parseing BondTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.d0 = atof( the_token );

      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg,
		 "Error parseing BondTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.k0 = atof( the_token );
    }

    else{
      sprintf( painCave.errMsg, 
	       "Unknown DUFF bond type \"%s\" at line %d\n",
	       bondType,
	       lineNum );
      painCave.isFatal = 1;
      simError();
    }		  
    
    return 1;
  }
  else return 0;
}


int DUFF_NS::parseBend( char *lineBuffer, int lineNum, bendStruct &info ){

  char* the_token;
  
  the_token = strtok( lineBuffer, " \n\t,;" );
  if( the_token != NULL ){
    
    strcpy( info.nameA, the_token );

    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing BendTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    strcpy( info.nameB, the_token );

    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing BendTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    strcpy( info.nameC, the_token );

    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing BendTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    strcpy( info.type, the_token );

    if( !strcmp( info.type, "quadratic" ) ){
      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg,
		 "Error parseing BendTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
	
      info.k1 = atof( the_token );
      
      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg,
		 "Error parseing BendTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.k2 = atof( the_token );
      
      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg,
		 "Error parseing BendTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
	
      info.k3 = atof( the_token );
      
      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg,
		 "Error parseing BendTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.t0 = atof( the_token );
    }
    
    else{
      sprintf( painCave.errMsg, 
	       "Unknown DUFF bend type \"%s\" at line %d\n",
	       info.type,
	       lineNum );
      painCave.isFatal = 1;
      simError();
    }		  
        
    return 1;
  }
  else return 0;
}

int DUFF_NS::parseTorsion( char *lineBuffer, int lineNum, torsionStruct &info ){
  
  char*  the_token;

  the_token = strtok( lineBuffer, " \n\t,;" );
  if( the_token != NULL ){
    
    strcpy( info.nameA, the_token );
	
    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg,
	       "Error parseing TorsionTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    strcpy( info.nameB, the_token );

    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg,
	       "Error parseing TorsionTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    strcpy( info.nameC, the_token );
    
    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg,
	       "Error parseing TorsionTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    strcpy( info.nameD, the_token );
    
    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg,
	       "Error parseing TorsionTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    strcpy( info.type, the_token );
    
    if( !strcmp( info.type, "cubic" ) ){
      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg,
		 "Error parseing TorsionTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.k1 = atof( the_token );
      
      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg,
		 "Error parseing TorsionTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.k2 = atof( the_token );
      
      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg,
		 "Error parseing TorsionTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.k3 = atof( the_token );
      
      if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
	sprintf( painCave.errMsg,
		 "Error parseing TorsionTypes: line %d\n", lineNum );
	painCave.isFatal = 1;
	simError();
      }
      
      info.k4 = atof( the_token );
    
    }
    
    else{
      sprintf( painCave.errMsg, 
	       "Unknown DUFF torsion type \"%s\" at line %d\n",
	       info.type,
	       lineNum );
      painCave.isFatal = 1;
      simError();
    }		  
    
    return 1;
  }
  
  else return 0;
}
