#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
using namespace std;

#ifdef IS_MPI
#include <mpi.h>
#endif //is_mpi

#include "ForceFields.hpp"
#include "SRI.hpp"
#include "simError.h"

#include "fortranWrappers.hpp"

#ifdef IS_MPI
#include "mpiForceField.h"
#endif // is_mpi



namespace LJ_NS{

  // Declare the structures that will be passed by the parser and  MPI
  
  typedef struct{
    char name[15];
    double mass;
    double epslon;
    double sigma;
    int ident;
    int last;      //  0  -> default
                   //  1  -> in MPI: tells nodes to stop listening
  } atomStruct;

  int parseAtom( char *lineBuffer, int lineNum, atomStruct &info );
  
#ifdef IS_MPI
  
  MPI_Datatype mpiAtomStructType;
  
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
    

    void add( atomStruct &info ){
    
      // check for duplicates
      
      if( !strcmp( info.name, name ) ){
	sprintf( painCave.errMsg,
		 "Duplicate LJ atom type \"%s\" found in "
		 "the LJFF param file./n",
		 name );
	painCave.isFatal = 1;
	simError();
      }
      
      if( next != NULL ) next->add(info);
      else{
	next = new LinkedAtomType();
	strcpy(next->name, info.name);
	next->mass     = info.mass;
	next->epslon   = info.epslon;
	next->sigma    = info.sigma;
	next->ident    = info.ident;
      }
    }
    

#ifdef IS_MPI
    
    void duplicate( atomStruct &info ){
      strcpy(info.name, name);
      info.mass     = mass;
      info.epslon   = epslon;
      info.sigma    = sigma;
      info.ident    = ident;
      info.last     = 0;
    }


#endif

    char name[15];
    double mass;
    double epslon;
    double sigma;
    int ident;
    LinkedAtomType* next;
  };

  LinkedAtomType* headAtomType; 
  LinkedAtomType* currentAtomType;

}

using namespace LJ_NS;

//****************************************************************
// begins the actual forcefield stuff.	
//****************************************************************


LJFF::LJFF(){

  char fileName[200];
  char* ffPath_env = "FORCE_PARAM_PATH";
  char* ffPath;
  char temp[200];

  headAtomType = NULL;
  currentAtomType = NULL;

  // do the funtion wrapping
  wrapMeFF( this );

#ifdef IS_MPI
  int i;
  
   // **********************************************************************
  // Init the atomStruct mpi type

  atomStruct atomProto; // mpiPrototype
  int atomBC[3] = {15,3,2};  // block counts
  MPI_Aint atomDspls[3];           // displacements
  MPI_Datatype atomMbrTypes[3];    // member mpi types

  MPI_Address(&atomProto.name, &atomDspls[0]);
  MPI_Address(&atomProto.mass, &atomDspls[1]);
  MPI_Address(&atomProto.ident, &atomDspls[2]);
  
  atomMbrTypes[0] = MPI_CHAR;
  atomMbrTypes[1] = MPI_DOUBLE;
  atomMbrTypes[2] = MPI_INT;
  
  for (i=2; i >= 0; i--) atomDspls[i] -= atomDspls[0];
  
  MPI_Type_struct(3, atomBC, atomDspls, atomMbrTypes, &mpiAtomStructType);
  MPI_Type_commit(&mpiAtomStructType);

  // ***********************************************************************
  
  if( worldRank == 0 ){
#endif
    
    // generate the force file name
    
    strcpy( fileName, "LJFF.frc" ); 
    // fprintf( stderr,"Trying to open %s\n", fileName );
    
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
  
  sprintf( checkPointMsg, "LJFF file opened sucessfully." );
  MPIcheckPoint();
  
#endif // is_mpi
}


LJFF::~LJFF(){

  if( headAtomType != NULL ) delete headAtomType;

#ifdef IS_MPI
  if( worldRank == 0 ){
#endif // is_mpi
    
    fclose( frcFile );
    
#ifdef IS_MPI
  }
#endif // is_mpi
}

void LJFF::initForceField( int ljMixRule ){
  
  initFortran( ljMixRule, 0 );
}

void LJFF::cleanMe( void ){

#ifdef IS_MPI
  
  // keep the linked list in the mpi version

#else // is_mpi

  // delete the linked list in the single processor version

  if( headAtomType != NULL ) delete headAtomType;

#endif // is_mpi
}

void LJFF::readParams( void ){

  atomStruct info;
  info.last = 1; // initialize last to have the last set. 
                 // if things go well, last will be set to 0

  int identNum;
  

  bigSigma = 0.0;
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
	if( parseAtom( readLine, lineNum, info ) ){
	  info.ident = identNum;
	  headAtomType->add( info );;
	  identNum++;
	}
      }
      eof_test = fgets( readLine, sizeof(readLine), frcFile );
      lineNum++;
    }

#ifdef IS_MPI
    
    // send out the linked list to all the other processes

    sprintf( checkPointMsg,
	     "LJFF atom structures read successfully." );
    MPIcheckPoint();

    currentAtomType = headAtomType->next; //skip the first element who is a place holder.
    while( currentAtomType != NULL ){
      currentAtomType->duplicate( info );

 

      sendFrcStruct( &info, mpiAtomStructType );

      sprintf( checkPointMsg,
	       "successfully sent lJ force type: \"%s\"\n",
	       info.name );
      MPIcheckPoint();

      currentAtomType = currentAtomType->next;
    }
    info.last = 1;
    sendFrcStruct( &info, mpiAtomStructType );
    
  }

  else{
    
    // listen for node 0 to send out the force params
    
    MPIcheckPoint();

    headAtomType = new LinkedAtomType;
    receiveFrcStruct( &info, mpiAtomStructType );
    
    while( !info.last ){



      headAtomType->add( info );
      
      MPIcheckPoint();

      receiveFrcStruct( &info, mpiAtomStructType );
    }
  }
#endif // is_mpi

  // call new A_types in fortran
  
  int isError;

  // dummy variables
  int isLJ = 1;
  int isDipole = 0;
  int isSSD = 0;
  int isGB = 0;
  int isEAM = 0;
  int isCharge = 0;
  double charge = 0.0;
  double dipole = 0.0;
  
  currentAtomType = headAtomType;
  while( currentAtomType != NULL ){
    
    if( currentAtomType->name[0] != '\0' ){
      isError = 0;
      makeAtype( &(currentAtomType->ident),
		 &isLJ,
		 &isSSD,
		 &isDipole,
		 &isGB,
		 &isEAM,
                 &isCharge,
		 &(currentAtomType->epslon),
		 &(currentAtomType->sigma),
                 &charge,
		 &dipole,
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
      
  entry_plug->useLJ = 1;

#ifdef IS_MPI
  sprintf( checkPointMsg,
	   "LJFF atom structures successfully sent to fortran\n" );
  MPIcheckPoint();
#endif // is_mpi

}

double LJFF::getAtomTypeMass (char* atomType) {

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

void LJFF::initializeAtoms( int nAtoms, Atom** the_atoms ){
  
  int i;

  // initialize the atoms
  

  for( i=0; i<nAtoms; i++ ){
    
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
  }
}

void LJFF::initializeBonds( int nBonds, Bond** BondArray,
			     bond_pair* the_bonds ){
  
    if( nBonds ){
      sprintf( painCave.errMsg,
	       "LJFF does not support bonds.\n" );
      painCave.isFatal = 1;
      simError();
    }
}

void LJFF::initializeBends( int nBends, Bend** bendArray,
			     bend_set* the_bends ){

    if( nBends ){
      sprintf( painCave.errMsg,
	       "LJFF does not support bends.\n" );
      painCave.isFatal = 1;
      simError();
    }
}

void LJFF::initializeTorsions( int nTorsions, Torsion** torsionArray,
				torsion_set* the_torsions ){

    if( nTorsions ){
      sprintf( painCave.errMsg,
	       "LJFF does not support torsions.\n" );
      painCave.isFatal = 1;
      simError();
    }
}

void LJFF::fastForward( char* stopText, char* searchOwner ){

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



int LJ_NS::parseAtom( char *lineBuffer, int lineNum,  atomStruct &info ){

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
    
    return 1;
  }
  else return 0;
}


