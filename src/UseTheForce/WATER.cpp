#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
using namespace std;

#ifdef IS_MPI
#include <mpi.h>
#endif //is_mpi

#include "UseTheForce/ForceFields.hpp"
#include "primitives/SRI.hpp"
#include "utils/simError.h"
#include "types/AtomType.hpp"
#include "types/DirectionalAtomType.hpp"
#include "UseTheForce/DarkSide/lj_interface.h"
#include "UseTheForce/DarkSide/charge_interface.h"
#include "UseTheForce/DarkSide/dipole_interface.h"
#include "UseTheForce/DarkSide/sticky_interface.h"

#ifdef IS_MPI
#include "UseTheForce/mpiForceField.h"
#endif // is_mpi



namespace WATER_NS{

  // Declare the structures that will be passed by the parser and MPI
  
  typedef struct{
    char name[15];
    double mass;
    double epslon;
    double sigma;
    double charge;
    int isDirectional;
    int isLJ;
    int isCharge;
    int ident;
    int last;      //  0  -> default
                   //  1  -> in MPI: tells nodes to stop listening
  } atomStruct;

  typedef struct{
    char name[15];
    double Ixx;
    double Iyy;
    double Izz;
    double dipole;
    double w0;
    double v0;
    double v0p;
    double rl;
    double ru;
    double rlp;
    double rup;
    int isDipole;
    int isSticky;
    int last;      //  0  -> default
                   //  1  -> in MPI: tells nodes to stop listening
  } directionalStruct;

  int parseAtom( char *lineBuffer, int lineNum, atomStruct &info );
  int parseDirectional( char *lineBuffer, int lineNum, directionalStruct &info );


#ifdef IS_MPI
  
  MPI_Datatype mpiAtomStructType;
  MPI_Datatype mpiDirectionalStructType;

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
		 "Duplicate WATER atom type \"%s\" found in "
		 "the WATER param file./n",
		 name );
	painCave.isFatal = 1;
	simError();
      }
      
      if( next != NULL ) next->add(info);
      else{
	next = new LinkedAtomType();
	strcpy(next->name, info.name);
	next->isDirectional  = info.isDirectional;
	next->isLJ     = info.isLJ;
	next->isCharge = info.isCharge;
	next->mass     = info.mass;
	next->epslon   = info.epslon;
	next->sigma    = info.sigma;
	next->charge   = info.charge;
	next->ident    = info.ident;
      }
    }

#ifdef IS_MPI
    
    void duplicate( atomStruct &info ){
      strcpy(info.name, name);
      info.isDirectional  = isDirectional;
      info.isLJ     = isLJ;
      info.isCharge = isCharge;
      info.mass     = mass;
      info.epslon   = epslon;
      info.sigma    = sigma;
      info.charge   = charge;
      info.ident    = ident;
      info.last     = 0;
    }

#endif

    char name[15];
    int isDirectional;
    int isLJ;
    int isCharge;
    double mass;
    double epslon;
    double sigma;
    double charge;
    int ident;
    LinkedAtomType* next;
  };

  class LinkedDirectionalType {
  public:
    LinkedDirectionalType(){ 
      next = NULL;
      name[0] = '\0';
    }
    ~LinkedDirectionalType(){ if( next != NULL ) delete next; }

    LinkedDirectionalType* find(char* key){ 
      if( !strcmp(name, key) ) return this;
      if( next != NULL ) return next->find(key);
      return NULL;
    }
    

    void add( directionalStruct &info ){
    
      // check for duplicates
      
      if( !strcmp( info.name, name ) ){
	sprintf( painCave.errMsg,
		 "Duplicate WATER directional type \"%s\" found in "
		 "the WATER param file./n",
		 name );
	painCave.isFatal = 1;
	simError();
      }
      
      if( next != NULL ) next->add(info);
      else{
	next = new LinkedDirectionalType();
	strcpy(next->name, info.name);
	next->isDipole = info.isDipole;
	next->isSticky = info.isSticky;
	next->Ixx      = info.Ixx;
	next->Iyy      = info.Iyy;
	next->Izz      = info.Izz;
	next->dipole   = info.dipole;
	next->w0       = info.w0;
	next->v0       = info.v0;
	next->v0p      = info.v0p;
	next->rl       = info.rl;
	next->ru       = info.ru;
	next->rlp      = info.rlp;
	next->rup      = info.rup;
      }
    }

#ifdef IS_MPI
    
    void duplicate( directionalStruct &info ){
      strcpy(info.name, name);
      info.isDipole = isDipole;
      info.isSticky = isSticky;
      info.Ixx      = Ixx;
      info.Iyy      = Iyy;
      info.Izz      = Izz;
      info.dipole   = dipole;
      info.w0       = w0;
      info.v0       = v0;
      info.v0p      = v0p;
      info.rl       = rl;
      info.ru       = ru;
      info.rlp      = rlp;
      info.rup      = rup;
      info.last     = 0;
    }

#endif

    char name[15];
    int isDipole;
    int isSticky;
    double Ixx;
    double Iyy;
    double Izz;
    double dipole;
    double w0;
    double v0;
    double v0p;
    double rl;
    double ru;
    double rlp;
    double rup;
    LinkedDirectionalType* next;
  };

  LinkedAtomType* headAtomType; 
  LinkedAtomType* currentAtomType;
  LinkedDirectionalType* headDirectionalType; 
  LinkedDirectionalType* currentDirectionalType;
} // namespace

using namespace WATER_NS;

//****************************************************************
// begins the actual forcefield stuff.	
//****************************************************************


WATER::WATER(){

  string fileName;
  string tempString;

  headAtomType = NULL;
  currentAtomType = NULL;
  headDirectionalType = NULL;
  currentDirectionalType = NULL;

#ifdef IS_MPI
  int i;
  
   // **********************************************************************
  // Init the atomStruct mpi type

  atomStruct atomProto; // mpiPrototype
  int atomBC[3] = {15,4,5};  // block counts
  MPI_Aint atomDspls[3];           // displacements
  MPI_Datatype atomMbrTypes[3];    // member mpi types

  MPI_Address(&atomProto.name, &atomDspls[0]);
  MPI_Address(&atomProto.mass, &atomDspls[1]);
  MPI_Address(&atomProto.isDirectional, &atomDspls[2]);
  
  atomMbrTypes[0] = MPI_CHAR;
  atomMbrTypes[1] = MPI_DOUBLE;
  atomMbrTypes[2] = MPI_INT;
  
  for (i=2; i >= 0; i--) atomDspls[i] -= atomDspls[0];
  
  MPI_Type_struct(3, atomBC, atomDspls, atomMbrTypes, &mpiAtomStructType);
  MPI_Type_commit(&mpiAtomStructType);

  // ***********************************************************************

   // **********************************************************************
  // Init the directionalStruct mpi type

  directionalStruct directionalProto; // mpiPrototype
  int directionalBC[3] = {15,11,3};  // block counts
  MPI_Aint directionalDspls[3];           // displacements
  MPI_Datatype directionalMbrTypes[3];    // member mpi types

  MPI_Address(&directionalProto.name, &directionalDspls[0]);
  MPI_Address(&directionalProto.Ixx, &directionalDspls[1]);
  MPI_Address(&directionalProto.isDipole, &directionalDspls[2]);
  
  directionalMbrTypes[0] = MPI_CHAR;
  directionalMbrTypes[1] = MPI_DOUBLE;
  directionalMbrTypes[2] = MPI_INT;
  
  for (i=2; i >= 0; i--) directionalDspls[i] -= directionalDspls[0];
  
  MPI_Type_struct(3, directionalBC, directionalDspls, directionalMbrTypes,
		  &mpiDirectionalStructType);
  MPI_Type_commit(&mpiDirectionalStructType);

  // ***********************************************************************

  if( worldRank == 0 ){
#endif
    
    // generate the force file name   

    fileName = "WATER.frc";

    //    fprintf( stderr,"Trying to open %s\n", fileName );
    
    // attempt to open the file in the current directory first.
    
    frcFile = fopen( fileName.c_str(), "r" );
    
    if( frcFile == NULL ){
      
      tempString = ffPath + "/" + fileName;
      fileName = tempString;
      
      frcFile = fopen( fileName.c_str(), "r" );
      
      if( frcFile == NULL ){
	
	sprintf( painCave.errMsg,
		 "Error opening the force field parameter file:\n"
                 "\t%s\n"
		 "\tHave you tried setting the FORCE_PARAM_PATH environment "
		 "variable?\n",
		 fileName.c_str() );
        painCave.severity = OOPSE_ERROR;
	painCave.isFatal = 1;
	simError();
      }
    }
    
#ifdef IS_MPI
  }
  
  sprintf( checkPointMsg, "WATER file opened sucessfully." );
  MPIcheckPoint();
  
#endif // is_mpi
}


WATER::~WATER(){

  if( headAtomType != NULL ) delete headAtomType;
  if( headDirectionalType != NULL ) delete headDirectionalType;

#ifdef IS_MPI
  if( worldRank == 0 ){
#endif // is_mpi
    
    fclose( frcFile );
    
#ifdef IS_MPI
  }
#endif // is_mpi
}

void WATER::cleanMe( void ){

#ifdef IS_MPI
  
  // keep the linked list in the mpi version

#else // is_mpi

  // delete the linked list in the single processor version

  if( headAtomType != NULL ) delete headAtomType;
  if( headDirectionalType != NULL ) delete headDirectionalType;

#endif // is_mpi
}


void WATER::initForceField(){
  
  initFortran(entry_plug->useReactionField );
}


void WATER::readParams( void ){

  int identNum, isError;
  int tempDirect0, tempDirect1;

  atomStruct atomInfo;
  directionalStruct directionalInfo;
  fpos_t *atomPos;
  AtomType* at;

  atomInfo.last = 1;         // initialize last to have the last set. 
  directionalInfo.last = 1;  // if things go well, last will be set to 0
  
  atomPos = new fpos_t;
  bigSigma = 0.0;
  
#ifdef IS_MPI
  if( worldRank == 0 ){
#endif
    
    // read in the atom types.

    headAtomType = new LinkedAtomType;
    headDirectionalType = new LinkedDirectionalType;
    
    fastForward( "AtomTypes", "initializeAtoms" );
    
    // we are now at the AtomTypes section.
    
    eof_test = fgets( readLine, sizeof(readLine), frcFile );
    lineNum++;
    
    
    // read a line, and start parsing out the atom types 
    
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
	  headAtomType->add( atomInfo );
	  if ( atomInfo.isDirectional ) {

	    // if the atom is directional, skip to the directional section
	    // and parse directional info.
	    fgetpos( frcFile, atomPos );
	    sectionSearch( "DirectionalTypes", atomInfo.name, 
			   "initializeDirectional" );
	    parseDirectional( readLine, lineNum, directionalInfo );
	    headDirectionalType->add( directionalInfo );

	    // return to the AtomTypes section
	    fsetpos( frcFile, atomPos ); 
	  }
	  identNum++;
	}
      }
      eof_test = fgets( readLine, sizeof(readLine), frcFile );
      lineNum++;
    }

#ifdef IS_MPI
  
    // send out the linked list to all the other processes
    
    sprintf( checkPointMsg,
	   "WATER atom and directional structures read successfully." );
    MPIcheckPoint();
    currentAtomType = headAtomType->next; //skip the first element place holder
    currentDirectionalType = headDirectionalType->next; // same w/ directional
    
    while( currentAtomType != NULL ){
      currentAtomType->duplicate( atomInfo );

      sendFrcStruct( &atomInfo, mpiAtomStructType );
      
      sprintf( checkPointMsg,
	       "successfully sent WATER force type: \"%s\"\n",
	       atomInfo.name );
      
      if ( atomInfo.isDirectional ){
	// send out the directional linked list to all the other processes
	
	currentDirectionalType->duplicate( directionalInfo );
	sendFrcStruct( &directionalInfo, mpiDirectionalStructType );
	
	sprintf( checkPointMsg,
		 "successfully sent WATER directional type: \"%s\"\n",
		 directionalInfo.name );
      }
      
      MPIcheckPoint();
      tempDirect0 = atomInfo.isDirectional;
      currentAtomType = currentAtomType->next;
      if( tempDirect0 ) 
	currentDirectionalType = currentDirectionalType->next;
    }
    
    atomInfo.last = 1;
    sendFrcStruct( &atomInfo, mpiAtomStructType );
    directionalInfo.last = 1;
    if ( atomInfo.isDirectional )
      sendFrcStruct( &directionalInfo, mpiDirectionalStructType ); 
  }
  
  else{
    // listen for node 0 to send out the force params
    
    MPIcheckPoint();

    headAtomType = new LinkedAtomType;
    headDirectionalType = new LinkedDirectionalType; 
    receiveFrcStruct( &atomInfo, mpiAtomStructType );

    if ( atomInfo.isDirectional )
      receiveFrcStruct( &directionalInfo, mpiDirectionalStructType );

    while( !atomInfo.last ){

      headAtomType->add( atomInfo );
      
      MPIcheckPoint();

      receiveFrcStruct( &atomInfo, mpiAtomStructType );

      if( atomInfo.isDirectional ){
	headDirectionalType->add( directionalInfo );

	receiveFrcStruct( &directionalInfo, mpiDirectionalStructType );
      }
    }
  }

#endif // is_mpi
  
  // dummy variables

  currentAtomType = headAtomType->next;
  currentDirectionalType = headDirectionalType->next;

  while( currentAtomType != NULL ){
    if( currentAtomType->name[0] != '\0' ){
      if (currentAtomType->isDirectional) 
	at = new DirectionalAtomType();         
      else 
	at = new AtomType();
      
      if (currentAtomType->isLJ) {
	at->setLennardJones();
      }

      if (currentAtomType->isCharge) {
	at->setCharge();
      }

      if (currentAtomType->isDirectional) {
	if (currentDirectionalType->isDipole) {
	  ((DirectionalAtomType*)at)->setDipole();
	}
      
	if (currentDirectionalType->isSticky) {
	  ((DirectionalAtomType*)at)->setSticky();
	}
      }
      
      at->setIdent(currentAtomType->ident);
      at->setName(currentAtomType->name);     
      at->complete();
    }
    currentAtomType = currentAtomType->next;
  }

  currentAtomType = headAtomType->next;
  currentDirectionalType = headDirectionalType->next;

  while( currentAtomType != NULL ){    

    currentDirectionalType = headDirectionalType->find(currentAtomType->name);
  
    if( currentAtomType->isLJ ){
      isError = 0;
      newLJtype( &(currentAtomType->ident), &(currentAtomType->sigma), 
                 &(currentAtomType->epslon), &isError);
      if( isError ){
        sprintf( painCave.errMsg,
                 "Error initializing the \"%s\" LJ type in fortran\n",
                 currentAtomType->name );
        painCave.isFatal = 1;
        simError();
      }
    }

    if (currentAtomType->isCharge) {
      newChargeType(&(currentAtomType->ident), &(currentAtomType->charge),
		    &isError);
      if( isError ){
	sprintf( painCave.errMsg,
		 "Error initializing the \"%s\" charge type in fortran\n",
		 currentAtomType->name );
	painCave.isFatal = 1;
	simError();
      }
    }

    if (currentAtomType->isDirectional){
      if (currentDirectionalType->isDipole) {
	newDipoleType(&(currentAtomType->ident), 
		      &(currentDirectionalType->dipole),
		      &isError);
	if( isError ){
	  sprintf( painCave.errMsg,
		   "Error initializing the \"%s\" dipole type in fortran\n",
		   currentDirectionalType->name );
	  painCave.isFatal = 1;
	  simError();
	}
      }
      
      if(currentDirectionalType->isSticky) {        
	makeStickyType( &(currentDirectionalType->w0), 
			&(currentDirectionalType->v0), 
			&(currentDirectionalType->v0p), 
			&(currentDirectionalType->rl), 
			&(currentDirectionalType->ru), 
			&(currentDirectionalType->rlp), 
			&(currentDirectionalType->rup));
      }
    }
    currentAtomType = currentAtomType->next;
  }
  
#ifdef IS_MPI
  sprintf( checkPointMsg,
	   "WATER atom and directional structures successfully"
	   "sent to fortran\n" );
  MPIcheckPoint();
#endif // is_mpi
  
}


void WATER::initializeAtoms( int nAtoms, Atom** the_atoms ){
  
  int i,j,k;

  // initialize the atoms
  DirectionalAtom* dAtom;
  double ji[3];
  double inertialMat[3][3];

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

    if (currentAtomType->isLJ) entry_plug->useLennardJones = 1;
    if (currentAtomType->isCharge) entry_plug->useCharges = 1;
    if (currentAtomType->isDirectional) {
      if (currentDirectionalType->isDipole) entry_plug->useDipoles = 1;      
      if (currentDirectionalType->isSticky) entry_plug->useSticky = 1;
    }

    if( bigSigma < currentAtomType->sigma ) bigSigma = currentAtomType->sigma;

    the_atoms[i]->setHasCharge(currentAtomType->isCharge);

    if( currentAtomType->isDirectional ){
      currentDirectionalType = 
	headDirectionalType->find( the_atoms[i]->getType() );
      if( currentDirectionalType == NULL ){
	sprintf( painCave.errMsg, 
		 "DirectionalType error, %s not found in force file.\n",
		 the_atoms[i]->getType() );
	painCave.isFatal = 1;
	simError();
      }

      // zero out the moments of inertia matrix
      for( j=0; j<3; j++ )
	for( k=0; k<3; k++ )
	  inertialMat[j][k] = 0.0;

      // load the force file moments of inertia
      inertialMat[0][0] = currentDirectionalType->Ixx;
      inertialMat[1][1] = currentDirectionalType->Iyy;
      inertialMat[2][2] = currentDirectionalType->Izz;
 
      dAtom = (DirectionalAtom *) the_atoms[i];
      dAtom->setHasDipole( currentDirectionalType->isDipole );

      ji[0] = 0.0;
      ji[1] = 0.0;
      ji[2] = 0.0;
      dAtom->setJ( ji );
      dAtom->setI( inertialMat );
 
      entry_plug->n_dipoles++;
    }
  }
}

void WATER::initializeBonds( int nBonds, Bond** BondArray,
			     bond_pair* the_bonds ){
  
    if( nBonds ){
      sprintf( painCave.errMsg,
	       "WATER does not support bonds.\n" );
      painCave.isFatal = 1;
      simError();
    }
}

void WATER::initializeBends( int nBends, Bend** bendArray,
			     bend_set* the_bends ){

    if( nBends ){
      sprintf( painCave.errMsg,
	       "WATER does not support bends.\n" );
      painCave.isFatal = 1;
      simError();
    }
}

void WATER::initializeTorsions( int nTorsions, Torsion** torsionArray,
				torsion_set* the_torsions ){

    if( nTorsions ){
      sprintf( painCave.errMsg,
	       "WATER does not support torsions.\n" );
      painCave.isFatal = 1;
      simError();
    }
}

void WATER::fastForward( char* stopText, char* searchOwner ){

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

void WATER::sectionSearch( char* secHead, char* stopText, char* searchOwner ){

  int foundSection = 0;
  int foundText = 0;
  char* the_token;
  fpos_t *tempPos;

  rewind( frcFile );
  lineNum = 0;
  tempPos = new fpos_t;

  eof_test = fgets( readLine, sizeof(readLine), frcFile );
  lineNum++;
  if( eof_test == NULL ){
    sprintf( painCave.errMsg, "Error fast forwarding force file for %s: "
	     " file is empty.\n",
	     searchOwner );
    painCave.isFatal = 1;
    simError();
  }
  
  while( !foundSection ){
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
    foundSection = !strcmp( secHead, the_token );
    
    if( !foundSection ){
      eof_test = fgets( readLine, sizeof(readLine), frcFile );
      lineNum++;
      
      if( eof_test == NULL ){
	sprintf( painCave.errMsg, 
		 "Error section searching force file for %s at "
		 "line %d: file ended unexpectedly.\n",
		 searchOwner,
		 lineNum );
	painCave.isFatal = 1;
	simError();
      } 
    }
  }

  while( !foundText ){
    
    fgetpos( frcFile, tempPos );
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
    
    the_token = strtok( readLine, " ,;\t#\n" );
    if( the_token != NULL ){
      foundText = !strcmp( stopText, the_token );
    }
  }  
  fsetpos( frcFile, tempPos );
  eof_test = fgets( readLine, sizeof(readLine), frcFile );
  lineNum++;
}


int WATER_NS::parseAtom( char *lineBuffer, int lineNum, atomStruct &info ){

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
    
    info.isDirectional = atoi( the_token );
    
    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing AtomTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }

    info.isLJ = atoi( the_token );
    
    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing AtomTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }

    info.isCharge = atoi( the_token );
    
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

    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing AtomTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }

    info.charge = atof( the_token );

    return 1;
  }
  else return 0;
}

int WATER_NS::parseDirectional( char *lineBuffer, int lineNum, directionalStruct &info ){

  char* the_token;
  
  the_token = strtok( lineBuffer, " \n\t,;" );
  if( the_token != NULL ){
    
    strcpy( info.name, the_token );
	  
    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing DirectionalTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }
    
    info.isDipole = atoi( the_token );
    
    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing DirectionalTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }

    info.isSticky = atoi( the_token );
    
    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing DirectionalTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }

    info.Ixx = atof( the_token );
    
    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing DirectionalTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }

    info.Iyy = atof( the_token );
    
    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing DirectionalTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }

    info.Izz = atof( the_token );
    
    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing DirectionalTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }


    info.dipole = atof( the_token );
    
    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing DirectionalTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }
	
    info.w0 = atof( the_token );
	  
    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing DirectionalTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }
	
    info.v0 = atof( the_token );

    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing DirectionalTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }

    info.v0p = atof( the_token );

    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing DirectionalTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }

    info.rl = atof( the_token );

    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing DirectionalTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }

    info.ru = atof( the_token );

    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing DirectionalTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }

    info.rlp = atof( the_token );

    if( ( the_token = strtok( NULL, " \n\t,;" ) ) == NULL ){
      sprintf( painCave.errMsg, 
	       "Error parseing DirectionalTypes: line %d\n", lineNum );
      painCave.isFatal = 1;
      simError();
    }

    info.rup = atof( the_token );

    return 1;
  }
  else return 0;
}
