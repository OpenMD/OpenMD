#ifdef IS_MPI
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define __mpiBASSEVENT
#include "mpiBASS.h"
#include "simError.h"


void mpiCatchEvent(void);




void mpiEventInit(void)
{
  int blockCounts[5] = {1,3,1,120,80};
  MPI_Aint dspls[5];
  MPI_Datatype types[5];
  mBEvent protoEvent;

  int i;  

  MPI_Address(&protoEvent.type, &dspls[0]);
  MPI_Address(&protoEvent.d1, &dspls[1]);
  MPI_Address(&protoEvent.i1, &dspls[2]);
  MPI_Address(&protoEvent.cArray, &dspls[3]);
  MPI_Address(&protoEvent.lhs   , &dspls[4]);


  types[0] = MPI_INT;
  types[1] = MPI_DOUBLE;
  types[2] = MPI_INT;
  types[3] = MPI_CHAR;
  types[4] = MPI_CHAR;


  for (i =4; i >= 0; i--) dspls[i] -= dspls[0];

  MPI_Type_struct(5,blockCounts,dspls,types,&mpiBASSEventType);
  MPI_Type_commit(&mpiBASSEventType);
}


void throwMPIEvent(event* the_event)
{
  mBEvent mpiEventContainer;
  int mpiStatus;


  if (the_event == NULL)     mpiStatus = MPI_INTERFACE_DONE;
  else                       mpiStatus = MPI_INTERFACE_CONTINUE;
  
  MPI_Bcast(&mpiStatus,1,MPI_INT,0,MPI_COMM_WORLD);
  
  if (!mpiStatus){
    switch (the_event->event_type){
    case MOLECULE:
      mpiEventContainer.type = mpiMOLECULE;
      mpiEventContainer.i1 = the_event->evt.blk_index; // pack block index into first int  
    break;

    case RIGIDBODY:
      mpiEventContainer.type = mpiRIGIDBODY;
      mpiEventContainer.i1 = the_event->evt.blk_index; // pack block index into first int  
    break;

    case CUTOFFGROUP:
      mpiEventContainer.type = mpiCUTOFFGROUP;
      mpiEventContainer.i1 = the_event->evt.blk_index; // pack block index into first int  
    break;

    case ATOM:
      mpiEventContainer.type = mpiATOM;
      mpiEventContainer.i1 = the_event->evt.blk_index; // pack block index into first int
      break;

    case BOND:
      mpiEventContainer.type = mpiBOND;
      mpiEventContainer.i1 = the_event->evt.blk_index; // pack block index into first int
      break;

    case BEND:
      mpiEventContainer.type = mpiBEND;
      mpiEventContainer.i1 = the_event->evt.blk_index; // pack block index into first int
      break;

    case TORSION:
      mpiEventContainer.type = mpiTORSION;
      mpiEventContainer.i1 = the_event->evt.blk_index; // pack block index into first int
      break;

    case ZCONSTRAINT:
      mpiEventContainer.type = mpiZCONSTRAINT;
      mpiEventContainer.i1 = the_event->evt.blk_index; // pack block index into first int
      break;

    case COMPONENT:
      mpiEventContainer.type = mpiCOMPONENT;
      mpiEventContainer.i1 = the_event->evt.blk_index; // pack block index into first int
      break;

    case POSITION:
      mpiEventContainer.type = mpiPOSITION;
      mpiEventContainer.d1 = the_event->evt.pos.x; // pack pos coord into d
      mpiEventContainer.d2 = the_event->evt.pos.y; 
      mpiEventContainer.d3 = the_event->evt.pos.z;   
      break;

    case ORIENTATION:
      mpiEventContainer.type = mpiORIENTATION;
      mpiEventContainer.d1 = the_event->evt.ornt.phi; // pack orientation coord into d
      mpiEventContainer.d2 = the_event->evt.ornt.theta; 
      mpiEventContainer.d3 = the_event->evt.ornt.psi;   
      break;
      
    case CONSTRAINT:
      mpiEventContainer.type = mpiCONSTRAINT;
      mpiEventContainer.d1 = the_event->evt.cnstr; // pack constraint coord into d
      break;
      
    case MEMBERS:
      mpiEventContainer.type = mpiMEMBERS;
      mpiEventContainer.i1 = the_event->evt.mbrs.nMembers ; // pack member ints into i
      break;
      
    case ASSIGNMENT:

      strcpy(mpiEventContainer.lhs,the_event->evt.asmt.lhs);

#ifdef MPIBASS_VERBOSE

      fprintf(stderr, 
	      "mpiDiag at node %d: evt Assignment: \"%s\" = ?\n"
	      "                    mpi Assignment: \"%s\" = ?\n",
	      worldRank, 
	      the_event->evt.asmt.lhs,
	      mpiEventContainer.lhs);
#endif

      switch (the_event->evt.asmt.asmt_type){
      case STRING:
	mpiEventContainer.type = mpiASSIGNMENT_s;
	strcpy(mpiEventContainer.cArray,the_event->evt.asmt.rhs.sval);
	break;

      case INT:
	mpiEventContainer.type = mpiASSIGNMENT_i;
	mpiEventContainer.i1 = the_event->evt.asmt.rhs.ival;
	break;

      case DOUBLE:
	mpiEventContainer.type = mpiASSIGNMENT_d;
	mpiEventContainer.d1 = the_event->evt.asmt.rhs.dval;
	break;
      }
      break;

    case BLOCK_END:
      mpiEventContainer.type = mpiBLOCK_END;
      break;
    }


    MPI_Bcast(&mpiEventContainer,1,mpiBASSEventType,0,MPI_COMM_WORLD);

    if (the_event->event_type == MEMBERS) {

      // For member lists, we need a separate broadcast to spew out the
      // membership array:
      MPI_Bcast(the_event->evt.mbrs.memberList, the_event->evt.mbrs.nMembers, 
                MPI_INT, 0, MPI_COMM_WORLD);
      
    }   

    sprintf( checkPointMsg, 
	     "BASS Event broadcast successful" );

    MPIcheckPoint();
  }
}


// Everybody but node 0 runs this
void mpiEventLoop(void)
{
  int mpiContinue;

#ifdef MPIBASS_VERBOSE
  fprintf(stderr, 
	  "event key List at node %d:\n"
	  "  MOLECULE    %d\n"
	  "  ATOM        %d\n"
	  "  BOND        %d\n"
	  "  BEND        %d\n"
	  "  TORSION     %d\n"
	  "  COMPONENT   %d\n"
	  "  POSITION    %d\n"
	  "  ASSIGNMENT  %d\n"
	  "  MEMBERS     %d\n"
	  "  CONSTRAINT  %d\n"
	  "  ORIENTATION %d\n"
	  "  ZCONSTRAINT %d\n"
	  "  RIGIDBODY   %d\n"
	  "  CUTOFFGROUP %d\n"
	  "  BLOCK_END   %d\n"
	  "\n",
	  worldRank,
	  MOLECULE, ATOM, BOND, BEND, TORSION, COMPONENT, 
	  POSITION, ASSIGNMENT, MEMBERS, CONSTRAINT, ORIENTATION,
	  ZCONSTRAINT, RIGIDBODY, CUTOFFGROUP, BLOCK_END );
#endif

  MPI_Bcast(&mpiContinue,1,MPI_INT,0,MPI_COMM_WORLD);

  while(!mpiContinue){
    
    mpiCatchEvent();

    MPI_Bcast(&mpiContinue,1,MPI_INT,0,MPI_COMM_WORLD);
  }

  if (mpiContinue == MPI_INTERFACE_ABORT){
    MPI_Finalize();
    exit (0);
  }
}

void mpiCatchEvent(void)
{
  event the_event;
  mBEvent mpiEventContainer;

  
  MPI_Bcast(&mpiEventContainer,1,mpiBASSEventType,0,MPI_COMM_WORLD);

  switch (mpiEventContainer.type){
  case mpiMOLECULE:
    the_event.event_type = MOLECULE;
    the_event.evt.blk_index = mpiEventContainer.i1;
    break;

  case mpiRIGIDBODY:
    the_event.event_type = RIGIDBODY;
    the_event.evt.blk_index = mpiEventContainer.i1;
    break;

  case mpiCUTOFFGROUP:
    the_event.event_type = CUTOFFGROUP;
    the_event.evt.blk_index = mpiEventContainer.i1;
    break;
    
  case mpiATOM:
    the_event.event_type = ATOM;
    the_event.evt.blk_index = mpiEventContainer.i1;
    break;
    
  case mpiBOND:
    the_event.event_type = BOND;
    the_event.evt.blk_index = mpiEventContainer.i1;
    break;
      
  case mpiBEND:
    the_event.event_type = BEND;
    the_event.evt.blk_index = mpiEventContainer.i1;
    break;

  case mpiTORSION:
    the_event.event_type = TORSION;
    the_event.evt.blk_index = mpiEventContainer.i1;
    break;

  case mpiZCONSTRAINT:
    the_event.event_type = ZCONSTRAINT;
    the_event.evt.blk_index = mpiEventContainer.i1;
    break;
    
  case mpiCOMPONENT:
    the_event.event_type = COMPONENT;
    the_event.evt.blk_index = mpiEventContainer.i1;
    break;


  case mpiPOSITION:
    the_event.event_type = POSITION;
    the_event.evt.pos.x = mpiEventContainer.d1; 
    the_event.evt.pos.y = mpiEventContainer.d2; 
    the_event.evt.pos.z = mpiEventContainer.d3;   
    break;

  case mpiORIENTATION:
    the_event.event_type = ORIENTATION;
    the_event.evt.ornt.phi   = mpiEventContainer.d1; 
    the_event.evt.ornt.theta = mpiEventContainer.d2; 
    the_event.evt.ornt.psi   = mpiEventContainer.d3;   
    break;
      
  case mpiCONSTRAINT:
    the_event.event_type = CONSTRAINT;
    the_event.evt.cnstr = mpiEventContainer.d1; 
    break;
    
  case mpiMEMBERS:
    the_event.event_type = MEMBERS;
    the_event.evt.mbrs.nMembers = mpiEventContainer.i1; 

    the_event.evt.mbrs.memberList = (int *) calloc(the_event.evt.mbrs.nMembers,
                                                    sizeof(int));
    
    // Grab the member list since we have a number of members:
    MPI_Bcast(the_event.evt.mbrs.memberList, the_event.evt.mbrs.nMembers, 
              MPI_INT, 0, MPI_COMM_WORLD);
    
    break;
    
  case mpiASSIGNMENT_s:
    the_event.event_type = ASSIGNMENT;
    the_event.evt.asmt.asmt_type = STRING;
    strcpy(the_event.evt.asmt.lhs,mpiEventContainer.lhs);
    strcpy(the_event.evt.asmt.rhs.sval,mpiEventContainer.cArray);
    
#ifdef MPIBASS_VERBOSE

    fprintf(stderr, "mpiDiag at node %d: Assignment:  %s = %s\n", worldRank, 
	    the_event.evt.asmt.lhs,
	    the_event.evt.asmt.rhs.sval );
#endif

    break;

  case mpiASSIGNMENT_i:
    the_event.event_type = ASSIGNMENT;
    the_event.evt.asmt.asmt_type = INT;
    strcpy(the_event.evt.asmt.lhs,mpiEventContainer.lhs);
    the_event.evt.asmt.rhs.ival = mpiEventContainer.i1;      

#ifdef MPIBASS_VERBOSE

    fprintf(stderr, "mpiDiag at node %d: Assignment:  %s = %d\n", worldRank, 
	    the_event.evt.asmt.lhs,
	    the_event.evt.asmt.rhs.ival );
#endif

    break;

  case mpiASSIGNMENT_d:
    the_event.event_type = ASSIGNMENT;
    the_event.evt.asmt.asmt_type = DOUBLE;
    strcpy(the_event.evt.asmt.lhs,mpiEventContainer.lhs);
    the_event.evt.asmt.rhs.dval = mpiEventContainer.d1;      

#ifdef MPIBASS_VERBOSE

    fprintf(stderr, "mpiDiag at node %d: Assignment:  %s = %lf\n", worldRank, 
	    the_event.evt.asmt.lhs,
	    the_event.evt.asmt.rhs.dval );
#endif

    break;

  case mpiBLOCK_END:
    the_event.event_type = BLOCK_END;
    break;
  }

#ifdef MPIBASS_VERBOSE

  fprintf(stderr, "mpiDiag at node %d: event type is %d\n", worldRank, 
	  the_event.event_type);
#endif

  if (!event_handler(&the_event)){ 

    sprintf(painCave.errMsg,
	    "MPI event handling error at node %d => %s\n",
	    worldRank,
	    the_event.err_msg);
    painCave.isFatal = 1;
    simError();
  }
  MPIcheckPoint();
}


void mpiInterfaceExit(void){
  int mpiStatus = MPI_INTERFACE_ABORT;
  
  MPI_Bcast(&mpiStatus,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Finalize();
  exit (0);

}
#endif
