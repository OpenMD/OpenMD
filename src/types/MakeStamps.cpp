#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "MakeStamps.hpp"
#include "MoleculeStamp.hpp"
#include "RigidBodyStamp.hpp"
#include "CutoffGroupStamp.hpp"
#include "simError.h"
#ifdef IS_MPI
#include "mpiBASS.h"
#endif // is_mpi

LinkedMolStamp::~LinkedMolStamp(){ 
  if( mol_stamp != NULL ) delete mol_stamp;
  if( next != NULL ) delete next; 
}

void LinkedMolStamp::add( LinkedMolStamp* newbie ){
  
  if( next != NULL ) next->add( newbie );
  else{
    next = newbie;
    next->setPrev( this );
  }
}

MoleculeStamp* LinkedMolStamp::match( char* id ){

  if( mol_stamp != NULL ){
    if(!strcmp( mol_stamp->getID(), id )){
      
      // make sure we aren't hiding somebody else with the same name
      if(next != NULL ){
	if( next->match( id ) != NULL){
	  sprintf( painCave.errMsg,
		   "Molecule Stamp Error. Two separate of declarations of "
		   "%s present.\n",
		   id );
	  painCave.isFatal = 1;
	  simError();
#ifdef IS_MPI
	  if( painCave.isEventLoop ){
	    if( worldRank == 0 ) mpiInterfaceExit();
	  }
#endif //is_mpi
	}
	else return mol_stamp;
      }
      else  return mol_stamp;
    }
  }
  
  if( next != NULL ) return next->match( id );
  
  return NULL;
}

LinkedMolStamp* LinkedMolStamp::extract( char* id ){

  if( mol_stamp != NULL ){
    if(!strcmp( mol_stamp->getID(), id )){
      
      // make sure we aren't hiding somebody else with the same name
      if(next != NULL ){
	if( next->match( id ) != NULL){
	  sprintf( painCave.errMsg,
		   "Molecule Stamp Error. Two separate of declarations of "
		   "%s present.\n",
		   id );
	  painCave.isFatal = 1;
	  simError();
#ifdef IS_MPI
	  if( painCave.isEventLoop ){
	    if( worldRank == 0 ) mpiInterfaceExit();
	  }
#endif //is_mpi
	}
      }
	  
      prev->setNext( next );
      if( next != NULL ) next->setPrev( prev );
      return this;
    }
  }
  
  if( next != NULL ) return next->extract( id );
  
  return NULL;
}

MakeStamps::MakeStamps(){
  
  int i;

  hash_size = 47;
  hash_shift = 4;
  
  my_mols = new LinkedMolStamp*[hash_size];
  for( i=0; i<hash_size; i++ ){
    my_mols[i] = new LinkedMolStamp();
  }

}

MakeStamps::~MakeStamps(){
  
  int i;

  for( i=0; i<hash_size; i++ ){
    if( my_mols[i] != NULL ) delete my_mols[i];
  }
  delete[] my_mols;
}

int MakeStamps::hash( char* text ){

  register unsigned short int i = 0; // loop counter
  int key = 0; // the hash key

  while( text[i] != '\0' ){

    key = ( ( key << hash_shift ) + text[i] ) % hash_size;
    
    i++;
  }
  
  if( key < 0 ){

    // if the key is less than zero, we've had an overflow error

    sprintf( painCave.errMsg,
	     "There has been an overflow error in the MakeStamps hash key.");
    painCave.isFatal = 1;
    simError();
#ifdef IS_MPI
    if( painCave.isEventLoop ){
      if( worldRank == 0 ) mpiInterfaceExit();
    }
#endif //is_mpi
  }
  
  return key;
}

LinkedMolStamp* MakeStamps::extractMolStamp( char* the_id ){
  int key;
  
  key = hash( the_id );
  if( my_mols[key] != NULL ){
    return my_mols[key]->extract( the_id );
  }
  
  return NULL;
}

void MakeStamps::addMolStamp( MoleculeStamp* the_stamp ){

  int key;
  LinkedMolStamp* linked_mol;
  
  key = hash( the_stamp->getID() );

  linked_mol = new LinkedMolStamp;
  linked_mol->setStamp( the_stamp );

  my_mols[key]->add( linked_mol );

}


int MakeStamps::newMolecule( event* the_event ){
  
  current_mol = new MoleculeStamp;
  return 1;
}

int MakeStamps::moleculeAssign( event* the_event ){

  switch( the_event->evt.asmt.asmt_type ){
    
  case STRING:
    the_event->err_msg = current_mol->assignString( the_event->evt.asmt.lhs,
				      the_event->evt.asmt.rhs.sval );
    break;

  case DOUBLE:
    the_event->err_msg = current_mol->assignDouble( the_event->evt.asmt.lhs,
				      the_event->evt.asmt.rhs.dval );
    break;
    
  case INT:
    the_event->err_msg = current_mol->assignInt( the_event->evt.asmt.lhs,
				   the_event->evt.asmt.rhs.ival );
    break;

  default:
    the_event->err_msg = strdup( "MakeStamp error. Invalid molecule"
				 " assignment type" );
    return 0;
    break;
  }
  if( the_event->err_msg != NULL ) return 0;
  return 1;
}

int MakeStamps::moleculeEnd( event* the_event ){

  the_event->err_msg = current_mol->checkMe();
  // if err_msg is set, then something is wrong
  if( the_event->err_msg != NULL ) return 0; 

  addMolStamp( current_mol );
  return 1;
}

int MakeStamps::newRigidBody( event* the_event ){
  
  current_rigidbody = new RigidBodyStamp;

  the_event->err_msg = current_mol->addRigidBody( current_rigidbody,
                                                  the_event->evt.blk_index );
  if( the_event->err_msg != NULL ) return 0;
  return 1;
}

int MakeStamps::rigidBodyAssign( event* the_event ){

  switch( the_event->evt.asmt.asmt_type ){
    
  case STRING:
    the_event->err_msg = 
      current_rigidbody->assignString( the_event->evt.asmt.lhs,
                                       the_event->evt.asmt.rhs.sval );
    if( the_event->err_msg != NULL ) return 0;
    return 1;
    break;
    
  case DOUBLE:
    the_event->err_msg = 
      current_rigidbody->assignDouble( the_event->evt.asmt.lhs,
                                       the_event->evt.asmt.rhs.dval );
    if( the_event->err_msg != NULL ) return 0;
    return 1;    
    break;

  case INT:
    the_event->err_msg = 
      current_rigidbody->assignInt( the_event->evt.asmt.lhs,
                                    the_event->evt.asmt.rhs.ival );
    if( the_event->err_msg != NULL ) return 0;
    return 1;
    break;
    
  default:
    the_event->err_msg = strdup( "MakeStamp error. Invalid rigidBody"
				 " assignment type" );
    return 0;
    break;
  }
  return 0;
}

int MakeStamps::rigidBodyMembers( event* the_event ){

  int i;

  if( the_event->evt.mbrs.nMembers > 0 ){

    for (i = 0; i < the_event->evt.mbrs.nMembers; i++) {      
      current_rigidbody->addMember(the_event->evt.mbrs.memberList[i]);
    }
    
    return 1;
    
  } else {
    the_event->err_msg = strdup( "MakeStamp error. No members in memberList "
                                 " for this rigidBody.");
    return 0;

  }
}

int MakeStamps::rigidBodyEnd( event* the_event ){

  the_event->err_msg = current_rigidbody->checkMe();
  if( the_event->err_msg != NULL ) return 0;
  
  return 1;
}

int MakeStamps::newCutoffGroup( event* the_event ){
  
  current_cutoffgroup = new CutoffGroupStamp;
  
  the_event->err_msg = current_mol->addCutoffGroup( current_cutoffgroup,
                                                    the_event->evt.blk_index );
  if( the_event->err_msg != NULL ) return 0;
  return 1;
}

int MakeStamps::cutoffGroupAssign( event* the_event ){
  
  switch( the_event->evt.asmt.asmt_type ){
    
  case STRING:
    the_event->err_msg = 
      current_cutoffgroup->assignString( the_event->evt.asmt.lhs,
                                         the_event->evt.asmt.rhs.sval );
    if( the_event->err_msg != NULL ) return 0;
    return 1;
    break;
    
  case DOUBLE:
    the_event->err_msg = 
      current_cutoffgroup->assignDouble( the_event->evt.asmt.lhs,
                                         the_event->evt.asmt.rhs.dval );
    if( the_event->err_msg != NULL ) return 0;
    return 1;    
    break;
    
  case INT:
    the_event->err_msg = 
      current_cutoffgroup->assignInt( the_event->evt.asmt.lhs,
                                      the_event->evt.asmt.rhs.ival );
    if( the_event->err_msg != NULL ) return 0;
    return 1;
    break;
    
  default:
    the_event->err_msg = strdup( "MakeStamp error. Invalid CutoffGroup"
				 " assignment type" );
    return 0;
    break;
  }
  return 0;
}

int MakeStamps::cutoffGroupMembers( event* the_event ){

  int i;
  
  if( the_event->evt.mbrs.nMembers > 0 ){
    
    for (i = 0; i < the_event->evt.mbrs.nMembers; i++) {      
      current_cutoffgroup->addMember(the_event->evt.mbrs.memberList[i]);
    }
    
    return 1;
    
  } else {
    the_event->err_msg = strdup( "MakeStamp error. No members in memberList "
                                 " for this CutoffGroup.");
    return 0;

  }
}

int MakeStamps::cutoffGroupEnd( event* the_event ){
  
  the_event->err_msg = current_cutoffgroup->checkMe();
  if( the_event->err_msg != NULL ) return 0;
  
  return 1;
}

int MakeStamps::newAtom( event* the_event ){
  
  current_atom = new AtomStamp;
  
  the_event->err_msg = current_mol->addAtom( current_atom,
					     the_event->evt.blk_index );
  
  if( the_event->err_msg != NULL ) return 0;
  return 1;
}

int MakeStamps::atomPosition( event* the_event ){
  
  current_atom->setPosition( the_event->evt.pos.x,
			     the_event->evt.pos.y,
			     the_event->evt.pos.z );
  return 1;
}


int MakeStamps::atomOrientation( event* the_event ){
  
  current_atom->setOrientation( the_event->evt.ornt.phi,
				the_event->evt.ornt.theta,
				the_event->evt.ornt.psi );
  return 1;
}

int MakeStamps::atomAssign( event* the_event ){

  switch( the_event->evt.asmt.asmt_type ){
    
  case STRING:
    the_event->err_msg = 
      current_atom->assignString( the_event->evt.asmt.lhs,
				  the_event->evt.asmt.rhs.sval );
    if( the_event->err_msg != NULL ) return 0;
    return 1;
    break;

  case DOUBLE:
    the_event->err_msg = 
      current_atom->assignDouble( the_event->evt.asmt.lhs,
				  the_event->evt.asmt.rhs.dval );
    if( the_event->err_msg != NULL ) return 0;
    return 1;    
    break;

  case INT:
    the_event->err_msg = 
      current_atom->assignInt( the_event->evt.asmt.lhs,
			       the_event->evt.asmt.rhs.ival );
    if( the_event->err_msg != NULL ) return 0;
    return 1;
    break;

  default:
    the_event->err_msg = strdup( "MakeStamp error. Invalid atom"
				 " assignment type" );
    return 0;
    break;
  }
  return 0;
}

int MakeStamps::atomEnd( event* the_event ){

  the_event->err_msg = current_atom->checkMe();
  if( the_event->err_msg != NULL ) return 0;
  
  return 1;
}

int MakeStamps::newBond( event* the_event ){
  
  current_bond = new BondStamp;
  
  the_event->err_msg = current_mol->addBond( current_bond,
					     the_event->evt.blk_index );
  if( the_event->err_msg != NULL ) return 0;

  return 1;
}

int MakeStamps::bondAssign( event* the_event ){

  switch( the_event->evt.asmt.asmt_type ){
    
  case STRING:
    current_bond->assignString( the_event->evt.asmt.lhs,
			       the_event->evt.asmt.rhs.sval );
    return 1;
    break;

  case DOUBLE:
    current_bond->assignDouble( the_event->evt.asmt.lhs,
			       the_event->evt.asmt.rhs.dval );
    return 1;
    break;

  case INT:
    current_bond->assignInt( the_event->evt.asmt.lhs,
			    the_event->evt.asmt.rhs.ival );
    return 1;
    break;

  default:
    the_event->err_msg = strdup( "MakeStamp error. Invalid bond"
				 " assignment type" );
    return 0;
    break;
  }
  return 0;
}

int MakeStamps::bondMembers( event* the_event ){
  
  if( the_event->evt.mbrs.nMembers == 2 ){
    
    current_bond->members( the_event->evt.mbrs.memberList[0],
                           the_event->evt.mbrs.memberList[1] );
    return 1;

  } else {
    the_event->err_msg = strdup( "MakeStamp error. Wrong number of members "
                                 " in bond");
    return 0;

  }

}

int MakeStamps::bondConstraint( event* the_event ){

  current_bond->constrain( the_event->evt.cnstr );
  return 1;
}

int MakeStamps::bondEnd( event* the_event ){

  the_event->err_msg = current_bond->checkMe();
  if( the_event->err_msg != NULL ) return 0;
  
  return 1;
}

int MakeStamps::newBend( event* the_event ){
  
  current_bend = new BendStamp;
  
  the_event->err_msg = current_mol->addBend( current_bend,
					     the_event->evt.blk_index );
  if( the_event->err_msg != NULL ) return 0;

  return 1;
}

int MakeStamps::bendAssign( event* the_event ){

  switch( the_event->evt.asmt.asmt_type ){
    
  case STRING:
    current_bend->assignString( the_event->evt.asmt.lhs,
			       the_event->evt.asmt.rhs.sval );
    return 1;
    break;

  case DOUBLE:
    current_bend->assignDouble( the_event->evt.asmt.lhs,
			       the_event->evt.asmt.rhs.dval );
    return 1;
    break;

  case INT:
    current_bend->assignInt( the_event->evt.asmt.lhs,
			    the_event->evt.asmt.rhs.ival );
    return 1;
    break;

  default:
    the_event->err_msg = strdup( "MakeStamp error. Invalid bend"
				 " assignment type" );
    return 0;
    break;
  }
  return 0;
}

int MakeStamps::bendMembers( event* the_event ){


  switch( the_event->evt.mbrs.nMembers ) {
  case 3:
    current_bend->members( the_event->evt.mbrs.memberList[0],
                           the_event->evt.mbrs.memberList[1],
                           the_event->evt.mbrs.memberList[2]);
    return 1;
    break;
  case 2:
    current_bend->members( the_event->evt.mbrs.memberList[0],
                           the_event->evt.mbrs.memberList[1],
                           0 );
    return 1;
    break;
  default: 	 
    the_event->err_msg = strdup( "MakeStamp error. Wrong number of members "
                                 "in bend.");
    return 0;
    break;
  }
  return 0;
}

int MakeStamps::bendConstraint( event* the_event ){

  current_bend->constrain( the_event->evt.cnstr );
  return 1;
}

int MakeStamps::bendEnd( event* the_event ){

  the_event->err_msg = current_bend->checkMe();
  if( the_event->err_msg != NULL ) return 0;
  
  return 1;
}

int MakeStamps::newTorsion( event* the_event ){
  
  current_torsion = new TorsionStamp;
  
  the_event->err_msg = current_mol->addTorsion( current_torsion,
						the_event->evt.blk_index );
  if( the_event->err_msg != NULL ) return 0;

  return 1;
}

int MakeStamps::torsionAssign( event* the_event ){

  switch( the_event->evt.asmt.asmt_type ){
    
  case STRING:
    current_torsion->assignString( the_event->evt.asmt.lhs,
			       the_event->evt.asmt.rhs.sval );
    return 1;
    break;

  case DOUBLE:
    current_torsion->assignDouble( the_event->evt.asmt.lhs,
			       the_event->evt.asmt.rhs.dval );
    return 1;
    break;

  case INT:
    current_torsion->assignInt( the_event->evt.asmt.lhs,
			    the_event->evt.asmt.rhs.ival );
    return 1;
    break;

  default:
    the_event->err_msg = strdup( "MakeStamp error. Invalid torsion"
				 " assignment type" );
    return 0;
    break;
  }
  return 0;
}

int MakeStamps::torsionMembers( event* the_event ){

  if( the_event->evt.mbrs.nMembers == 4 ){
    
    current_torsion->members( the_event->evt.mbrs.memberList[0],
                              the_event->evt.mbrs.memberList[1],
                              the_event->evt.mbrs.memberList[2],
                              the_event->evt.mbrs.memberList[3]);
    return 1;
    
  } else {
    the_event->err_msg = strdup( "MakeStamp error. Wrong number of members "
                                 " in torsion");
    return 0;

  }
}

int MakeStamps::torsionConstraint( event* the_event ){

  current_torsion->constrain( the_event->evt.cnstr );
  return 1;
}

int MakeStamps::torsionEnd( event* the_event ){

  the_event->err_msg = current_torsion->checkMe();
  if( the_event->err_msg != NULL ) return 0;
  
  return 1;
}
