#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Globals.hpp"
#include "simError.h"
#ifdef IS_MPI
#include "mpiBASS.h"
#endif // is_mpi

/*
 * The following section lists all of the defined tokens for the
 * gloabal assignment statements. All are prefixed with a G_ to avoid
 * stepping on any previously defined enumerations. 
 *
 * NOTE: tokens start at 1, 0 is a resrved token number
 */

//required parameters
#define G_FORCEFIELD         1
#define G_NCOMPONENTS        2
#define G_TARGETTEMP         3
#define G_ENSEMBLE           4
#define G_DT                 5
#define G_RUNTIME            6

//optional parameters
#define G_INITIALCONFIG      7
#define G_FINALCONFIG        8
#define G_NMOL               9
#define G_DENSITY           10
#define G_BOX               11
#define G_BOXX              12
#define G_BOXY              13
#define G_BOXZ              14
#define G_SAMPLETIME        15
#define G_STATUSTIME        16
#define G_RCUT              17
#define G_RSW               18
#define G_DIELECTRIC        19
#define G_TEMPSET           20
#define G_THERMALTIME       21
#define G_USEPBC            22
#define G_MIXINGRULE        23
#define G_USERF             24
#define G_TARGETPRESSURE    25
#define G_TAUTHERMOSTAT     26
#define G_TAUBAROSTAT       27
#define G_ZCONSTIME         28
#define G_NZCONSTRAINTS     29
#define G_ZCONSTOL          30
#define G_ZCONSFORCEPOLICY  31
#define G_SEED              32
#define G_RESETTIME         33
#define G_USEINITTIME       34
#define G_USEINIT_XS_STATE  35
#define G_ORTHOBOXTOLERANCE 36
#define G_MINIMIZER         37
#define G_MIN_MAXITER       38
#define G_MIN_WRITEFRQ      39
#define G_MIN_STEPSIZE      40
#define G_MIN_FTOL          41
#define G_MIN_GTOL          42
#define G_MIN_LSTOL         43
#define G_MIN_LSMAXITER     44
#define G_ZCONSGAP          45
#define G_ZCONSFIXTIME      46
#define G_ZCONSUSINGSMD     47
#define G_USE_SOLID_THERM_INT     48
#define G_USE_LIQUID_THERM_INT    49
#define G_THERM_INT_LAMBDA  50
#define G_THERM_INT_K       51
#define G_FORCEFIELD_VARIANT 52

Globals::Globals(){
  initalize();
}

Globals::~Globals(){
  int i;

  for( i=0; i<hash_size; i++ ){
    if( command_table[i] != NULL ) delete command_table[i];
  }
  delete[] command_table;

  if( components != NULL ){
    for( i=0; i<n_components; i++ ) delete components[i];
    delete[] components;
  }
}

void Globals::initalize(){
  int i;
  
  hash_size = 23;
  hash_shift = 4;
  
  components = NULL;
  
  command_table = new LinkedCommand*[hash_size];
  for( i=0; i<hash_size; i++ ) command_table[i] = NULL;
  
  addHash( "forceField",    G_FORCEFIELD );
  addHash( "nComponents",   G_NCOMPONENTS );
  addHash( "targetTemp",    G_TARGETTEMP );
  addHash( "ensemble",      G_ENSEMBLE );
  
  addHash( "dt",            G_DT );
  addHash( "runTime",       G_RUNTIME );
  
  addHash( "initialConfig", G_INITIALCONFIG );
  addHash( "finalConfig",   G_FINALCONFIG );
  addHash( "nMol",          G_NMOL );
  addHash( "density",       G_DENSITY );
  addHash( "box",           G_BOX );
  addHash( "boxX",          G_BOXX );
  addHash( "boxY",          G_BOXY );
  addHash( "boxZ",          G_BOXZ );
  addHash( "sampleTime",    G_SAMPLETIME );
  addHash( "resetTime",     G_RESETTIME );
  addHash( "statusTime",    G_STATUSTIME );
  addHash( "cutoffRadius",  G_RCUT );
  addHash( "switchingRadius",  G_RSW );
  addHash( "dielectric",    G_DIELECTRIC );
  addHash( "tempSet",       G_TEMPSET );
  addHash( "thermalTime",   G_THERMALTIME );
  addHash( "mixingRule",    G_MIXINGRULE);
  addHash( "usePeriodicBoundaryConditions",        G_USEPBC);
  addHash( "useReactionField",                     G_USERF );
  addHash( "targetPressure",                       G_TARGETPRESSURE);
  addHash( "tauThermostat",                        G_TAUTHERMOSTAT);
  addHash( "tauBarostat",                          G_TAUBAROSTAT);
  addHash( "zconsTime",                            G_ZCONSTIME);
  addHash( "nZconstraints",                        G_NZCONSTRAINTS);
  addHash( "zconsTol",                             G_ZCONSTOL);
  addHash( "zconsForcePolicy",                     G_ZCONSFORCEPOLICY);
  addHash( "seed",                                 G_SEED);
  addHash( "useInitialTime",                       G_USEINITTIME);
  addHash( "useInitialExtendedSystemState",        G_USEINIT_XS_STATE);
  addHash( "orthoBoxTolerance",                    G_ORTHOBOXTOLERANCE);
  addHash( "minimizer",                            G_MINIMIZER);
  addHash( "minimizerMaxIter",                     G_MIN_MAXITER);
  addHash( "minimizerWriteFrq",                    G_MIN_WRITEFRQ);
  addHash( "minimizerStepSize",                    G_MIN_STEPSIZE);
  addHash( "minimizerFTol",                        G_MIN_FTOL);
  addHash( "minimizerGTol",                        G_MIN_GTOL);
  addHash( "minimizerLSTol",                       G_MIN_LSTOL);
  addHash( "minimizerLSMaxIter",                   G_MIN_LSMAXITER);
  addHash( "zconsGap",                             G_ZCONSGAP);
  addHash( "zconsFixtime",                         G_ZCONSFIXTIME);
  addHash( "zconsUsingSMD",                        G_ZCONSUSINGSMD);
  addHash( "useSolidThermInt",                     G_USE_SOLID_THERM_INT);
  addHash( "useLiquidThermInt",                    G_USE_LIQUID_THERM_INT);
  addHash( "thermodynamicIntegrationLambda",       G_THERM_INT_LAMBDA);
  addHash( "thermodynamicIntegrationK",            G_THERM_INT_K);
  addHash( "forceFieldVariant",                    G_FORCEFIELD_VARIANT);

  strcpy( mixingRule,"standard");  //default mixing rules to standard.
  usePBC = 1; //default  periodic boundry conditions to on
  useRF  = 0;
  useInitTime = 0; // default to pull init time from the init file
  useInitXSstate = 0; // default to pull the extended state from the init file
  orthoBoxTolerance = 1E-6;
  useSolidThermInt = 0; // default solid-state thermodynamic integration to off
  useLiquidThermInt = 0; // default liquid thermodynamic integration to off

  have_force_field =  0;
  have_n_components = 0;
  have_target_temp =  0;
  have_ensemble =     0;
  have_dt =           0;
  have_run_time =     0;
  
  have_initial_config = 0;
  have_final_config =   0;
  have_n_mol =          0;
  have_density =        0;
  have_box =            0;
  have_box_x =          0;
  have_box_y =          0;
  have_box_y =          0;
  have_box_z =          0;
  have_sample_time =    0;
  have_status_time =    0;
  have_reset_time =     0;
  have_thermal_time =   0;
  have_rcut =           0;
  have_rsw =            0;
  have_dielectric =     0;
  have_tempSet =        0;
  have_target_pressure =0;
  have_q_mass =         0;
  have_tau_thermostat = 0;
  have_tau_barostat   = 0;
  have_zcons_time     = 0;
  have_n_zConstraints = 0;
  have_zConstraints   = 0;
  have_zcons_tol = 0;
  have_zcons_gap = 0;
  have_zcons_fixtime = 0;
  have_zcons_using_smd = 0;  
  have_seed = 0;
  have_minimizer = 0;
  have_minimizer_maxiteration = 0;
  have_minimizer_writefrq = 0;
  have_minimizer_stepsize = 0;
  have_minimizer_ftol = 0;
  have_minimizer_gtol = 0;
  have_minimizer_ls_tol = 0;
  have_minimizer_ls_maxiteration = 0;
  have_thermodynamic_integration_lambda = 0;
  have_thermodynamic_integration_k = 0;
  have_forcefield_variant = 0;

}

int Globals::newComponent( event* the_event ){
  
  current_component = new Component;
  int index = the_event->evt.blk_index;
  char err[200];
  
  if( have_n_components && index < n_components ) 
    components[index] = current_component;
  else{
    if( have_n_components ){
      sprintf( err, "meta-data parsing error: %d out of nComponents range", 
	       index );
      the_event->err_msg = strdup( err );
      return 0;
    }
    else{
      the_event->err_msg = strdup("meta-data parsing error: nComponents not given before"
				  " first component declaration." );
      return 0;
    }
  }  
 
  return 1;
}



int Globals::componentAssign( event* the_event ){

  switch( the_event->evt.asmt.asmt_type ){
    
  case STRING:
    return current_component->assignString( the_event->evt.asmt.lhs,
					    the_event->evt.asmt.rhs.sval,
					    &(the_event->err_msg));
    break;
    
  case DOUBLE:
    return current_component->assignDouble( the_event->evt.asmt.lhs,
					    the_event->evt.asmt.rhs.dval,
					    &(the_event->err_msg));
    break;
    
  case INT:
    return current_component->assignInt( the_event->evt.asmt.lhs,
					 the_event->evt.asmt.rhs.ival,
					 &(the_event->err_msg));
    break;
    
  default:
    the_event->err_msg = strdup( "Globals error. Invalid component"
				 " assignment type" );
    return 0;
    break;
  }
  return 0;
}

int Globals::componentEnd( event* the_event ){

  the_event->err_msg = current_component->checkMe();
  if( the_event->err_msg != NULL ) return 0;

  return 1;
}

int Globals::newZconstraint( event* the_event ){
  

  int index = the_event->evt.blk_index;
  char err[200];
  current_zConstraint = new ZconStamp( index );
  
  have_zConstraints = 1;

  if( have_n_zConstraints && index < n_zConstraints ) 
    zConstraints[index] = current_zConstraint;
  else{
    if( have_n_zConstraints ){
      sprintf( err, "meta-data parsing error: %d out of nZconstraints range", 
	       index );
      the_event->err_msg = strdup( err );
      return 0;
    }
    else{
      the_event->err_msg = strdup("meta-data parsing error: nZconstraints"
				  " not given before"
				  " first zConstraint declaration." );
      return 0;
    }
  }  
 
  return 1;
}



int Globals::zConstraintAssign( event* the_event ){

  switch( the_event->evt.asmt.asmt_type ){
    
  case STRING:
    return current_zConstraint->assignString( the_event->evt.asmt.lhs,
					      the_event->evt.asmt.rhs.sval,
					      &(the_event->err_msg));
    break;
    
  case DOUBLE:
    return current_zConstraint->assignDouble( the_event->evt.asmt.lhs,
					      the_event->evt.asmt.rhs.dval,
					      &(the_event->err_msg));
    break;
    
  case INT:
    return current_zConstraint->assignInt( the_event->evt.asmt.lhs,
					   the_event->evt.asmt.rhs.ival,
					   &(the_event->err_msg));
    break;
    
  default:
    the_event->err_msg = strdup( "Globals error. Invalid zConstraint"
				 " assignment type" );
    return 0;
    break;
  }
  return 0;
}

int Globals::zConstraintEnd( event* the_event ){

  the_event->err_msg = current_zConstraint->checkMe();
  if( the_event->err_msg != NULL ) return 0;

  return 1;
}

int Globals::globalAssign( event* the_event ){
  
  int key;
  int token;
  interface_assign_type the_type =  the_event->evt.asmt.asmt_type;
  char* lhs = the_event->evt.asmt.lhs;
  char err[300];
  
  token = 0;
  key = hash( lhs );
  if( command_table[key] != NULL ) token = command_table[key]->match( lhs );
  
  if( token ){
    
    switch( token ){
      
    case G_FORCEFIELD:
      if( the_type == STRING ){
	strcpy( force_field, the_event->evt.asmt.rhs.sval );
	have_force_field = 1;
	return 1;
      }
      
      the_event->err_msg = 
	strdup( "Error in parsing meta-data file!\n\tforceField was not a string assignment.\n" );
      return 0;
      break;
      
    case G_NCOMPONENTS:
      if( the_type == STRING ){
	the_event->err_msg = 
	  strdup("Error in parsing meta-data file!\n\tnComponents is not a double or an int.\n" );
	return 0;
      }
      
      else if( the_type == DOUBLE ){
	n_components = (int)the_event->evt.asmt.rhs.dval;
	components = new Component*[n_components];
	have_n_components = 1;
	return 1;
      }
      
      else{
	n_components = the_event->evt.asmt.rhs.ival;
	components = new Component*[n_components];
	have_n_components = 1;
	return 1;
      }
      break;

    case G_NZCONSTRAINTS:
      if( the_type == STRING ){
	the_event->err_msg = 
	  strdup("Error in parsing meta-data file!\n\tnZconstraints is not a double or an int.\n" );
	return 0;
      }
      
      else if( the_type == DOUBLE ){
	n_zConstraints = (int)the_event->evt.asmt.rhs.dval;
	zConstraints = new ZconStamp*[n_zConstraints];
	have_n_zConstraints = 1;
	return 1;
      }
      
      else{
	n_zConstraints = the_event->evt.asmt.rhs.ival;
	zConstraints = new ZconStamp*[n_zConstraints];
	have_n_zConstraints = 1;
	return 1;
      }
      break;
      
    case G_TARGETTEMP:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\ttargetTemp is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	target_temp = the_event->evt.asmt.rhs.dval;
	have_target_temp = 1;
	return 1;
	break;
	
      case INT:
	target_temp = (double)the_event->evt.asmt.rhs.ival;
	have_target_temp = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\ttargetTemp unrecognized.\n" );
	return 0;
	break;
      }
      break;

    case G_ORTHOBOXTOLERANCE:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\torthoBoxTolerance is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	orthoBoxTolerance = the_event->evt.asmt.rhs.dval;
	have_target_temp = 1;
	return 1;
	break;
	
      case INT:
	orthoBoxTolerance = (double)the_event->evt.asmt.rhs.ival;
	have_target_temp = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Global error.orthoBoxTolerance unrecognized.\n" );
	return 0;
	break;
      }
      break;
      
    case G_ENSEMBLE:
      if( the_type == STRING ){
	strcpy( ensemble, the_event->evt.asmt.rhs.sval );
	have_ensemble = 1;
	return 1;
      }
      
      the_event->err_msg = 
	strdup( "Error in parsing meta-data file!\n\tensemble was not assigned to a string\n" );
      return 0;
      break;      

    case G_MIXINGRULE:
      if( the_type == STRING ){
	strcpy( mixingRule, the_event->evt.asmt.rhs.sval );
	return 1;
      }
      
      the_event->err_msg = 
	strdup( "Error in parsing meta-data file!\n\tmixing rule was not assigned to a string\n" );
      return 0;
      break;      
      
    case G_DT:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tdt is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	dt = the_event->evt.asmt.rhs.dval;
	have_dt = 1;
	return 1;
	break;
	
      case INT:
	dt = (double)the_event->evt.asmt.rhs.ival;
	have_dt = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tdt unrecognized.\n" );
	return 0;
	break;
      }
      break;
      
    case G_RUNTIME:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\trunTime is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	run_time = the_event->evt.asmt.rhs.dval;
	have_run_time = 1;
	return 1;
	break;
	
      case INT:
	run_time = (double)the_event->evt.asmt.rhs.ival;
	have_run_time = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\trunTime unrecognized.\n" );
	return 0;
	break;
      }
      break;
      
    case G_INITIALCONFIG:
      if( the_type == STRING ){
	strcpy( initial_config, the_event->evt.asmt.rhs.sval );
	have_initial_config = 1;
	return 1;
      }
      
      the_event->err_msg = 
	strdup( "Error in parsing meta-data file!\n\tinitialConfig was not a string assignment.\n" );
      return 0;
      break;
      
    case G_FINALCONFIG:
      if( the_type == STRING ){
	strcpy( final_config, the_event->evt.asmt.rhs.sval );
	have_final_config = 1;
	return 1;
      }
      
      the_event->err_msg = 
	strdup( "Error in parsing meta-data file!\n\tfinalConfig was not a string assignment.\n" );
      return 0;
      break;
      
    case G_NMOL:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tnMol is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	n_mol = (int)the_event->evt.asmt.rhs.dval;
	have_n_mol = 1;
	return 1;
	break;
	
      case INT:
	n_mol = the_event->evt.asmt.rhs.ival;
	have_n_mol = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tnMol unrecognized.\n" );
	return 0;
	break;
      }
      break;
      
    case G_DENSITY:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tdensity is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	density = the_event->evt.asmt.rhs.dval;
	have_density = 1;
	return 1;
	break;
	
      case INT:
	density = (double)the_event->evt.asmt.rhs.ival;
	have_density = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tdensity unrecognized.\n" );
	return 0;
	break;
      }
      break;
      
    case G_BOX:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tbox is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	box = the_event->evt.asmt.rhs.dval;
	have_box = 1;
	return 1;
	break;
	
      case INT:
	box = (double)the_event->evt.asmt.rhs.ival;
	have_box = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tbox unrecognized.\n" );
	return 0;
	break;
      }
      break;
      
    case G_BOXX:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tboxX is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	box_x = the_event->evt.asmt.rhs.dval;
	have_box_x = 1;
	return 1;
	break;
	
      case INT:
	box_x = (double)the_event->evt.asmt.rhs.ival;
	have_box_x = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tboxX unrecognized.\n" );
	return 0;
	break;
      }
      break;
      
    case G_BOXY:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tboxY is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	box_y = the_event->evt.asmt.rhs.dval;
	have_box_y = 1;
	return 1;
	break;
	
      case INT:
	box_y = (double)the_event->evt.asmt.rhs.ival;
	have_box_y = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tboxY unrecognized.\n" );
	return 0;
	break;
      }
      break;
      
    case G_BOXZ:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tboxZ is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	box_z = the_event->evt.asmt.rhs.dval;
	have_box_z = 1;
	return 1;
	break;
	
      case INT:
	box_z = (double)the_event->evt.asmt.rhs.ival;
	have_box_z = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tboxZ unrecognized.\n" );
	return 0;
	break;
      }
      break;
      
    case G_SAMPLETIME:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tsampleTime is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	sample_time = the_event->evt.asmt.rhs.dval;
	have_sample_time = 1;
	return 1;
	break;
	
      case INT:
	sample_time = (double)the_event->evt.asmt.rhs.ival;
	have_sample_time = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tsampleTime unrecognized.\n" );
	return 0;
	break;
      }
      break;
      
    case G_STATUSTIME:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tstatusTime is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	status_time = the_event->evt.asmt.rhs.dval;
	have_status_time = 1;
	return 1;
	break;
	
      case INT:
	status_time = (double)the_event->evt.asmt.rhs.ival;
	have_status_time = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tstatusTime unrecognized.\n" );
	return 0;
	break;
      }
      break;

    case G_RESETTIME:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tresetTime is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	resetTime = the_event->evt.asmt.rhs.dval;
	have_reset_time = 1;
	return 1;
	break;
	
      case INT:
	resetTime = (double)the_event->evt.asmt.rhs.ival;
	have_reset_time = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tresetTime unrecognized.\n" );
	return 0;
	break;
      }
      break;

    case G_THERMALTIME:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tthermalTime is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	thermal_time = the_event->evt.asmt.rhs.dval;
	have_thermal_time = 1;
	return 1;
	break;
	
      case INT:
	thermal_time = (double)the_event->evt.asmt.rhs.ival;
	have_thermal_time = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tthermalTime unrecognized.\n" );
	return 0;
	break;
      }
      break;
      
    case G_RCUT:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tcutoffRadius is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	rcut = the_event->evt.asmt.rhs.dval;
	have_rcut = 1;
	return 1;
	break;
	
      case INT:
	rcut = (double)the_event->evt.asmt.rhs.ival;
	have_rcut = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tcutoffRadius unrecognized.\n" );
	return 0;
	break;
      }
      break;

    case G_RSW:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tswitchingRadius is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	rsw = the_event->evt.asmt.rhs.dval;
	have_rsw = 1;
	return 1;
	break;
	
      case INT:
	rsw = (double)the_event->evt.asmt.rhs.ival;
	have_rsw = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tswitchingRadius unrecognized.\n" );
	return 0;
	break;
      }
      break;
      
    case G_DIELECTRIC:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tdielectric is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	dielectric = the_event->evt.asmt.rhs.dval;
	have_dielectric = 1;
	return 1;
	break;
	
      case INT:
	dielectric = (double)the_event->evt.asmt.rhs.ival;
	have_dielectric = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tdielectric unrecognized.\n" );
	return 0;
	break;
      }
      break;
      
    case G_TEMPSET:
      if( the_type == STRING ){
	
	if( !strcasecmp( "true", the_event->evt.asmt.rhs.sval )) tempSet = 1;
	else if( !strcasecmp( "false", the_event->evt.asmt.rhs.sval )) tempSet = 0;
	else{
	  the_event->err_msg = 
	    strdup( "Error in parsing meta-data file!\n\ttempSet was not \"true\" or \"false\".\n" );
	  return 0;
	}
	have_tempSet = 1;
	return 1;
      }
      
      the_event->err_msg = 
	strdup( "Error in parsing meta-data file!\n\ttempSet was not \"true\" or \"false\".\n" );
      return 0;
      break;

    case G_USEINITTIME:
      if( the_type == STRING ){
	
	if( !strcasecmp( "true", the_event->evt.asmt.rhs.sval )) useInitTime = 1;
	else if( !strcasecmp( "false", the_event->evt.asmt.rhs.sval )) useInitTime = 0;
	else{
	  the_event->err_msg = 
	    strdup( "Error in parsing meta-data file!\n\tuseInitTime was not \"true\" or \"false\".\n" );
	  return 0;
	}
	return 1;
      }
      
      the_event->err_msg = 
	strdup( "Error in parsing meta-data file!\n\tuseInitTime was not \"true\" or \"false\".\n" );
      return 0;
      break;

    case G_USEINIT_XS_STATE:
      if( the_type == STRING ){
	
	if( !strcasecmp( "true", the_event->evt.asmt.rhs.sval )) 
	  useInitXSstate = 1;
	else if( !strcasecmp( "false", the_event->evt.asmt.rhs.sval )) 
	  useInitXSstate = 0;
	else{
	  the_event->err_msg = 
	    strdup( "Error in parsing meta-data file!\n\tuseInitExtendedSystemState was not \"true\" or \"false\".\n" );
	  return 0;
	}
	return 1;
      }
      
      the_event->err_msg = 
	strdup( "Error in parsing meta-data file!\n\tuseInitExtendedSystemState was not \"true\" or \"false\".\n" );
      return 0;
      break;
      
    case G_USEPBC:
      if( the_type == STRING ){
	
	if( !strcasecmp( "true", the_event->evt.asmt.rhs.sval )) usePBC = 1;
	else if( !strcasecmp( "false", the_event->evt.asmt.rhs.sval )) usePBC = 0;
	else{
	  the_event->err_msg = 
	    strdup( "Error in parsing meta-data file!\n\tusePeriodicBoundaryConditions was not \"true\" or \"false\".\n" );
	  return 0;
	}
	return 1;
      }
      
      the_event->err_msg = 
	strdup( "Error in parsing meta-data file!\n\tusePeriodicBoundaryConditions was not \"true\" or \"false\".\n" );
      return 0;
      break;

    case G_USERF:
      if( the_type == STRING ){
	
	if( !strcasecmp( "true", the_event->evt.asmt.rhs.sval )) useRF = 1;
	else if( !strcasecmp( "false", the_event->evt.asmt.rhs.sval )) useRF = 0;
	else{
	  the_event->err_msg = 
	    strdup( "Error in parsing meta-data file!\n\tuseReactionField was not \"true\" or \"false\".\n" );
	  return 0;
	}
	return 1;
      }
      
      the_event->err_msg = 
	strdup( "Error in parsing meta-data file!\n\tuseReactionField was not \"true\" or \"false\".\n" );
      return 0;
      break;

    case G_TARGETPRESSURE:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\ttargetPressure is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	target_pressure = the_event->evt.asmt.rhs.dval;
	have_target_pressure = 1;
	return 1;
	break;
	
      case INT:
	target_pressure = (double)the_event->evt.asmt.rhs.ival;
	have_target_pressure = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\ttargetPressure unrecognized.\n" );
	return 0;
	break;
      }
      break;

    case G_TAUTHERMOSTAT:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\ttauThermostat is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	tau_thermostat = the_event->evt.asmt.rhs.dval;
	have_tau_thermostat = 1;
	return 1;
	break;
	
      case INT:
	tau_thermostat = (double)the_event->evt.asmt.rhs.ival;
	have_tau_thermostat = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\ttauThermostat unrecognized.\n" );
	return 0;
	break;
      }
      break;

    case G_TAUBAROSTAT:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\ttauBarostat is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	tau_barostat = the_event->evt.asmt.rhs.dval;
	have_tau_barostat = 1;
	return 1;
	break;
	
      case INT:
	tau_barostat = (double)the_event->evt.asmt.rhs.ival;
	have_tau_barostat = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\ttauBarostat unrecognized.\n" );
	return 0;
	break;
      }
      break;
      
    case G_ZCONSTIME:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tzcons_time is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	zcons_time = the_event->evt.asmt.rhs.dval;
	have_zcons_time = 1;
	return 1;
	break;
	
      case INT:
	zcons_time = (double)the_event->evt.asmt.rhs.ival;
	have_zcons_time = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tzcons_time unrecognized.\n" );
	return 0;
	break;
      }
      break;

    case G_ZCONSTOL:
      switch( the_type ){
	
      case STRING:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tzcons_tol is not a double or int.\n" );
	return 0;
	break;
	
      case DOUBLE:
	zcons_tol = the_event->evt.asmt.rhs.dval;
	have_zcons_tol = 1;
	return 1;
	break;
	
      case INT:
	zcons_tol = (double)the_event->evt.asmt.rhs.ival;
	have_zcons_tol = 1;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tzcons_ol unrecognized.\n" );
	return 0;
	break;
      }
      break;
   
    case G_ZCONSFORCEPOLICY:
      switch( the_type ){
	
      case STRING:
   strcpy(zconsForcePolicy, the_event->evt.asmt.rhs.sval);

   for(int i = 0; zconsForcePolicy[i] != '\0'; i++)
	{
      zconsForcePolicy[i] = toupper(zconsForcePolicy[i]);
   }
	have_zcons_force_policy = 1;
   return 1;
	break;
	
      case DOUBLE:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tzconsForcePolicy is not a double or int.\n" );
	return 0;
	break;
	
      case INT:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tzconsForcePolicy is not a double or int.\n" );
	return 0;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tzconsForcePolicy unrecognized.\n" );
	return 0;
	break;
      }
      break;
      
    case G_ZCONSGAP:
      switch( the_type ){
  
      case STRING:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tzcons_gap is not a double or int.\n" );
        return 0;
        break;
  
      case DOUBLE:
        zcons_gap = the_event->evt.asmt.rhs.dval;
        have_zcons_gap= 1;
        return 1;
        break;
  
      case INT:
        zcons_gap= (double)the_event->evt.asmt.rhs.ival;
        have_zcons_gap= 1;
        return 1;
        break;
  
      default:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tzcons_gap unrecognized.\n" );
        return 0;
        break;
      }
      break;

    case G_ZCONSFIXTIME:
      switch( the_type ){
  
      case STRING:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tzcons_fixtime is not a double or int.\n" );
        return 0;
        break;
  
      case DOUBLE:
        zcons_fixtime= the_event->evt.asmt.rhs.dval;
        have_zcons_fixtime= 1;
        return 1;
        break;
  
      case INT:
        zcons_fixtime= (double)the_event->evt.asmt.rhs.ival;
        have_zcons_fixtime= 1;
        return 1;
        break;
  
      default:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tzcons_fixtime unrecognized.\n" );
        return 0;
        break;
      }
      break;

    case G_ZCONSUSINGSMD:
      switch( the_type ){
  
      case STRING:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tzcons_fixtime is not an  int.\n" );
        return 0;
        break;
  
      case DOUBLE:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tzcons_fixtime is not an  int.\n" );
        return 0;
        break;
  
      case INT:
        zcons_using_smd= the_event->evt.asmt.rhs.ival;
        have_zcons_using_smd= 1;
        return 1;
        break;
  
      default:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tzcons_usingsmd unrecognized.\n" );
        return 0;
        break;
      }
      break;
      
    case G_MINIMIZER:
      switch( the_type ){

      case STRING:
        strcpy(minimizer_name, the_event->evt.asmt.rhs.sval);

        for(int i = 0; zconsForcePolicy[i] != '\0'; i++){
          zconsForcePolicy[i] = toupper(zconsForcePolicy[i]);
        }
        have_minimizer= 1;
        return 1;
        break;
        
      case DOUBLE:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_name is not a double or int.\n" );
        return 0;
        break;
        
      case INT:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_name is not a double or int.\n" );
        return 0;
        break;
        
      default:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_name unrecognized.\n" );
        return 0;
        break;
      }
      break;

    case G_MIN_MAXITER:
      switch( the_type ){

      case STRING:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_maxiteration is not a double or int.\n" );
        return 1;
        break;
        
      case DOUBLE:
        minimizer_maxiteration = (int)the_event->evt.asmt.rhs.dval;
        have_minimizer_maxiteration = 1;
        return 1;
        break;
        
      case INT:
        minimizer_maxiteration = the_event->evt.asmt.rhs.ival;
        have_minimizer_maxiteration = 1;
        return 1;
        break;
        
      default:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_maxiteration unrecognized.\n" );
        return 0;
        break;
      }
      break;
      
    case G_MIN_WRITEFRQ:
      switch( the_type ){

      case STRING:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_writefrq is not a double or int.\n" );
        return 1;
        break;
        
      case DOUBLE:
        minimizer_writefrq= the_event->evt.asmt.rhs.dval;
        have_minimizer_writefrq = 1;
        return 1;
        break;
        
      case INT:
        minimizer_writefrq= the_event->evt.asmt.rhs.ival;
        have_minimizer_writefrq = 1;
        return 1;
        break;
        
      default:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_writefrq unrecognized.\n" );
        return 0;
        break;
      }
      break;

    case G_MIN_STEPSIZE:
      switch( the_type ){

      case STRING:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_resetfrq is not a double or int.\n" );
        return 1;
        break;
        
      case DOUBLE:
        minimizer_stepsize= the_event->evt.asmt.rhs.dval;
        have_minimizer_stepsize = 1;
        return 1;
        break;
        
      case INT:
        minimizer_stepsize= the_event->evt.asmt.rhs.ival;
        have_minimizer_stepsize = 1;
        return 1;
        break;
        
      default:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_resetfrq unrecognized.\n" );
        return 0;
        break;
      }
      break;      

    case G_MIN_FTOL:
      switch( the_type ){

      case STRING:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_ftol is not a double or int.\n" );
        return 1;
        break;
        
      case DOUBLE:
        minimizer_ftol= the_event->evt.asmt.rhs.dval;
        have_minimizer_ftol = 1;
        return 1;
        break;
        
      case INT:
        minimizer_ftol= the_event->evt.asmt.rhs.ival;
        have_minimizer_ftol = 1;
        return 1;
        break;
        
      default:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_ftol unrecognized.\n" );
        return 0;
        break;
      }
      break; 
      
    case G_MIN_GTOL:
      switch( the_type ){

      case STRING:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_gtol is not a double or int.\n" );
        return 1;
        break;
        
      case DOUBLE:
        minimizer_gtol= the_event->evt.asmt.rhs.dval;
        have_minimizer_gtol = 1;
        return 1;
        break;
        
      case INT:
        minimizer_gtol= the_event->evt.asmt.rhs.ival;
        have_minimizer_gtol = 1;
        return 1;
        break;
        
      default:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_gtol unrecognized.\n" );
        return 0;
        break;
      }
      break; 
      
    case G_MIN_LSMAXITER:
      switch( the_type ){

      case STRING:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_ls_maxiteration is not a double or int.\n" );
        return 1;
        break;
        
      case DOUBLE:
        minimizer_ls_maxiteration = the_event->evt.asmt.rhs.dval;
        have_minimizer_ls_maxiteration = 1;
        return 1;
        break;
        
      case INT:
        minimizer_ls_maxiteration = the_event->evt.asmt.rhs.ival;
        have_minimizer_ls_maxiteration = 1;
        return 1;
        break;
        
      default:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_ls_maxiteration unrecognized.\n" );
        return 0;
        break;
      }
      break;      

    case G_MIN_LSTOL:
      switch( the_type ){

      case STRING:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_ls_tol is not a double or int.\n" );
        return 1;
        break;
        
      case DOUBLE:
        minimizer_ls_tol= the_event->evt.asmt.rhs.dval;
        have_minimizer_ls_tol = 1;
        return 1;
        break;
        
      case INT:
        minimizer_ls_tol= the_event->evt.asmt.rhs.ival;
        have_minimizer_ls_tol = 1;
        return 1;
        break;
        
      default:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tminimizer_ls_tol unrecognized.\n" );
        return 0;
        break;
      }
      break; 
      
      // add more token cases here.
    case G_SEED:
      switch( the_type ){
	
      case STRING:
   the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tseed is not a string.\n" );
	return 0;
   return 0;
	break;
	
      case DOUBLE:
   have_seed = 1;
   seed = (int)the_event->evt.asmt.rhs.dval;
	return 1;
	break;
	
      case INT:
   have_seed = 1;
   seed =  the_event->evt.asmt.rhs.ival ;
	return 1;
	break;
	
      default:
	the_event->err_msg = 
	  strdup( "Error in parsing meta-data file!\n\tseed unrecognized.\n" );
	return 0;
	break;
      }
      break;

    case G_USE_SOLID_THERM_INT:
      if( the_type == STRING ){
	
	if( !strcasecmp( "true", the_event->evt.asmt.rhs.sval )) useSolidThermInt = 1;
	else if( !strcasecmp( "false", the_event->evt.asmt.rhs.sval )) useSolidThermInt = 0;
	else{
	  the_event->err_msg = 
	    strdup( "Error in parsing meta-data file!\n\tuseSolidThermInt was not \"true\" or \"false\".\n" );
	  return 0;
	}
	return 1;
      }
      
      the_event->err_msg = 
	strdup( "Error in parsing meta-data file!\n\tuseSolidThermInt was not \"true\" or \"false\".\n" );
      return 0;
      break;

    case G_USE_LIQUID_THERM_INT:
      if( the_type == STRING ){
	
	if( !strcasecmp( "true", the_event->evt.asmt.rhs.sval )) useLiquidThermInt = 1;
	else if( !strcasecmp( "false", the_event->evt.asmt.rhs.sval )) useLiquidThermInt = 0;
	else{
	  the_event->err_msg = 
	    strdup( "Error in parsing meta-data file!\n\tuseLiquidThermInt was not \"true\" or \"false\".\n" );
	  return 0;
	}
	return 1;
      }
      
      the_event->err_msg = 
	strdup( "Error in parsing meta-data file!\n\tuseLiquidThermInt was not \"true\" or \"false\".\n" );
      return 0;
      break;

    case G_THERM_INT_LAMBDA:
      switch( the_type ){
	
      case STRING:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tthermodynamicIntegrationLambda is not a double or int.\n" );
        return 1;
        break;
        
      case DOUBLE:
        thermodynamic_integration_lambda = the_event->evt.asmt.rhs.dval;
        have_thermodynamic_integration_lambda = 1;
        return 1;
        break;
        
      case INT:
        thermodynamic_integration_lambda = (double)the_event->evt.asmt.rhs.dval;
        have_thermodynamic_integration_lambda = 1;
        return 1;
        break;
        
      default:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tthermodynamicIntegrationLambda unrecognized.\n" );
        return 0;
        break;
      }
      break;      

    case G_THERM_INT_K:
      switch( the_type ){
	
      case STRING:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tthermodynamicIntegrationK is not a double or int.\n" );
        return 1;
        break;
        
      case DOUBLE:
        thermodynamic_integration_k = the_event->evt.asmt.rhs.dval;
        have_thermodynamic_integration_k = 1;
        return 1;
        break;
        
      case INT:
        thermodynamic_integration_k = (double)the_event->evt.asmt.rhs.dval;
        have_thermodynamic_integration_k = 1;
        return 1;
        break;
        
      default:
        the_event->err_msg = 
          strdup( "Error in parsing meta-data file!\n\tthermodynamicIntegrationK unrecognized.\n" );
        return 0;
        break;
      }
      break;      
    case G_FORCEFIELD_VARIANT:
      if( the_type == STRING ){
	strcpy( forcefield_variant, the_event->evt.asmt.rhs.sval );
	have_forcefield_variant = 1;
	return 1;
      }
      
      the_event->err_msg = 
	strdup( "Error in parsing meta-data file!\n\tforceFieldVariant was not a string assignment.\n" );
      return 0;
      break;      
      // add more token cases here.      
    }
  }
  
  switch( the_type ){
    
  case STRING:
    sprintf( err,
	     "\tUnrecognized assignment:\n"
	     "\t\t-> %s = %s\n",
	     lhs, the_event->evt.asmt.rhs.sval );
    break;
    
  case DOUBLE:
    sprintf( err,
	     "\tUnrecognized assignment:\n"
	     "\t\t-> %s = %lf\n",
	     lhs, the_event->evt.asmt.rhs.dval );
    break;
    
  case INT:
    sprintf( err,
	     "\tUnrecognized assignment:\n"
	     "\t\t-> %s = %d\n",
	     lhs, the_event->evt.asmt.rhs.ival );
    break;
    
  default:
    sprintf( err,
	     "\tUnrecognized assignment:\n"
	     "\t\t-> %s = ?\n",
	     lhs );
    break;
  }
  
  the_event->err_msg = strdup( err );
  return 0;
}

char* Globals::checkMe( void ){
  
  char err[300];
  short int have_err = 0;

  strcpy( err, "The following required keywords are missing:\n" );
  
  if( !have_force_field ){
    strcat( err, "\t\t->forceField\n" );
    have_err= 1;
  }

  if( !have_n_components ){
    strcat( err, "\t\t->nComponents\n" );
    have_err= 1;
  }

  
  if ( !have_ensemble ) {
    // I'm not doing MD:
    if ( !have_minimizer ) {
      // I'm not doing MD or minimization:
      strcat( err, "\t\t->either ensemble or minimizer must be set!\n" );
      have_err = 1;      
    } else {
      // I'm a minimizer:
    }
  } else {
    // I *am* doing MD:
    if( !have_dt ){
      strcat( err, "\t\t->dt (timestep in fs)\n" );
      have_err= 1;
    } 
    if( !have_run_time ){
      strcat( err, "\t\t->runTime (total run time in fs)\n" );
      have_err= 1;
    }    
  }
      
  if( have_err ) return strdup( err );
  
  return NULL;
}

int Globals::globalEnd( event* the_event ){
  
  the_event->err_msg = checkMe();
  if( the_event->err_msg != NULL ) return 0;

  return 1;
}

int Globals::hash( char* text ){

  register unsigned short int i = 0; // loop counter
  int key = 0; // the hash key
  
  while( text[i] != '\0' ){
    
    key = ( ( key << hash_shift ) + text[i] ) % hash_size;
    
    i++;
  }
  
  if( key < 0 ){

    // if the key is less than zero, we've had an overflow error

    sprintf( painCave.errMsg,
	     "There has been an overflow error in the Globals' hash key.");
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

void Globals::addHash( char* text, int token ){

  int key;
  LinkedCommand* the_element;

  the_element = new LinkedCommand;
  the_element->setValues( text, token );

  key = hash( text );

  the_element->setNext( command_table[key] );
  command_table[key] = the_element;
}
