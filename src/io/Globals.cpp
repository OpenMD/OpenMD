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
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "io/Globals.hpp"
#include "utils/simError.h"
#ifdef IS_MPI
#include "io/mpiBASS.h"
#endif // is_mpi

#include "io/ParamConstraint.hpp"


Globals::Globals(){
 
  DefineParameter(ForceField, "forceField")
  DefineParameter(NComponents, "nComponents")
 
  DefineOptionalParameter(TargetTemp, "targetTemp");
  DefineOptionalParameter(Ensemble, "ensemble");
  DefineOptionalParameter(Dt, "dt");
  DefineOptionalParameter(RunTime, "runTime");
  DefineOptionalParameter(InitialConfig, "initialConfig");
  DefineOptionalParameter(FinalConfig, "finalConfig");
  DefineOptionalParameter(NMol, "nMol");
  DefineOptionalParameter(Density, "density");
  DefineOptionalParameter(Box, "box");
  DefineOptionalParameter(BoxX, "boxX");
  DefineOptionalParameter(BoxY, "boxY");
  DefineOptionalParameter(BoxZ, "boxZ");
  DefineOptionalParameter(SampleTime, "sampleTime");
  DefineOptionalParameter(ResetTime, "resetTime");
  DefineOptionalParameter(StatusTime, "statusTime");
  DefineOptionalParameter(CutoffRadius, "cutoffRadius");
  DefineOptionalParameter(SwitchingRadius, "switchingRadius");
  DefineOptionalParameter(Dielectric, "dielectric");
  DefineOptionalParameter(TempSet, "tempSet");
  DefineOptionalParameter(ThermalTime, "thermalTime");
  DefineOptionalParameter(TargetPressure, "targetPressure");
  DefineOptionalParameter(TauThermostat, "tauThermostat");
  DefineOptionalParameter(TauBarostat, "tauBarostat");
  DefineOptionalParameter(ZconsTime, "zconsTime");
  DefineOptionalParameter(NZconstraints, "nZconstraints");  
  DefineOptionalParameter(ZconsTol, "zconsTol");
  DefineOptionalParameter(ZconsForcePolicy, "zconsForcePolicy");
  DefineOptionalParameter(Seed, "seed");
  DefineOptionalParameter(Minimizer, "minimizer");
  DefineOptionalParameter(MinimizerMaxIter,"minimizerMaxIter");
  DefineOptionalParameter(MinimizerWriteFrq, "minimizerWriteFrq");
  DefineOptionalParameter(MinimizerStepSize, "minimizerStepSize");
  DefineOptionalParameter(MinimizerFTol, "minimizerFTol");
  DefineOptionalParameter(MinimizerGTol, "minimizerGTol");
  DefineOptionalParameter(MinimizerLSTol, "minimizerLSTol");
  DefineOptionalParameter(MinimizerLSMaxIter, "minimizerLSMaxIter");
  DefineOptionalParameter(ZconsGap, "zconsGap");
  DefineOptionalParameter(ZconsFixtime, "zconsFixtime");
  DefineOptionalParameter(ZconsUsingSMD, "zconsUsingSMD");
  DefineOptionalParameter(ThermodynamicIntegrationLambda, "thermodynamicIntegrationLambda");
  DefineOptionalParameter(ThermodynamicIntegrationK, "thermodynamicIntegrationK");
  DefineOptionalParameter(ForceFieldVariant, "forceFieldVariant");
  DefineOptionalParameter(ForceFieldFileName, "forceFieldFileName");
  DefineOptionalParameter(ThermIntDistSpringConst, "thermIntDistSpringConst");
  DefineOptionalParameter(ThermIntThetaSpringConst, "thermIntThetaSpringConst");
  DefineOptionalParameter(ThermIntOmegaSpringConst, "thermIntOmegaSpringConst");
  DefineOptionalParameter(SurfaceTension, "surfaceTension");
  DefineOptionalParameter(PrintPressureTensor, "printPressureTensor");
  DefineOptionalParameter(ElectrostaticSummationMethod, "electrostaticSummationMethod");
  DefineOptionalParameter(ElectrostaticScreeningMethod, "electrostaticScreeningMethod");
  DefineOptionalParameter(CutoffPolicy, "cutoffPolicy");
  DefineOptionalParameter(SwitchingFunctionType, "switchingFunctionType");
  DefineOptionalParameterWithDefaultValue(MixingRule, "mixingRule", "standard");
  DefineOptionalParameterWithDefaultValue(UsePeriodicBoundaryConditions, "usePeriodicBoundaryConditions", true);
  DefineOptionalParameterWithDefaultValue(UseInitalTime, "useInitialTime", false);
  DefineOptionalParameterWithDefaultValue(UseIntialExtendedSystemState, "useInitialExtendedSystemState", false);
  DefineOptionalParameterWithDefaultValue(OrthoBoxTolerance, "orthoBoxTolerance", 1E-6);  
  DefineOptionalParameterWithDefaultValue(UseSolidThermInt, "useSolidThermInt", false);
  DefineOptionalParameterWithDefaultValue(UseLiquidThermInt, "useLiquidThermInt", false);
  DefineOptionalParameterWithDefaultValue(ThermIntDistSpringConst, "thermIntDistSpringConst", 6.0);
  DefineOptionalParameterWithDefaultValue(ThermIntThetaSpringConst, "thermIntThetaSpringConst", 7.5);
  DefineOptionalParameterWithDefaultValue(ThermIntOmegaSpringConst, "thermIntOmegaSpringConst", 13.5);
  DefineOptionalParameterWithDefaultValue(DampingAlpha, "dampingAlpha", 0.2);
  DefineOptionalParameterWithDefaultValue(CompressDumpFile, "compressDumpFile", 0);
  DefineOptionalParameterWithDefaultValue(OutputForceVector, "outputForceVector", 0);
  DefineOptionalParameterWithDefaultValue(SkinThickness, "skinThickness", 1.0);
  DefineOptionalParameterWithDefaultValue(StatFileFormat, "statFileFormat", "TIME|TOTAL_ENERGY|POTENTIAL_ENERGY|KINETIC_ENERGY|TEMPERATURE|PRESSURE|VOLUME|CONSERVED_QUANTITY");    

  
}

int Globals::globalAssign( event* the_event ){
  char errorMessage[65535]; 
  int key;
  interface_assign_type the_type =  the_event->evt.asmt.asmt_type;
  char* lhs = the_event->evt.asmt.lhs;
  std::string keyword(lhs);

  bool result = false;

  /**@todo fix memory leak */  
  ParamMap::iterator i =parameters_.find(keyword);
  if (i != parameters_.end()) {
    if( the_type == STRING ){
       result = i->second->setData(std::string(the_event->evt.asmt.rhs.sval));
       if (!result ) {
  	    sprintf(errorMessage, "Error in parsing %s: expect %s, but get a string \"%s\".\n", keyword.c_str(), i->second->getParamType().c_str(), the_event->evt.asmt.rhs.sval);
            the_event->err_msg = strdup(errorMessage);
       }
    } else if( the_type == DOUBLE ){
      result = i->second->setData(the_event->evt.asmt.rhs.dval);
       if (!result )
         sprintf(errorMessage, "Error in parsing %s: expect %s, but get a double %f.\n", keyword.c_str(), i->second->getParamType().c_str(), the_event->evt.asmt.rhs.dval );
       the_event->err_msg = strdup(errorMessage);
    }      
    else if (the_type == INT ){
      result = i->second->setData(the_event->evt.asmt.rhs.ival);
       if (!result )
         sprintf(errorMessage,  "Error in parsing %s: expect %s, but get an int %d.\n", keyword.c_str(), i->second->getParamType().c_str(), the_event->evt.asmt.rhs.ival );
       the_event->err_msg = strdup(errorMessage);
      
    } else {
        sprintf(errorMessage,  "Internal error of parser\n");
	the_event->err_msg = strdup(errorMessage);
    }
  } else {
    sprintf(errorMessage,  "%s is an unrecognized keyword\n", keyword.c_str() );
    the_event->err_msg = strdup(errorMessage);
  }

  if (keyword == "nComponents" && getNComponents() > 0) {
    components = new Component*[getNComponents()];
  }else if (keyword == "nZconstraints" && getNZconstraints() > 0) {
    zConstraints = new ZconStamp*[getNZconstraints()];
  }
  
  return result;
}

int Globals::newComponent( event* the_event ){
  
  current_component = new Component;
  int index = the_event->evt.blk_index;
  char err[200];
  
  if( haveNComponents() && index < getNComponents() ) 
    components[index] = current_component;
  else{
    if( haveNComponents()  ){
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
  
  if( haveNZconstraints() && index < getNZconstraints() ) 
    zConstraints[index] = current_zConstraint;
  else{
    if( haveNZconstraints() ){
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

char* Globals::checkMe( void ){


  std::string err("The following required keywords are missing:\n");
  short int have_err = 0;

  ParamMap::iterator i;
  for (i = parameters_.begin(); i != parameters_.end(); ++i) {
    if (!i->second->isOptional() && i->second->empty()) {
        err +=  i->second->getKeyword() + "\n";
    }
  }

  CheckParameter(ForceField, isNotEmpty());
  CheckParameter(NComponents,isPositive()); 
  CheckParameter(TargetTemp, isPositive());
  CheckParameter(Ensemble, isEqualIgnoreCase(std::string("NVE")) || 
		 isEqualIgnoreCase(std::string("NVT")) ||
		 isEqualIgnoreCase(std::string("NPTi")) || 
		 isEqualIgnoreCase(std::string("NPTf")) || 
		 isEqualIgnoreCase(std::string("NPTxyz")) );
  CheckParameter(Dt, isPositive());
  CheckParameter(RunTime, isPositive());
  CheckParameter(InitialConfig, isNotEmpty());
  CheckParameter(FinalConfig, isNotEmpty());
  CheckParameter(NMol, isPositive());
  CheckParameter(Density, isPositive());
  CheckParameter(Box, isPositive());
  CheckParameter(BoxX, isPositive());
  CheckParameter(BoxY, isPositive());
  CheckParameter(BoxZ, isPositive());
  CheckParameter(SampleTime, isNonNegative());
  CheckParameter(ResetTime, isNonNegative());
  CheckParameter(StatusTime, isNonNegative());
  CheckParameter(CutoffRadius, isPositive());
  CheckParameter(SwitchingRadius, isNonNegative());
  CheckParameter(Dielectric, isPositive());
  CheckParameter(ThermalTime,  isNonNegative());
  CheckParameter(TargetPressure,  isPositive());
  CheckParameter(TauThermostat, isPositive());
  CheckParameter(TauBarostat, isPositive());
  CheckParameter(ZconsTime, isPositive());
  CheckParameter(NZconstraints, isPositive());  
  CheckParameter(ZconsTol, isPositive());
  //CheckParameter(ZconsForcePolicy,);
  CheckParameter(Seed, isPositive());
  CheckParameter(Minimizer, isEqualIgnoreCase(std::string("SD")) || 
		 isEqualIgnoreCase(std::string("CG")));
  CheckParameter(MinimizerMaxIter, isPositive());
  CheckParameter(MinimizerWriteFrq, isPositive());
  CheckParameter(MinimizerStepSize, isPositive());
  CheckParameter(MinimizerFTol, isPositive());
  CheckParameter(MinimizerGTol, isPositive());
  CheckParameter(MinimizerLSTol, isPositive());
  CheckParameter(MinimizerLSMaxIter, isPositive());
  CheckParameter(ZconsGap, isPositive());
  CheckParameter(ZconsFixtime, isPositive());
  CheckParameter(ThermodynamicIntegrationLambda, isPositive());
  CheckParameter(ThermodynamicIntegrationK, isPositive());
  CheckParameter(ForceFieldVariant, isNotEmpty());
  CheckParameter(ForceFieldFileName, isNotEmpty());
  CheckParameter(ThermIntDistSpringConst, isPositive());
  CheckParameter(ThermIntThetaSpringConst, isPositive());
  CheckParameter(ThermIntOmegaSpringConst, isPositive());
  CheckParameter(SurfaceTension, isPositive());
  CheckParameter(ElectrostaticSummationMethod, 
		 isEqualIgnoreCase(std::string("NONE")) || 
		 isEqualIgnoreCase(std::string("SHIFTED_POTENTIAL")) || 
		 isEqualIgnoreCase(std::string("SHIFTED_FORCE")) || 
		 isEqualIgnoreCase(std::string("REACTION_FIELD")));
  CheckParameter(ElectrostaticScreeningMethod, 
		 isEqualIgnoreCase(std::string("UNDAMPED")) || 
		 isEqualIgnoreCase(std::string("DAMPED"))); 
  CheckParameter(CutoffPolicy, isEqualIgnoreCase(std::string("MIX")) || 
		 isEqualIgnoreCase(std::string("MAX")) || 
		 isEqualIgnoreCase(std::string("TRADITIONAL")));
  CheckParameter(SwitchingFunctionType, 
		 isEqualIgnoreCase(std::string("CUBIC")) || 
		 isEqualIgnoreCase(std::string("FIFTH_ORDER_POLYNOMIAL")));
  //CheckParameter(StatFileFormat,);     
  //CheckParameter(MixingRule,);
  CheckParameter(OrthoBoxTolerance, isPositive());  
  CheckParameter(ThermIntDistSpringConst, isPositive());
  CheckParameter(ThermIntThetaSpringConst, isPositive());
  CheckParameter(ThermIntOmegaSpringConst, isPositive());
  CheckParameter(DampingAlpha,isNonNegative());
  CheckParameter(SkinThickness, isPositive());
  
  //@todo memory leak
  if( have_err )
    return strdup( err.c_str() );
  
  return NULL;


}

int Globals::globalEnd( event* the_event ){
  
  the_event->err_msg = checkMe();
  if( the_event->err_msg != NULL ) return 0;

  return 1;
}
