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
 
#ifndef IO_GLOBALS_HPP
#define IO_GLOBALS_HPP

#include <iostream>

#include <stdlib.h>
#include <vector>
#include <string>
#include <map>

#include "io/BASS_interface.h"
#include "types/Component.hpp"
#include "types/MakeStamps.hpp"
#include "types/ZconStamp.hpp"
#include "utils/CaseConversion.hpp"

template<typename T> 
struct ParameterTraits;

//string
template<>                     
struct ParameterTraits<std::string>{
  typedef std::string RepType;       // Representation type of the value

template<typename T> static bool    convert(T v, RepType& r){return false;} // !NB everything is ok
template<typename T> static RepType convert(T v)            {RepType tmp; convert(v,tmp);return tmp;} 
  static bool convert(RepType v, RepType& r) { r = v; return true;}
  static std::string getParamType() { return "string";}  
}; 
//bool
template<>                     
struct ParameterTraits<bool>{
  typedef bool RepType;
  template<typename T> static bool    convert(T, RepType&){return false;} 
  template<typename T> static RepType convert(T v)        {RepType tmp; convert(v,tmp);return tmp;} 
  static bool convert(std::string v, RepType& r) { 
    oopse::toLower(v); 
    bool result = false;
    if (v == "true") {
        r = true;
        result = true;
    } else if (v == "false") {
        r = false;
        result = true;
    }
    
     return result;
  }
  static std::string getParamType() { return "bool";}  
};

//int   
template<>
struct ParameterTraits<int>{
    typedef int RepType;
    template<typename T> static bool    convert(T, RepType&){return false;} 
    template<typename T> static RepType convert(T v)        {RepType tmp; convert(v,tmp);return tmp;} 
    static bool convert(RepType v, RepType& r)            { r=v; return true;}
    static std::string getParamType() { return "int";}  
};

//double
template<>                     
struct ParameterTraits<double>{
    typedef double RepType;
    template<typename T> static bool    convert(T, RepType&){return false;} 
    template<typename T> static RepType convert(T v)        {RepType tmp; convert(v,tmp);return tmp;} 
    static bool convert(RepType v, RepType& r)            {r=v; return true;}
    static bool convert(int v, RepType& r)                {r = static_cast<double>(v); return true;}
    static std::string getParamType() { return "double";}    
};


class ParameterBase {
  public:    
    ParameterBase() : keyword_(), optional_(false), defaultValue_(false), empty_(true) {}
    virtual ~ParameterBase() {}
    bool isOptional() {return optional_;}
    void setOptional(bool optional) {optional_ = optional;}
    bool hasDefaultValue() {return defaultValue_;}
    virtual bool isValid() { return true;}
    const std::string& getKeyword() {return keyword_;}
    void setKeyword(const std::string& keyword) { keyword_ = keyword;}
    bool empty() {return empty_;}
    virtual bool setData(std::string) = 0;
    virtual bool setData(int) = 0;
    virtual bool setData(double) = 0;
    virtual std::string getParamType() = 0;
  protected:
    std::string keyword_;
    bool optional_;
    bool defaultValue_;
    bool empty_;
};

template<class ParamType>
class Parameter : public ParameterBase{
   public:    
    typedef ParameterTraits<ParamType> ValueType;
     void setDefaultValue(const ParamType& value) {data_ = value; defaultValue_ = true;}
     ParamType getData() { return data_;}

    virtual bool setData(std::string sval) {
        return internalSetData<std::string>(sval);
    }
    virtual bool setData(int ival) {
        return internalSetData<int>(ival);
    }
    
    virtual bool setData(double dval) {
        return internalSetData<double>(dval);
    }

    virtual std::string getParamType() { return ParameterTraits<ParamType>::getParamType();}
   private: 
     template<class T> bool internalSetData(T data) {
        ParamType tmp;
        bool result = ValueType::convert(data, tmp);
        if (result) {
            empty_ = false;
            data_ = tmp;
        }
        return result;
     }
     
   private:
     ParamType data_;
      
};

#define DeclareParameter(NAME, TYPE)         \
    private:                                                   \
      Parameter<TYPE> NAME;                                     \
    public:                                                      \
      bool have##NAME() { return !NAME.empty();}  \
      TYPE get##NAME() { return NAME.getData();}


class Globals {
  public:
    Globals();
    
  DeclareParameter(ForceField, std::string);
  DeclareParameter(NComponents, int);  
  DeclareParameter(TargetTemp, double);
  DeclareParameter(Ensemble, std::string);
  DeclareParameter(Dt, double);
  DeclareParameter(RunTime, double);
  DeclareParameter(InitialConfig, std::string);
  DeclareParameter(FinalConfig, std::string);
  DeclareParameter(NMol, int);
  DeclareParameter(Density, double);
  DeclareParameter(Box, double);
  DeclareParameter(BoxX, double);
  DeclareParameter(BoxY, double);
  DeclareParameter(BoxZ, double);
  DeclareParameter(SampleTime, double);
  DeclareParameter(ResetTime, double);
  DeclareParameter(StatusTime, double);
  DeclareParameter(CutoffRadius, double);
  DeclareParameter(SwitchingRadius, double);
  DeclareParameter(Dielectric, double);
  DeclareParameter(TempSet, bool);
  DeclareParameter(ThermalTime, double);
  DeclareParameter(MixingRule, std::string);
  DeclareParameter(UsePeriodicBoundaryConditions, bool);
  DeclareParameter(TargetPressure, double);
  DeclareParameter(TauThermostat, double);
  DeclareParameter(TauBarostat, double);
  DeclareParameter(ZconsTime, double);
  DeclareParameter(NZconstraints, int);
  DeclareParameter(ZconsTol, double);
  DeclareParameter(ZconsForcePolicy, std::string);
  DeclareParameter(Seed, int);
  DeclareParameter(UseInitalTime, bool);
  DeclareParameter(UseIntialExtendedSystemState, bool);
  DeclareParameter(OrthoBoxTolerance, double);
  DeclareParameter(Minimizer, std::string);
  DeclareParameter(MinimizerMaxIter, double);
  DeclareParameter(MinimizerWriteFrq, int);
  DeclareParameter(MinimizerStepSize, double);
  DeclareParameter(MinimizerFTol, double);
  DeclareParameter(MinimizerGTol, double);
  DeclareParameter(MinimizerLSTol, double);
  DeclareParameter(MinimizerLSMaxIter, int);
  DeclareParameter(ZconsGap, double);
  DeclareParameter(ZconsFixtime, double);
  DeclareParameter(ZconsUsingSMD, bool);
  DeclareParameter(UseSolidThermInt, bool);
  DeclareParameter(UseLiquidThermInt, bool);
  DeclareParameter(ThermodynamicIntegrationLambda, double);
  DeclareParameter(ThermodynamicIntegrationK, double);
  DeclareParameter(ForceFieldVariant, std::string);
  DeclareParameter(ForceFieldFileName, std::string);
  DeclareParameter(ThermIntDistSpringConst, double);
  DeclareParameter(ThermIntThetaSpringConst, double);
  DeclareParameter(ThermIntOmegaSpringConst, double);
  DeclareParameter(SurfaceTension, double);
  DeclareParameter(PrintPressureTensor, bool);
  DeclareParameter(ElectrostaticSummationMethod, std::string);
  DeclareParameter(ElectrostaticScreeningMethod, std::string);
  DeclareParameter(DampingAlpha, double);
  DeclareParameter(CutoffPolicy, std::string);
  DeclareParameter(CompressDumpFile, bool);
  DeclareParameter(SkinThickness, double);
  DeclareParameter(StatFileFormat, std::string);    

  private:
    typedef std::map<std::string, ParameterBase*> ParamMap;
    ParamMap parameters_;
    
    Component* current_component;
    Component** components; // the array of components

    ZconStamp* current_zConstraint;
    ZconStamp** zConstraints; // the array of zConstraints

    char* checkMe();

  public:
    int newComponent( event* the_event );
    int componentAssign( event* the_event );
    int componentEnd( event* the_event );

    int newZconstraint( event* the_event );
    int zConstraintAssign( event* the_event );
    int zConstraintEnd( event* the_event );
  
    int globalAssign( event* the_event );
    int globalEnd( event* the_event );    

    ZconStamp** getZconStamp() {return zConstraints;}
    Component** getComponents() {return components;}
};
#endif

