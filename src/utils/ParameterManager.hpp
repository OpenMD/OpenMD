/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
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
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *                                                                      
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [4] , Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). *
 *  ParameterManager.hpp
 *
 *  Created by Charles F. Vardeman II on 11/16/05.
 *  @author  Charles F. Vardeman II 
 *  @version $Id$
 *
 */

#ifndef UTILS_PARAMETERMANAGER_HPP
#define UTILS_PARAMETERMANAGER_HPP

#include <iostream>
#include <cstdio>

#include <stdlib.h>
#include <vector>
#include <string>
#include <map>
#include "config.h"


#include "utils/simError.h"
#include "utils/StringTokenizer.hpp"
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
    OpenMD::toLower(v); 
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

//int   
template<>
struct ParameterTraits<unsigned long int>{
  typedef unsigned long int RepType;
  template<typename T> static bool    convert(T, RepType&){return false;} 
  template<typename T> static RepType convert(T v)        {RepType tmp; convert(v,tmp);return tmp;} 
  static bool convert(RepType v, RepType& r)            { r=v; return true;}
  static bool convert(int v, RepType& r)                {r = static_cast<unsigned long int>(v); return true;}
  static std::string getParamType() { return "unsigned long int";}  
};

//RealType
template<>                     
struct ParameterTraits<RealType>{
  typedef RealType RepType;
  template<typename T> static bool    convert(T, RepType&){return false;} 
  template<typename T> static RepType convert(T v)        {RepType tmp; convert(v,tmp);return tmp;} 
  static bool convert(RepType v, RepType& r)            {r=v; return true;}
  static bool convert(int v, RepType& r)                {r = static_cast<RealType>(v); return true;}
  static bool convert(unsigned long int v, RepType& r)  {r = static_cast<RealType>(v); return true;}
  static std::string getParamType() { return "RealType";}    
};

//Pair of ints
template<>                     
struct ParameterTraits<std::pair<int, int> >{
  typedef std::pair<int, int>  RepType;
  template<typename T> static bool    convert(T, RepType&){return false;} 
  template<typename T> static RepType convert(T v)        {RepType tmp; convert(v,tmp);return tmp;} 
  static bool convert(RepType v, RepType& r)            {r=v; return true;}
  static bool convert(std::string v, RepType& r) { 
    OpenMD::StringTokenizer tokenizer(v," ;,\t\n\r");
    if (tokenizer.countTokens() == 2) {
      int atom1 = tokenizer.nextTokenAsInt();
      int atom2 = tokenizer.nextTokenAsInt();
      r = std::make_pair(atom1, atom2);
      return true;
    } else {
      sprintf(painCave.errMsg, 
              "ParameterManager Error: "
              "Not enough tokens to make pair!\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();    
    }
    return false;
  }
  static std::string getParamType() { return "std::pair<int, int>";}    
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
  virtual bool setData(unsigned long int) = 0;
  virtual bool setData(RealType) = 0;
  virtual bool setData(std::pair<int, int>) = 0;
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

  virtual bool setData(unsigned long int lival) {
    return internalSetData<unsigned long int>(lival);
  }
  
  virtual bool setData(RealType dval) {
    return internalSetData<RealType>(dval);
  }

  virtual bool setData(std::pair<int, int> pval) {
    return internalSetData<std::pair<int, int> >(pval);
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

#define DeclareAlterableParameter(NAME, TYPE)         \
private:                                                   \
Parameter<TYPE> NAME;                                     \
public:                                                      \
bool have##NAME() { return !NAME.empty();}  \
TYPE get##NAME() { return NAME.getData();}  \
bool set##NAME(TYPE s) { return NAME.setData(s);}  \



#define DefineParameter(NAME,KEYWORD)                              \
NAME.setKeyword(KEYWORD);                  \
parameters_.insert(std::map<std::string, ParameterBase*>::value_type(std::string(KEYWORD), static_cast<ParameterBase*>(&NAME)));

#define DefineOptionalParameter(NAME,KEYWORD)                              \
NAME.setKeyword(KEYWORD); NAME.setOptional(true);                    \
parameters_.insert(std::map<std::string, ParameterBase*>::value_type(std::string(KEYWORD), static_cast<ParameterBase*>(&NAME)));

#define DefineOptionalParameterWithDefaultValue(NAME,KEYWORD, DEFAULTVALUE)                              \
NAME.setKeyword(KEYWORD); NAME.setOptional(true); NAME.setDefaultValue(DEFAULTVALUE);                      \
parameters_.insert(std::map<std::string, ParameterBase*>::value_type(std::string(KEYWORD), static_cast<ParameterBase*>(&NAME)));

#define CheckParameter(NAME, CONSTRAINT)                              \
if (!NAME.empty()) { if (!(CONSTRAINT)(NAME.getData())) { sprintf(painCave.errMsg,"Error in checking %s : should be %s\n",NAME.getKeyword().c_str(),(CONSTRAINT).getConstraintDescription().c_str()); painCave.isFatal = 1; painCave.severity = OPENMD_ERROR; simError();} }                 



#endif
