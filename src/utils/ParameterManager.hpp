/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#ifndef UTILS_PARAMETERMANAGER_HPP
#define UTILS_PARAMETERMANAGER_HPP

#include <config.h>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "utils/CaseConversion.hpp"
#include "utils/StringTokenizer.hpp"
#include "utils/simError.h"

template<typename T>
struct ParameterTraits;

// string
template<>
struct ParameterTraits<std::string> {
  using RepType = std::string;  // Representation type of the value

  template<typename T>
  static bool convert(T, RepType&) {
    return false;
  }  // !NB everything is ok

  template<typename T>
  static RepType convert(T v) {
    RepType tmp;
    convert(v, tmp);
    return tmp;
  }

  static bool convert(RepType v, RepType& r) {
    r = v;
    return true;
  }

  static std::string getParamType() { return "string"; }
};
// bool
template<>
struct ParameterTraits<bool> {
  using RepType = bool;

  template<typename T>
  static bool convert(T, RepType&) {
    return false;
  }

  template<typename T>
  static RepType convert(T v) {
    RepType tmp;
    convert(v, tmp);
    return tmp;
  }

  static bool convert(std::string v, RepType& r) {
    OpenMD::toLower(v);
    bool result = false;
    if (v == "true") {
      r      = true;
      result = true;
    } else if (v == "false") {
      r      = false;
      result = true;
    }

    return result;
  }

  static std::string getParamType() { return "bool"; }
};

// int
template<>
struct ParameterTraits<int> {
  using RepType = int;

  template<typename T>
  static bool convert(T, RepType&) {
    return false;
  }

  template<typename T>
  static RepType convert(T v) {
    RepType tmp;
    convert(v, tmp);
    return tmp;
  }

  static bool convert(RepType v, RepType& r) {
    r = v;
    return true;
  }

  static std::string getParamType() { return "int"; }
};

// int
template<>
struct ParameterTraits<unsigned long int> {
  using RepType = unsigned long int;

  template<typename T>
  static bool convert(T, RepType&) {
    return false;
  }

  template<typename T>
  static RepType convert(T v) {
    RepType tmp;
    convert(v, tmp);
    return tmp;
  }

  static bool convert(RepType v, RepType& r) {
    r = v;
    return true;
  }

  static bool convert(int v, RepType& r) {
    r = static_cast<unsigned long int>(v);
    return true;
  }

  static std::string getParamType() { return "unsigned long int"; }
};

// RealType
template<>
struct ParameterTraits<RealType> {
  using RepType = RealType;

  template<typename T>
  static bool convert(T, RepType&) {
    return false;
  }

  template<typename T>
  static RepType convert(T v) {
    RepType tmp;
    convert(v, tmp);
    return tmp;
  }

  static bool convert(RepType v, RepType& r) {
    r = v;
    return true;
  }

  static bool convert(int v, RepType& r) {
    r = static_cast<RealType>(v);
    return true;
  }

  static bool convert(unsigned long int v, RepType& r) {
    r = static_cast<RealType>(v);
    return true;
  }

  static std::string getParamType() { return "RealType"; }
};

// Pair of ints
template<>
struct ParameterTraits<std::pair<int, int>> {
  using RepType = std::pair<int, int>;

  template<typename T>
  static bool convert(T, RepType&) {
    return false;
  }

  template<typename T>
  static RepType convert(T v) {
    RepType tmp;
    convert(v, tmp);
    return tmp;
  }

  static bool convert(RepType v, RepType& r) {
    r = v;
    return true;
  }

  static bool convert(std::string v, RepType& r) {
    OpenMD::StringTokenizer tokenizer(v, " ;,\t\n\r");
    if (tokenizer.countTokens() == 2) {
      int atom1 = tokenizer.nextTokenAsInt();
      int atom2 = tokenizer.nextTokenAsInt();
      r         = std::make_pair(atom1, atom2);
      return true;
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "ParameterManager Error: "
               "Incorrect number of tokens to make a pair!\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
    return false;
  }

  static std::string getParamType() { return "std::pair<int, int>"; }
};

template<>
struct ParameterTraits<std::vector<RealType>> {
  using RepType = std::vector<RealType>;

  template<typename T>
  static bool convert(T, RepType&) {
    return false;
  }

  template<typename T>
  static RepType convert(T v) {
    RepType tmp;
    convert(v, tmp);
    return tmp;
  }

  static bool convert(RepType v, RepType& r) {
    r = v;
    return true;
  }

  static bool convert(std::string v, RepType& r) {
    std::cerr << "calling tokenizer\n";
    OpenMD::StringTokenizer tokenizer(v, " ();,\t\n\r");
    unsigned int size = tokenizer.countTokens();
    r                 = std::vector<RealType>(size, 0.0);
    for (unsigned int i = 0; i < size; i++) {
      RealType v = tokenizer.nextTokenAsDouble();
      r[i]       = v;
    }
    return true;
  }

  static std::string getParamType() { return "std::vector<RealType>"; }
};

class ParameterBase {
public:
  ParameterBase() :
      keyword_(), optional_(false), defaultValue_(false), empty_(true) {}
  virtual ~ParameterBase() {}
  bool isOptional() { return optional_; }
  void setOptional(bool optional) { optional_ = optional; }
  bool hasDefaultValue() { return defaultValue_; }
  virtual bool isValid() { return true; }
  const std::string& getKeyword() { return keyword_; }
  void setKeyword(const std::string& keyword) { keyword_ = keyword; }
  bool empty() { return empty_; }
  virtual bool setData(std::string)           = 0;
  virtual bool setData(int)                   = 0;
  virtual bool setData(unsigned long int)     = 0;
  virtual bool setData(RealType)              = 0;
  virtual bool setData(std::pair<int, int>)   = 0;
  virtual bool setData(std::vector<RealType>) = 0;
  virtual std::string getParamType()          = 0;

protected:
  std::string keyword_;
  bool optional_;
  bool defaultValue_;
  bool empty_;
};

template<class ParamType>
class Parameter : public ParameterBase {
public:
  using ValueType = ParameterTraits<ParamType>;

  void setDefaultValue(const ParamType& value) {
    data_         = value;
    defaultValue_ = true;
    empty_        = false;
  }
  ParamType getData() { return data_; }

  virtual bool setData(std::string sval) {
    return internalSetData<std::string>(sval);
  }
  virtual bool setData(int ival) { return internalSetData<int>(ival); }
  virtual bool setData(unsigned long int lival) {
    return internalSetData<unsigned long int>(lival);
  }
  virtual bool setData(RealType dval) {
    return internalSetData<RealType>(dval);
  }
  virtual bool setData(std::pair<int, int> pval) {
    return internalSetData<std::pair<int, int>>(pval);
  }
  virtual bool setData(std::vector<RealType> pval) {
    return internalSetData<std::vector<RealType>>(pval);
  }

  virtual std::string getParamType() {
    return ParameterTraits<ParamType>::getParamType();
  }

private:
  template<class T>
  bool internalSetData(T data) {
    ParamType tmp;
    bool result = ValueType::convert(data, tmp);
    if (result) {
      empty_ = false;
      data_  = tmp;
    }
    return result;
  }

private:
  ParamType data_;
};

#define DeclareParameter(NAME, TYPE)          \
private:                                      \
  Parameter<TYPE> NAME;                       \
                                              \
public:                                       \
  bool have##NAME() { return !NAME.empty(); } \
  TYPE get##NAME() { return NAME.getData(); }

#define DeclareAlterableParameter(NAME, TYPE) \
private:                                      \
  Parameter<TYPE> NAME;                       \
                                              \
public:                                       \
  bool have##NAME() { return !NAME.empty(); } \
  TYPE get##NAME() { return NAME.getData(); } \
  bool set##NAME(TYPE s) { return NAME.setData(s); }

#define DefineParameter(NAME, KEYWORD)                                  \
  NAME.setKeyword(KEYWORD);                                             \
  parameters_.insert(std::map<std::string, ParameterBase*>::value_type( \
      std::string(KEYWORD), static_cast<ParameterBase*>(&NAME)));

#define DefineOptionalParameter(NAME, KEYWORD)                          \
  NAME.setKeyword(KEYWORD);                                             \
  NAME.setOptional(true);                                               \
  parameters_.insert(std::map<std::string, ParameterBase*>::value_type( \
      std::string(KEYWORD), static_cast<ParameterBase*>(&NAME)));

#define DefineOptionalParameterWithDefaultValue(NAME, KEYWORD, DEFAULTVALUE) \
  NAME.setKeyword(KEYWORD);                                                  \
  NAME.setOptional(true);                                                    \
  NAME.setDefaultValue(DEFAULTVALUE);                                        \
  parameters_.insert(std::map<std::string, ParameterBase*>::value_type(      \
      std::string(KEYWORD), static_cast<ParameterBase*>(&NAME)));

#define CheckParameter(NAME, CONSTRAINT)                         \
  if (!NAME.empty()) {                                           \
    if (!(CONSTRAINT)(NAME.getData())) {                         \
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,        \
               "Error in checking %s : should be %s\n",          \
               NAME.getKeyword().c_str(),                        \
               (CONSTRAINT).getConstraintDescription().c_str()); \
      painCave.isFatal  = 1;                                     \
      painCave.severity = OPENMD_ERROR;                          \
      simError();                                                \
    }                                                            \
  }

#endif
