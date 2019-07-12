/**********************************************************************

This basic Periodic Table class was originally taken from the data.cpp
file in OpenBabel. The code has been modified to match the OpenMD coding style.

We have retained the OpenBabel copyright and GPL license on this class:  

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

/**
 * @file ElementsTable.cpp
 * @author gezelter
 * @date 12/21/2007
 * @version 1.0
 */

#include "config.h"

#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "utils/ElementsTable.hpp"
#include "utils/simError.h"
#include "io/ifstrstream.hpp"

#ifdef WIN32
#define FILE_SEP_CHAR "\\"
#else
#define FILE_SEP_CHAR "/"
#endif

#ifndef BUFF_SIZE
#define BUFF_SIZE 32768
#endif

namespace OpenMD {

  ElementsTable etab;

  ElementsTable::ElementsTable() {
    init_ = false;
    dir_ = std::string("TO_STRING(FRC_PATH)");
    envvar_ = "FORCE_PARAM_PATH";
    filename_ = "element.txt";
  }
  
  ElementsTable::~ElementsTable() {
    std::vector<Element*>::iterator i;
    for (i = elements_.begin(); i != elements_.end(); ++i)
      delete *i;
  }
  
  void ElementsTable::ParseLine(const char *buffer) {
    int num, maxbonds;
    char symbol[4];
    char name[256];
    RealType Rcov,Rvdw,mass, elNeg, ARENeg, ionize, elAffin;
    RealType red, green, blue;
    
    // skip comment line (at the top)
    if (buffer[0] != '#')  {
      sscanf(buffer,"%d %3s %lf %lf %*f %lf %d %lf %lf %lf %lf %lf %lf %lf %255s",
             &num,
             symbol,
             &ARENeg,
             &Rcov,
             &Rvdw,
             &maxbonds,
             &mass,
             &elNeg,
             &ionize,
             &elAffin,
             &red,
             &green,
             &blue,
             name);
      
      Element *ele = new Element(num, symbol, ARENeg, Rcov, Rvdw, maxbonds,
                                 mass, elNeg, ionize, elAffin,
                                 red, green, blue, name);
      elements_.push_back(ele);

    }
  }

  unsigned int ElementsTable::GetNumberOfElements() { 
    if (!init_)
      Init();
    
    return elements_.size();
  }

  const char *ElementsTable::GetSymbol(int atomicnum) {
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum >= static_cast<int>(elements_.size()))
      return("\0");
    
    return(elements_[atomicnum]->GetSymbol());
  }

  int ElementsTable::GetMaxBonds(int atomicnum) {
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum >= static_cast<int>(elements_.size()))
      return(0);
    
    return(elements_[atomicnum]->GetMaxBonds());
  }
  
  RealType ElementsTable::GetElectroNeg(int atomicnum) {
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum >= static_cast<int>(elements_.size()))
      return(0.0);
    
    return(elements_[atomicnum]->GetElectroNeg());
  }

  RealType ElementsTable::GetAllredRochowElectroNeg(int atomicnum) {
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum >= static_cast<int>(elements_.size()))
      return(0.0);

    return(elements_[atomicnum]->GetAllredRochowElectroNeg());
  }
  
  RealType ElementsTable::GetIonization(int atomicnum) {
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum >= static_cast<int>(elements_.size()))
      return(0.0);
    
    return(elements_[atomicnum]->GetIonization());
  }
  

  RealType ElementsTable::GetElectronAffinity(int atomicnum) {
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum >= static_cast<int>(elements_.size()))
      return(0.0);

    return(elements_[atomicnum]->GetElectronAffinity());
  }
  
  std::vector<RealType> ElementsTable::GetRGB(int atomicnum) {
    if (!init_)
      Init();
    
    std::vector <RealType> colors;
    colors.reserve(3);
    
    if (atomicnum < 0 || atomicnum >= static_cast<int>(elements_.size())) {
      colors.push_back(0.0);
      colors.push_back(0.0);
      colors.push_back(0.0);
      return(colors);
    }

    colors.push_back(elements_[atomicnum]->GetRed());
    colors.push_back(elements_[atomicnum]->GetGreen());
    colors.push_back(elements_[atomicnum]->GetBlue());
    
    return (colors);
  }
  
  std::string ElementsTable::GetName(int atomicnum) {
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum >= static_cast<int>(elements_.size()))
      return("Unknown");
    
    return(elements_[atomicnum]->GetName());
  }

  RealType ElementsTable::GetVdwRad(int atomicnum) {
    if (!init_)
      Init();

    if (atomicnum < 0 || atomicnum >= static_cast<int>(elements_.size()))
      return(0.0);

    return(elements_[atomicnum]->GetVdwRad());
  }
  
  RealType ElementsTable::CorrectedBondRad(int atomicnum, int hyb) {
    RealType rad;
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum >= static_cast<int>(elements_.size()))
      return(1.0);
    
    rad = elements_[atomicnum]->GetCovalentRad();
    
    if (hyb == 2)
      rad *= 0.95;
    else if (hyb == 1)
      rad *= 0.90;

    return(rad);
  }

  RealType ElementsTable::CorrectedVdwRad(int atomicnum, int hyb) {
    RealType rad;
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum >= static_cast<int>(elements_.size()))
      return(1.95);
    
    rad = elements_[atomicnum]->GetVdwRad();
    
    if (hyb == 2)
      rad *= 0.95;
    else if (hyb == 1)
      rad *= 0.90;
    
    return(rad);
  }

  RealType ElementsTable::GetCovalentRad(int atomicnum) {
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum >= static_cast<int>(elements_.size()))
      return(0.0);
    
    return(elements_[atomicnum]->GetCovalentRad());
  }
  
  RealType ElementsTable::GetMass(int atomicnum) {
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum >= static_cast<int>(elements_.size()))
      return(0.0);
    
    return(elements_[atomicnum]->GetMass());
  }
  
  int ElementsTable::GetAtomicNum(const char *sym) {
    int temp;
    return GetAtomicNum(sym, temp);
  }

  int ElementsTable::GetAtomicNum(const char *identifier, int &iso) {
    if (!init_)
      Init();
    
    // Compare to symbol    
    std::vector<Element*>::iterator i;
    for (i = elements_.begin();i != elements_.end(); ++i)
      if (!strncasecmp(identifier,(*i)->GetSymbol(),3))
        return((*i)->GetAtomicNum());

    // Compare to IUPAC name (an abbreviated name will also work if 5
    // letters or more)
    int numCharsToTest = std::max<int>(strlen(identifier), 5);
    for (i = elements_.begin();i != elements_.end();++i)
      if (strncasecmp(identifier,(*i)->GetName().c_str(),numCharsToTest) == 0)
        return((*i)->GetAtomicNum());
    
    if (strcasecmp(identifier, "D") == 0 ||
        (strcasecmp(identifier, "Deuterium") == 0) ) {
      iso = 2;
      return(1);
    } else if (strcasecmp(identifier, "T") == 0 ||
               (strcasecmp(identifier, "Tritium") == 0) ) {
      iso = 3;
      return(1);
    } else if (strcasecmp(identifier, "Hl") == 0) {
      // ligand hydrogen -- found in some CIF PR#3048959.
      sprintf( painCave.errMsg,
               "ElementsTable warning.\n"
               "\tCannot understand the element label %s\n"
               "\tGuessing it's hydrogen\n", identifier);
      painCave.isFatal = 0;
      painCave.severity = OPENMD_WARNING;
      simError();            
      return(1);
    } else
      iso = 0;

    if(identifier[0]!='*') {
      // Quiet down the error messages for now:
      //sprintf( painCave.errMsg,
      //         "ElementsTable warning.\n"
      //         "\tCannot understand the element label %s\n", identifier);
      //painCave.isFatal = 0;
      //painCave.severity = OPENMD_WARNING;
      //simError();            
    }
    return(0);
  }

  int ElementsTable::GetAtomicNum(std::string name, int &iso) {
    return GetAtomicNum(name.c_str(), iso);
  }


  void ElementsTable::Init() {
    if (init_)
      return;
    init_ = true;
    
    std::string buffer, subbuffer;
    ifstrstream ifs;
    char charBuffer[BUFF_SIZE];

    // First, look for an environment variable
    if (getenv(envvar_.c_str()) != NULL) {
      buffer = getenv(envvar_.c_str());
      buffer += FILE_SEP_CHAR;
     
      if (!subdir_.empty()) {
        subbuffer = buffer;
        subbuffer += subdir_;
        subbuffer += FILE_SEP_CHAR;
      }
            
      buffer += filename_;
      subbuffer += filename_;
    }
    
    ifs.clear();
    ifs.open(subbuffer.c_str());
    
    if( !(&ifs)->is_open() ) {
      
      ifs.clear();
      ifs.open(buffer.c_str());
      
      if ( !(&ifs)->is_open() ) {
        sprintf( painCave.errMsg,
                 "ElementsTable error.\n"
                 "\tunable to open datafile %s \n", filename_.c_str());
        painCave.isFatal = 0;
        simError();
      }
    }
    
    if (ifs) {      
      while(ifs.getline(charBuffer,BUFF_SIZE))
        ParseLine(charBuffer);
    }
    
    if (ifs)
      ifs.close();
    
    if (GetSize() == 0) {
      sprintf( painCave.errMsg,
               "ElementsTable error.\n"
               "\tCannot initialize database %s \n", filename_.c_str());
      painCave.isFatal = 0;
      simError();       
    } 
  }
}
