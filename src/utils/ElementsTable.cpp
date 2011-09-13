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
 * @time 11:30am
 * @version 1.0
 */

#include <iostream>
#include "config.h"
#include <cstdlib>
#include <string>
#include <fstream>
#include <cstdlib>
#include "utils/ElementsTable.hpp"
#include "utils/simError.h"
#include "io/ifstrstream.hpp"

#if !HAVE_STRNCASECMP
extern "C" int strncasecmp(const char *s1, const char *s2, size_t n);
#endif

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
    for (i = elements_.begin(); i != elements_.end(); i++)
      delete *i;
  }
  
  void ElementsTable::ParseLine(const char *line) {
    int num, maxbonds;
    char symbol[5];
    char name[256];
    RealType Rcov,Rvdw,mass, elNeg, ionize, elAffin;
    RealType red, green, blue;

    // skip comment line (at the top)
    if (line[0] != '#')  {
      sscanf(line,"%d %5s %lf %*f %lf %d %lf %lf %lf %lf %lf %lf %lf %255s",
             &num,
             symbol,
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
    
      Element *ele = new Element(num, symbol, Rcov, Rvdw, maxbonds, mass, 
                                 elNeg, ionize, elAffin, red, green, blue, 
                                 name);
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
    
    if (atomicnum < 0 || atomicnum > static_cast<int>(elements_.size()))
      return("\0");
    
    return(elements_[atomicnum]->GetSymbol());
  }

  int ElementsTable::GetMaxBonds(int atomicnum) {
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum > static_cast<int>(elements_.size()))
      return(0);
    
    return(elements_[atomicnum]->GetMaxBonds());
  }
  
  RealType ElementsTable::GetElectroNeg(int atomicnum) {
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum > static_cast<int>(elements_.size()))
      return(0.0);
    
    return(elements_[atomicnum]->GetElectroNeg());
  }
  
  RealType ElementsTable::GetIonization(int atomicnum) {
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum > static_cast<int>(elements_.size()))
      return(0.0);
    
    return(elements_[atomicnum]->GetIonization());
  }
  

  RealType ElementsTable::GetElectronAffinity(int atomicnum) {
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum > static_cast<int>(elements_.size()))
      return(0.0);

    return(elements_[atomicnum]->GetElectronAffinity());
  }
  
  std::vector<RealType> ElementsTable::GetRGB(int atomicnum) {
    if (!init_)
      Init();
    
    std::vector <RealType> colors;
    colors.reserve(3);
    
    if (atomicnum < 0 || atomicnum > static_cast<int>(elements_.size())) {
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
    
    if (atomicnum < 0 || atomicnum > static_cast<int>(elements_.size()))
      return("Unknown");
    
    return(elements_[atomicnum]->GetName());
  }

  RealType ElementsTable::GetVdwRad(int atomicnum) {
    if (!init_)
      Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(elements_.size()))
      return(0.0);

    return(elements_[atomicnum]->GetVdwRad());
  }
  
  RealType ElementsTable::CorrectedBondRad(int atomicnum, int hyb) {
    RealType rad;
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum > static_cast<int>(elements_.size()))
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
    
    if (atomicnum < 0 || atomicnum > static_cast<int>(elements_.size()))
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
    
    if (atomicnum < 0 || atomicnum > static_cast<int>(elements_.size()))
      return(0.0);
    
    return(elements_[atomicnum]->GetCovalentRad());
  }
  
  RealType ElementsTable::GetMass(int atomicnum) {
    if (!init_)
      Init();
    
    if (atomicnum < 0 || atomicnum > static_cast<int>(elements_.size()))
      return(0.0);
    
    return(elements_[atomicnum]->GetMass());
  }
  
  int ElementsTable::GetAtomicNum(const char *sym) {
    int temp;
    return GetAtomicNum(sym, temp);
  }

  int ElementsTable::GetAtomicNum(const char *sym, int &iso) {
    if (!init_)
      Init();
    
    std::vector<Element*>::iterator i;
    for (i = elements_.begin();i != elements_.end();i++)
      if (!strncasecmp(sym,(*i)->GetSymbol(),2))
        return((*i)->GetAtomicNum());

    if (strcasecmp(sym, "D") == 0) {
      iso = 2;
      return(1);
    } else if (strcasecmp(sym, "T") == 0) {
      iso = 3;
      return(1);
    } else 
      iso = 0;
    return(0);
  }

  void ElementsTable::Init() {
    if (init_)
      return;
    init_ = true;
    
    std::string buffer, subbuffer;
    ifstrstream ifs1, ifs2, ifs3, ifs4, *ifsP;
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

      
      ifs1.open(subbuffer.c_str());
      ifsP= &ifs1;
      if (!(ifsP->is_open())) {
        ifs2.open(buffer.c_str());
        ifsP = &ifs2;
      }
      
    } else {
      sprintf( painCave.errMsg,
               "ElementsTable error.\n"
               "\tunable to open datafile %s \n", filename_.c_str());
      painCave.isFatal = 0;
      simError();
    }
      
    char charBuffer[BUFF_SIZE];
    if ((*ifsP)) {
      while(ifsP->getline(charBuffer,BUFF_SIZE))
        ParseLine(charBuffer);

      if (ifs1)
	ifs1.close();
      if (ifs2)
	ifs2.close();
      if (ifs3)
	ifs3.close();
      if (ifs4)
	ifs4.close();
      
      if (GetSize() == 0) {
	sprintf( painCave.errMsg,
		 "ElementsTable error.\n"
		 "\tCannot initialize database %s \n", filename_.c_str());
	painCave.isFatal = 0;
	simError();
      }
    
    }
  
  }
}
