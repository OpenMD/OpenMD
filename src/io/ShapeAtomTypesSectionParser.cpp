 /*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include <memory>
 
#include "brains/ForceField.hpp"
#include "io/ShapeAtomTypesSectionParser.hpp"
#include "math/RealSphericalHarmonic.hpp"
#include "math/SquareMatrix3.hpp"
#include "types/AtomType.hpp"
#include "types/DirectionalAdapter.hpp"
#include "types/ShapeAtomType.hpp"
#include "utils/StringUtils.hpp"
#include "utils/simError.h"
#include "utils/CaseConversion.hpp"

namespace OpenMD {
  
  ShapeAtomTypesSectionParser::ShapeAtomTypesSectionParser(ForceFieldOptions& options) : options_(options) {
    setSectionName("ShapeAtomTypes");
  }
  
  void ShapeAtomTypesSectionParser::parseLine(ForceField& ff, 
                                              const std::string& line, 
                                              int lineNo){
    StringTokenizer tokenizer(line);
      
    if (tokenizer.countTokens() >= 2) {
      std::string shapeTypeName = tokenizer.nextToken();
      std::string shapeFile = tokenizer.nextToken();

      AtomType* atomType = ff.getAtomType(shapeTypeName);
      if (atomType == NULL) {
        atomType = new AtomType();
        int ident = ff.getNAtomType() + 1;
        atomType->setIdent(ident); 
        atomType->setName(shapeTypeName);
        ff.addAtomType(shapeTypeName, atomType);
      } 
      
      parseShapeFile(ff, shapeFile, atomType);

    } else {
      sprintf(painCave.errMsg, 
              "ShapesAtomTypesSectionParser Error: "
              "Not enough tokens at line %d\n",
              lineNo);
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();    
    }
  }
  
  void ShapeAtomTypesSectionParser::parseShapeFile(ForceField& ff,
                                                   std::string& shapeFileName, 
                                                   AtomType* at) {
    
    const int bufferSize = 65535;
    char buffer[bufferSize];
    std::string token;

    Mat3x3d momInert;
    RealSphericalHarmonic* rsh;
    std::vector<RealSphericalHarmonic*> functionVector;
    ifstrstream shapeStream;
    std::string tempString;
    std::string ffPath;
    char* tempPath; 
    
    tempPath = getenv("FORCE_PARAM_PATH");
    
    if (tempPath == NULL) {
      //convert a macro from compiler to a string in c++
      STR_DEFINE(ffPath, FRC_PATH );
    } else {
      ffPath = tempPath;
    }
    
    shapeStream.open( shapeFileName.c_str() );
    
    if( !shapeStream.is_open() ){
      
      tempString = ffPath;
      tempString += "/";
      tempString += shapeFileName;
      shapeFileName = tempString;
      
      shapeStream.open( shapeFileName.c_str() );
      
      if( !shapeStream.is_open() ){
        
        sprintf( painCave.errMsg,
                 "Error opening the shape file:\n"
                 "\t%s\n"
                 "\tHave you tried setting the FORCE_PARAM_PATH environment "
                 "variable?\n",
                 shapeFileName.c_str() );
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();
      }
    }
    
    

    ShapeAtomType* st = new ShapeAtomType();
    
    // first parse the info. in the ShapeInfo section
    findBegin(shapeStream, "ShapeInfo");
    shapeStream.getline(buffer, bufferSize);
    
    // loop over the interior of the ShapeInfo section
    while( !shapeStream.eof() ) {
      // toss comment lines
      if( buffer[0] != '!' && buffer[0] != '#' ){
        // end marks section completion
        if (isEndLine(buffer)) break;
        StringTokenizer tokenInfo(buffer);
        // blank lines are ignored
        if (tokenInfo.countTokens() != 0) {
          if (tokenInfo.countTokens() < 5) {
            sprintf(painCave.errMsg, 
                    "ShapesAtomTypesSectionParser Error: Not enough "
                    "information on a ShapeInfo line in file: %s\n", 
                    shapeFileName.c_str() );
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal = 1;
            simError();  
          } else {
            tokenInfo.skipToken();
            at->setMass( tokenInfo.nextTokenAsDouble() );
            DirectionalAdapter da = DirectionalAdapter(at);
            Mat3x3d I;
            I(0,0) = tokenInfo.nextTokenAsDouble();
            I(1,1) = tokenInfo.nextTokenAsDouble();
            I(2,2) = tokenInfo.nextTokenAsDouble();
            da.makeDirectional(I);
          }
        }
      }
      shapeStream.getline(buffer, bufferSize);
    }
    
    // now grab the contact functions
    findBegin(shapeStream, "ContactFunctions");
    functionVector.clear();
    
    shapeStream.getline(buffer, bufferSize);
    while( !shapeStream.eof() ) {
      // toss comment lines
      if( buffer[0] != '!' && buffer[0] != '#' ){
        // end marks section completion
        if (isEndLine(buffer)) break;
        StringTokenizer tokenInfo1(buffer);
        // blank lines are ignored
        if (tokenInfo1.countTokens() != 0) {
          if (tokenInfo1.countTokens() < 4) {
            sprintf( painCave.errMsg,
                     "ShapesAtomTypesSectionParser Error: Not enough "
                     "information on a ContactFunctions line in file: %s\n", 
                     shapeFileName.c_str() );
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal = 1;
            simError();
          } else {
            // read in a spherical harmonic function
            rsh = new RealSphericalHarmonic();
            rsh->setL( tokenInfo1.nextTokenAsInt() );
            rsh->setM( tokenInfo1.nextTokenAsInt() );
            token = tokenInfo1.nextToken();
            toLower(token);
            if (token == "sin")
              rsh->makeSinFunction();
            else
              rsh->makeCosFunction();
            rsh->setCoefficient( tokenInfo1.nextTokenAsDouble() );
            
            functionVector.push_back(rsh);
          }
        }
      }
      shapeStream.getline(buffer, bufferSize);
    }
    
    // pass contact functions to ShapeType
    st->setContactFuncs(functionVector);
    
    // now grab the range functions
    findBegin(shapeStream, "RangeFunctions");
    functionVector.clear();
    
    shapeStream.getline(buffer, bufferSize);
    while( !shapeStream.eof() ) {
      // toss comment lines
      if( buffer[0] != '!' && buffer[0] != '#' ){
        // end marks section completion
        if (isEndLine(buffer)) break;
        StringTokenizer tokenInfo2(buffer);
        // blank lines are ignored
        if (tokenInfo2.countTokens() != 0) {
          if (tokenInfo2.countTokens() < 4) {
            sprintf( painCave.errMsg,
                     "ShapesAtomTypesSectionParser Error: Not enough "
                     "information on a RangeFunctions line in file: %s\n", 
                     shapeFileName.c_str() );
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal = 1;
            simError();
          } else {
            // read in a spherical harmonic function
            rsh = new RealSphericalHarmonic();
            rsh->setL( tokenInfo2.nextTokenAsInt() );
            rsh->setM( tokenInfo2.nextTokenAsInt() );
            token = tokenInfo2.nextToken();
            toLower(token);
            if (token == "sin")
              rsh->makeSinFunction();
            else
              rsh->makeCosFunction();
            rsh->setCoefficient( tokenInfo2.nextTokenAsDouble() );
            
            functionVector.push_back(rsh);
          }
        }
      }
      shapeStream.getline(buffer, bufferSize);
    }
    
    // pass range functions to ShapeType
    st->setRangeFuncs(functionVector);
    
    // finally grab the strength functions
    findBegin(shapeStream, "StrengthFunctions");
    functionVector.clear();
    
    shapeStream.getline(buffer, bufferSize);
    while( !shapeStream.eof() ) {
      // toss comment lines
      if( buffer[0] != '!' && buffer[0] != '#' ){
        // end marks section completion
        if (isEndLine(buffer)) break;
        StringTokenizer tokenInfo3(buffer);
        // blank lines are ignored
        if (tokenInfo3.countTokens() != 0) {
          if (tokenInfo3.countTokens() < 4) {
            sprintf( painCave.errMsg,
                     "ShapesAtomTypesSectionParser Error: Not enough "
                     "information on a StrengthFunctions line in file: %s\n", 
                     shapeFileName.c_str() );
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal = 1;
            simError();
          } else {
            // read in a spherical harmonic function
            rsh = new RealSphericalHarmonic();
            rsh->setL( tokenInfo3.nextTokenAsInt() );
            rsh->setM( tokenInfo3.nextTokenAsInt() );
            token = tokenInfo3.nextToken();
            toLower(token);
            if (token == "sin")
              rsh->makeSinFunction();
            else
              rsh->makeCosFunction();
            rsh->setCoefficient( tokenInfo3.nextTokenAsDouble() );
            
            functionVector.push_back(rsh);
          }
        }
      }
      shapeStream.getline(buffer, bufferSize);
    }
    
    // pass strength functions to ShapeType
    st->setStrengthFuncs(functionVector);
    at->addProperty(std::shared_ptr<GenericData>(new ShapeAtypeData("Shape", st)));     
  //  delete shapeStream;
  }
} //end namespace OpenMD

