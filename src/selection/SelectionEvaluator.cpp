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
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */

#include <stack>
#include "selection/SelectionEvaluator.hpp"
#include "primitives/Atom.hpp"
#include "primitives/DirectionalAtom.hpp"
#include "primitives/RigidBody.hpp"
#include "primitives/Molecule.hpp"
#include "io/ifstrstream.hpp"

namespace OpenMD {


  SelectionEvaluator::SelectionEvaluator(SimInfo* si) 
    : info(si), nameFinder(info), distanceFinder(info), hullFinder(info),
      indexFinder(info), 
      isLoaded_(false){    
      nStuntDouble = info->getNGlobalAtoms() + info->getNGlobalRigidBodies();
    }            

  bool SelectionEvaluator::loadScript(const std::string& filename, 
                                      const std::string& script) {
    clearDefinitionsAndLoadPredefined();
    this->filename = filename;
    this->script = script;
    if (! compiler.compile(filename, script)) {
      error = true;
      errorMessage = compiler.getErrorMessage();

      sprintf( painCave.errMsg,
               "SelectionCompiler Error: %s\n", errorMessage.c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
      return false;
    }

    pc = 0;
    aatoken = compiler.getAatokenCompiled();
    linenumbers = compiler.getLineNumbers();
    lineIndices = compiler.getLineIndices();

    std::vector<std::vector<Token> >::const_iterator i;   

    isDynamic_ = false;
    for (i = aatoken.begin(); i != aatoken.end(); ++i) {
      if (containDynamicToken(*i)) {
	isDynamic_ = true;
	break;
      }
    }

    isLoaded_ = true;
    return true;
  }

  void SelectionEvaluator::clearState() {
    error = false;
    errorMessage = "";
  }

  bool SelectionEvaluator::loadScriptString(const std::string& script) {
    clearState();
    return loadScript("", script);
  }

  bool SelectionEvaluator::loadScriptFile(const std::string& filename) {
    clearState();
    return loadScriptFileInternal(filename);
  }

  bool SelectionEvaluator::loadScriptFileInternal(const std::string & filename) {
    ifstrstream ifs(filename.c_str());
    if (!ifs.is_open()) {
      return false;
    }
    
    const int bufferSize = 65535;
    char buffer[bufferSize];
    std::string script;
    while(ifs.getline(buffer, bufferSize)) {
      script += buffer;
    }
    return loadScript(filename, script);
  }
  
  void SelectionEvaluator::instructionDispatchLoop(OpenMDBitSet& bs){
    
    while ( pc < aatoken.size()) {
      statement = aatoken[pc++];
      statementLength = statement.size();
      Token token = statement[0];
      switch (token.tok) {
      case Token::define:
	define();
	break;
      case Token::select:
	select(bs);
	break;
      default:
	unrecognizedCommand(token);
	return;
      }
    }

  }

  OpenMDBitSet SelectionEvaluator::expression(const std::vector<Token>& code, 
                                             int pcStart) {
    OpenMDBitSet bs;
    std::stack<OpenMDBitSet> stack; 
   
    for (int pc = pcStart; pc < code.size(); ++pc) {
      Token instruction = code[pc];

      switch (instruction.tok) {
      case Token::expressionBegin:
        break;
      case Token::expressionEnd:
        break;
      case Token::all:
        bs = OpenMDBitSet(nStuntDouble);
        bs.setAll();
        stack.push(bs);            
        break;
      case Token::none:
        bs = OpenMDBitSet(nStuntDouble);
        stack.push(bs);            
        break;
      case Token::opOr:
        bs = stack.top();
        stack.pop();
        stack.top() |= bs;
        break;
      case Token::opAnd:
        bs = stack.top();
        stack.pop();
        stack.top() &= bs;
        break;
      case Token::opNot:
        stack.top().flip();
        break;
      case Token::within:
        withinInstruction(instruction, stack.top());
        break;
      case Token::hull:
        stack.push(hull());
        break;
	//case Token::selected:
	//  stack.push(getSelectionSet());
	//  break;
      case Token::name:
        stack.push(nameInstruction(boost::any_cast<std::string>(instruction.value)));
        break;
      case Token::index:
        stack.push(indexInstruction(instruction.value));
        break;
      case Token::identifier:
        stack.push(lookupValue(boost::any_cast<std::string>(instruction.value)));
        break;
      case Token::opLT:
      case Token::opLE:
      case Token::opGE:
      case Token::opGT:
      case Token::opEQ:
      case Token::opNE:
        stack.push(comparatorInstruction(instruction));
        break;
      default:
        unrecognizedExpression();
      }
    }
    if (stack.size() != 1)
      evalError("atom expression compiler error - stack over/underflow");
          
    return stack.top();
  }



  OpenMDBitSet SelectionEvaluator::comparatorInstruction(const Token& instruction) {
    int comparator = instruction.tok;
    int property = instruction.intValue;
    float comparisonValue = boost::any_cast<float>(instruction.value);
    float propertyValue;
    OpenMDBitSet bs(nStuntDouble);
    bs.clearAll();
    
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::AtomIterator ai;
    Atom* atom;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;

    for (mol = info->beginMolecule(mi); mol != NULL; 
         mol = info->nextMolecule(mi)) {

      for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
	compareProperty(atom, bs, property, comparator, comparisonValue);
      }
     
      for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
           rb = mol->nextRigidBody(rbIter)) {
     	compareProperty(rb, bs, property, comparator, comparisonValue);
      }
    }

    return bs;
  }

  void SelectionEvaluator::compareProperty(StuntDouble* sd, OpenMDBitSet& bs, 
                                           int property, int comparator, 
                                           float comparisonValue) {
    RealType propertyValue = 0.0;
    switch (property) {
    case Token::mass:
      propertyValue = sd->getMass();
      break;
    case Token::charge:
      if (sd->isAtom()){
	Atom* atom = static_cast<Atom*>(sd);
	propertyValue = getCharge(atom);
      } else if (sd->isRigidBody()) {
	RigidBody* rb = static_cast<RigidBody*>(sd);
	RigidBody::AtomIterator ai;
	Atom* atom;
	for (atom = rb->beginAtom(ai); atom != NULL; atom = rb->nextAtom(ai)) {
	  propertyValue+=  getCharge(atom);
	}
      }
      break;
    case Token::x:
      propertyValue = sd->getPos().x();
      break;
    case Token::y:
      propertyValue = sd->getPos().y();
      break;
    case Token::z:
      propertyValue = sd->getPos().z();
      break;
    case Token::r:
      propertyValue = sd->getPos().length();
      break;
    default:
      unrecognizedAtomProperty(property);
    }
        
    bool match = false;
    switch (comparator) {
    case Token::opLT:
      match = propertyValue < comparisonValue;
      break;
    case Token::opLE:
      match = propertyValue <= comparisonValue;
      break;
    case Token::opGE:
      match = propertyValue >= comparisonValue;
      break;
    case Token::opGT:
      match = propertyValue > comparisonValue;
      break;
    case Token::opEQ:
      match = propertyValue == comparisonValue;
      break;
    case Token::opNE:
      match = propertyValue != comparisonValue;
      break;
    }
    if (match)
      bs.setBitOn(sd->getGlobalIndex());
    
  }

  void SelectionEvaluator::withinInstruction(const Token& instruction, 
                                             OpenMDBitSet& bs){
    
    boost::any withinSpec = instruction.value;
    float distance;
    if (withinSpec.type() == typeid(float)){
      distance = boost::any_cast<float>(withinSpec);
    } else if (withinSpec.type() == typeid(int)) {
      distance = boost::any_cast<int>(withinSpec);    
    } else {
      evalError("casting error in withinInstruction");
      bs.clearAll();
    }
    
    bs = distanceFinder.find(bs, distance);            
  }
  
  void SelectionEvaluator::define() {
    assert(statement.size() >= 3);
    
    std::string variable = boost::any_cast<std::string>(statement[1].value);
    
    variables.insert(VariablesType::value_type(variable, 
                                               expression(statement, 2)));
  }
  

  /** @todo */
  void SelectionEvaluator::predefine(const std::string& script) {
    
    if (compiler.compile("#predefine", script)) {
      std::vector<std::vector<Token> > aatoken = compiler.getAatokenCompiled();
      if (aatoken.size() != 1) {
	evalError("predefinition does not have exactly 1 command:"
		  + script);
	return;
      }
      std::vector<Token> statement = aatoken[0];
      if (statement.size() > 2) {
	int tok = statement[1].tok;
	if (tok == Token::identifier || 
            (tok & Token::predefinedset) == Token::predefinedset) {
	  std::string variable = boost::any_cast<std::string>(statement[1].value);
	  variables.insert(VariablesType::value_type(variable, statement));
          
	} else {
	  evalError("invalid variable name:" + script);
	}
      }else {
	evalError("bad predefinition length:" + script);
      }      
        
    } else {
      evalError("predefined set compile error:" + script +
		"\ncompile error:" + compiler.getErrorMessage());
    }
  }

  void SelectionEvaluator::select(OpenMDBitSet& bs){
    bs = expression(statement, 1);
  }
  
  OpenMDBitSet SelectionEvaluator::lookupValue(const std::string& variable){
    
    OpenMDBitSet bs(nStuntDouble);
    std::map<std::string, boost::any>::iterator i = variables.find(variable);
    
    if (i != variables.end()) {
      if (i->second.type() == typeid(OpenMDBitSet)) {
	return boost::any_cast<OpenMDBitSet>(i->second);
      } else if (i->second.type() ==  typeid(std::vector<Token>)){
	bs = expression(boost::any_cast<std::vector<Token> >(i->second), 2);
	i->second =  bs; /**@todo fixme */
	return bs;
      }
    } else {
      unrecognizedIdentifier(variable);
    }
    
    return bs;
  }
  
  OpenMDBitSet SelectionEvaluator::nameInstruction(const std::string& name){    
    return nameFinder.match(name);    
  }    

  bool SelectionEvaluator::containDynamicToken(const std::vector<Token>& tokens){
    std::vector<Token>::const_iterator i;
    for (i = tokens.begin(); i != tokens.end(); ++i) {
      if (i->tok & Token::dynamic) {
	return true;
      }
    }
    
    return false;
  }    

  void SelectionEvaluator::clearDefinitionsAndLoadPredefined() {
    variables.clear();
    //load predefine
    //predefine();
  }

  OpenMDBitSet SelectionEvaluator::evaluate() {
    OpenMDBitSet bs(nStuntDouble);
    if (isLoaded_) {
      pc = 0;
      instructionDispatchLoop(bs);
    }

    return bs;
  }

  OpenMDBitSet SelectionEvaluator::indexInstruction(const boost::any& value) {
    OpenMDBitSet bs(nStuntDouble);

    if (value.type() == typeid(int)) {
      int index = boost::any_cast<int>(value);
      if (index < 0 || index >= bs.size()) {
	invalidIndex(index);
      } else {
	bs = indexFinder.find(index);
      }
    } else if (value.type() == typeid(std::pair<int, int>)) {
      std::pair<int, int> indexRange= boost::any_cast<std::pair<int, int> >(value);
      assert(indexRange.first <= indexRange.second);
      if (indexRange.first < 0 || indexRange.second >= bs.size()) {
	invalidIndexRange(indexRange);
      }else {
	bs = indexFinder.find(indexRange.first, indexRange.second);
      }
    }

    return bs;
  }


  OpenMDBitSet SelectionEvaluator::hull() {
    OpenMDBitSet bs(nStuntDouble);
    
    bs = hullFinder.findHull();
    
    return bs;
  }



  RealType SelectionEvaluator::getCharge(Atom* atom) {
    RealType charge =0.0;
    AtomType* atomType = atom->getAtomType();
    if (atomType->isCharge()) {
      GenericData* data = atomType->getPropertyByName("Charge");
      if (data != NULL) {
	DoubleGenericData* doubleData= dynamic_cast<DoubleGenericData*>(data);

	if (doubleData != NULL) {
	  charge = doubleData->getData();

	} else {
	  sprintf( painCave.errMsg,
		   "Can not cast GenericData to DoubleGenericData\n");
	  painCave.severity = OPENMD_ERROR;
	  painCave.isFatal = 1;
	  simError();          
	}
      } 
    }

    return charge;
  }

}
