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

#include <stack>
#include "selection/SelectionEvaluator.hpp"
#include "primitives/Atom.hpp"
#include "primitives/DirectionalAtom.hpp"
#include "primitives/RigidBody.hpp"
#include "primitives/Molecule.hpp"
#include "io/ifstrstream.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"


namespace OpenMD {

  SelectionEvaluator::SelectionEvaluator(SimInfo* si) 
    : info(si), nameFinder(info), distanceFinder(info), hullFinder(info),
      alphaHullFinder(info), indexFinder(info), isLoaded_(false),
      hasSurfaceArea_(false), hasVolume_(false) {
    nObjects.push_back(info->getNGlobalAtoms() + info->getNGlobalRigidBodies());
    nObjects.push_back(info->getNGlobalBonds());
    nObjects.push_back(info->getNGlobalBends());
    nObjects.push_back(info->getNGlobalTorsions());
    nObjects.push_back(info->getNGlobalInversions());
    nObjects.push_back(info->getNGlobalMolecules());    
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
  
  void SelectionEvaluator::instructionDispatchLoop(SelectionSet& bs){
    
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

  void SelectionEvaluator::instructionDispatchLoop(SelectionSet& bs, int frame){
    
    while ( pc < aatoken.size()) {
      statement = aatoken[pc++];
      statementLength = statement.size();
      Token token = statement[0];
      switch (token.tok) {
      case Token::define:
	define();
	break;
      case Token::select:
	select(bs, frame);
	break;
      default:
	unrecognizedCommand(token);
	return;
      }
    }

  }

  SelectionSet SelectionEvaluator::expression(const std::vector<Token>& code, 
                                             int pcStart) {
    SelectionSet bs = createSelectionSets();
    std::stack<SelectionSet> stack; 
    vector<int> bsSize = bs.size();
   
    for (unsigned int pc = pcStart; pc < code.size(); ++pc) {
      Token instruction = code[pc];

      switch (instruction.tok) {
      case Token::expressionBegin:
        break;
      case Token::expressionEnd:
        break;
      case Token::all:
        bs = allInstruction();
        stack.push(bs);
        break;
      case Token::none:
        bs = createSelectionSets();
        stack.push(bs);
        break;
      case Token::opOr:
          bs= stack.top();
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
      case Token::alphahull:
        stack.push(alphaHullInstruction(instruction));
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


  SelectionSet SelectionEvaluator::expression(const std::vector<Token>& code, 
                                              int pcStart, int frame) {
    SelectionSet bs = createSelectionSets();
    std::stack<SelectionSet> stack; 
   
    for (unsigned int pc = pcStart; pc < code.size(); ++pc) {
      Token instruction = code[pc];

      switch (instruction.tok) {
      case Token::expressionBegin:
        break;
      case Token::expressionEnd:
        break;
      case Token::all:
        bs = allInstruction();
        stack.push(bs);            
        break;
      case Token::none:
        bs = SelectionSet(nObjects);
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
        withinInstruction(instruction, stack.top(), frame);
        break;
      case Token::alphahull:
        stack.push(alphaHullInstruction(instruction, frame));
        break;        
      case Token::hull:
        stack.push(hull(frame));
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
        stack.push(comparatorInstruction(instruction, frame));
        break;
      default:
        unrecognizedExpression();
      }
    }
    if (stack.size() != 1)
      evalError("atom expression compiler error - stack over/underflow");
          
    return stack.top();
  }



  SelectionSet SelectionEvaluator::comparatorInstruction(const Token& instruction) {
    int comparator = instruction.tok;
    int property = instruction.intValue;
    float comparisonValue = boost::any_cast<float>(instruction.value);
    SelectionSet bs = createSelectionSets();
    bs.clearAll();
    
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::AtomIterator ai;
    Atom* atom;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;

    for (mol = info->beginMolecule(mi); mol != NULL; 
         mol = info->nextMolecule(mi)) {

      compareProperty(mol, bs, property, comparator, comparisonValue);

      for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
	compareProperty(atom, bs, property, comparator, comparisonValue);
      }
     
      for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
           rb = mol->nextRigidBody(rbIter)) {
     	compareProperty(rb, bs, property, comparator, comparisonValue);
      }
    }

    return bs.parallelReduce();
  }

  SelectionSet SelectionEvaluator::comparatorInstruction(const Token& instruction, int frame) {
    int comparator = instruction.tok;
    int property = instruction.intValue;
    float comparisonValue = boost::any_cast<float>(instruction.value);
    SelectionSet bs = createSelectionSets();
    bs.clearAll();
    
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::AtomIterator ai;
    Atom* atom;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;

    for (mol = info->beginMolecule(mi); mol != NULL; 
         mol = info->nextMolecule(mi)) {

      compareProperty(mol, bs, property, comparator, comparisonValue, frame);

      for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
	compareProperty(atom, bs, property, comparator, comparisonValue, frame);
      }
     
      for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
           rb = mol->nextRigidBody(rbIter)) {
     	compareProperty(rb, bs, property, comparator, comparisonValue, frame);
      }
    }

    return bs.parallelReduce();
  }

  void SelectionEvaluator::compareProperty(StuntDouble* sd, SelectionSet& bs, 
                                           int property, int comparator, 
                                           float comparisonValue) {
    RealType propertyValue = 0.0;
    Vector3d pos;

    switch (property) {
    case Token::atomno:
      if (sd->isAtom()) {
	Atom* atom = static_cast<Atom*>(sd);
	propertyValue = atom->getGlobalIndex();
      }
      break;
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
    case Token::wrappedX:    
      pos = sd->getPos();
      info->getSnapshotManager()->getCurrentSnapshot()->wrapVector(pos);
      propertyValue = pos.x();
      break;
    case Token::wrappedY:
      pos = sd->getPos();
      info->getSnapshotManager()->getCurrentSnapshot()->wrapVector(pos);
      propertyValue = pos.y();
      break;
    case Token::wrappedZ:
      pos = sd->getPos();
      info->getSnapshotManager()->getCurrentSnapshot()->wrapVector(pos);
      propertyValue = pos.z();
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
      match = fabs(propertyValue - comparisonValue) <= OpenMD::epsilon;
      break;
    case Token::opNE:
      match = propertyValue != comparisonValue;
      break;
    }

    if (match) 
      bs.bitsets_[STUNTDOUBLE].setBitOn(sd->getGlobalIndex());

    
  }

  void SelectionEvaluator::compareProperty(Molecule* mol, SelectionSet& bs, 
                                           int property, int comparator, 
                                           float comparisonValue) {
    RealType propertyValue = 0.0;
    Vector3d pos;

    switch (property) {
    case Token::mass:
      propertyValue = mol->getMass();
      break;
    case Token::charge:
      {
        Molecule::AtomIterator ai;
        Atom* atom;
        for (atom = mol->beginAtom(ai); atom != NULL;
             atom = mol->nextAtom(ai)) {
          propertyValue+=  getCharge(atom);
        }
      }
      break;
    case Token::x:
      propertyValue = mol->getCom().x();
      break;
    case Token::y:
      propertyValue = mol->getCom().y();
      break;
    case Token::z:
      propertyValue = mol->getCom().z();
      break;
    case Token::wrappedX:    
      pos = mol->getCom();
      info->getSnapshotManager()->getCurrentSnapshot()->wrapVector(pos);
      propertyValue = pos.x();
      break;
    case Token::wrappedY:
      pos = mol->getCom();
      info->getSnapshotManager()->getCurrentSnapshot()->wrapVector(pos);
      propertyValue = pos.y();
      break;
    case Token::wrappedZ:
      pos = mol->getCom();
      info->getSnapshotManager()->getCurrentSnapshot()->wrapVector(pos);
      propertyValue = pos.z();
      break;
    case Token::r:
      propertyValue = mol->getCom().length();
      break;
    default:
      unrecognizedMoleculeProperty(property);
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
      match = fabs(propertyValue - comparisonValue) <= OpenMD::epsilon;
      break;
    case Token::opNE:
      match = propertyValue != comparisonValue;
      break;
    }
    
    if (match) 
      bs.bitsets_[MOLECULE].setBitOn(mol->getGlobalIndex());    
  }

  void SelectionEvaluator::compareProperty(Molecule* mol, SelectionSet& bs, 
                                           int property, int comparator, 
                                           float comparisonValue, int frame) {
    RealType propertyValue = 0.0;
    Vector3d pos;
    switch (property) {
    case Token::mass:
      propertyValue = mol->getMass();
      break;
    case Token::charge:
      {
        Molecule::AtomIterator ai;
        Atom* atom;
        for (atom = mol->beginAtom(ai); atom != NULL;
             atom = mol->nextAtom(ai)) {
          propertyValue+=  getCharge(atom,frame);
        }
      }      
      break;
    case Token::x:
      propertyValue = mol->getCom(frame).x();
      break;
    case Token::y:
      propertyValue = mol->getCom(frame).y();
      break;
    case Token::z:
      propertyValue = mol->getCom(frame).z();
      break;
    case Token::wrappedX:    
      pos = mol->getCom(frame);
      info->getSnapshotManager()->getSnapshot(frame)->wrapVector(pos);
      propertyValue = pos.x();
      break;
    case Token::wrappedY:
      pos = mol->getCom(frame);
      info->getSnapshotManager()->getSnapshot(frame)->wrapVector(pos);
      propertyValue = pos.y();
      break;
    case Token::wrappedZ:
      pos = mol->getCom(frame);
      info->getSnapshotManager()->getSnapshot(frame)->wrapVector(pos);
      propertyValue = pos.z();
      break;

    case Token::r:
      propertyValue = mol->getCom(frame).length();
      break;
    default:
      unrecognizedMoleculeProperty(property);
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
      match = fabs(propertyValue - comparisonValue) <= OpenMD::epsilon;
      break;
    case Token::opNE:
      match = propertyValue != comparisonValue;
      break;
    }
    if (match) 
      bs.bitsets_[MOLECULE].setBitOn(mol->getGlobalIndex());
    
    
  }
  void SelectionEvaluator::compareProperty(StuntDouble* sd, SelectionSet& bs, 
                                           int property, int comparator, 
                                           float comparisonValue, int frame) {
    RealType propertyValue = 0.0;
    Vector3d pos;
    switch (property) {
    case Token::atomno:
      if (sd->isAtom()) {
	Atom* atom = static_cast<Atom*>(sd);
	propertyValue = atom->getGlobalIndex();
      }
      break;      
    case Token::mass:
      propertyValue = sd->getMass();
      break;
    case Token::charge:
      if (sd->isAtom()){
	Atom* atom = static_cast<Atom*>(sd);
	propertyValue = getCharge(atom,frame);
      } else if (sd->isRigidBody()) {
	RigidBody* rb = static_cast<RigidBody*>(sd);
	RigidBody::AtomIterator ai;
	Atom* atom;
	for (atom = rb->beginAtom(ai); atom != NULL; atom = rb->nextAtom(ai)) {
	  propertyValue+=  getCharge(atom,frame);
	}
      }
      break;
    case Token::x:
      propertyValue = sd->getPos(frame).x();
      break;
    case Token::y:
      propertyValue = sd->getPos(frame).y();
      break;
    case Token::z:
      propertyValue = sd->getPos(frame).z();
      break;
    case Token::wrappedX:    
      pos = sd->getPos(frame);
      info->getSnapshotManager()->getSnapshot(frame)->wrapVector(pos);
      propertyValue = pos.x();
      break;
    case Token::wrappedY:
      pos = sd->getPos(frame);
      info->getSnapshotManager()->getSnapshot(frame)->wrapVector(pos);
      propertyValue = pos.y();
      break;
    case Token::wrappedZ:
      pos = sd->getPos(frame);
      info->getSnapshotManager()->getSnapshot(frame)->wrapVector(pos);
      propertyValue = pos.z();
      break;

    case Token::r:
      propertyValue = sd->getPos(frame).length();
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
      match = fabs(propertyValue - comparisonValue) <= OpenMD::epsilon;
      break;
    case Token::opNE:
      match = propertyValue != comparisonValue;
      break;
    }
    if (match) 
      bs.bitsets_[STUNTDOUBLE].setBitOn(sd->getGlobalIndex());    
  }
  
  
  void SelectionEvaluator::withinInstruction(const Token& instruction, 
                                             SelectionSet& bs){
    
    boost::any withinSpec = instruction.value;
    float distance(0.0);
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

  void SelectionEvaluator::withinInstruction(const Token& instruction, 
                                             SelectionSet& bs, int frame){
    
    boost::any withinSpec = instruction.value;
    float distance(0.0);
    if (withinSpec.type() == typeid(float)){
      distance = boost::any_cast<float>(withinSpec);
    } else if (withinSpec.type() == typeid(int)) {
      distance = boost::any_cast<int>(withinSpec);    
    } else {
      evalError("casting error in withinInstruction");
      bs.clearAll();
    }
    
    bs = distanceFinder.find(bs, distance, frame);            
  }

  SelectionSet SelectionEvaluator::alphaHullInstruction(const Token& instruction){


    SelectionSet bs = createSelectionSets();

    boost::any alphaSpec = instruction.value;
    float alpha(0.0);
    if (alphaSpec.type() == typeid(float)){
      alpha = boost::any_cast<float>(alphaSpec);
    } else if (alphaSpec.type() == typeid(int)) {
      alpha = boost::any_cast<int>(alphaSpec);    
    } else {
      evalError("casting error in alphaHullInstruction");
      bs.clearAll();
    }

    alphaHullFinder.setAlpha(alpha);
    bs = alphaHullFinder.findHull();
    surfaceArea_ = alphaHullFinder.getSurfaceArea();
    hasSurfaceArea_ = true;
    volume_ = alphaHullFinder.getVolume();
    hasVolume_ = true;

    return bs.parallelReduce();
  }

  SelectionSet SelectionEvaluator::alphaHullInstruction(const Token& instruction, int frame){
    SelectionSet bs = createSelectionSets();

    boost::any alphaSpec = instruction.value;
    float alpha(0.0);
    if (alphaSpec.type() == typeid(float)){
      alpha = boost::any_cast<float>(alphaSpec);
    } else if (alphaSpec.type() == typeid(int)) {
      alpha = boost::any_cast<int>(alphaSpec);    
    } else {
      evalError("casting error in alphaHullInstruction");
      bs.clearAll();
    }
    
    alphaHullFinder.setAlpha(alpha);
    bs = alphaHullFinder.findHull(frame);
    surfaceArea_ = alphaHullFinder.getSurfaceArea();
    hasSurfaceArea_ = true;
    volume_ = alphaHullFinder.getVolume();
    hasVolume_ = true;

    return bs.parallelReduce();
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

  void SelectionEvaluator::select(SelectionSet& bs){
    bs = expression(statement, 1);
  }

  void SelectionEvaluator::select(SelectionSet& bs, int frame){
    bs = expression(statement, 1, frame);
  }
  
  SelectionSet SelectionEvaluator::lookupValue(const std::string& variable){
    
    SelectionSet bs = createSelectionSets();
    std::map<std::string, boost::any>::iterator i = variables.find(variable);
    
    if (i != variables.end()) {
      if (i->second.type() == typeid(SelectionSet)) {
	return boost::any_cast<SelectionSet>(i->second);
      } else if (i->second.type() ==  typeid(std::vector<Token>)){
	bs = expression(boost::any_cast<std::vector<Token> >(i->second), 2);
	i->second =  bs; /**@todo fixme */
	return bs.parallelReduce();
      }
    } else {
      unrecognizedIdentifier(variable);
    }
    
    return bs.parallelReduce();
  }
  
  SelectionSet SelectionEvaluator::nameInstruction(const std::string& name){    
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

  SelectionSet SelectionEvaluator::createSelectionSets() {
    SelectionSet ss(nObjects);
    return ss;
  }

  SelectionSet SelectionEvaluator::evaluate() {
    SelectionSet bs = createSelectionSets();
    if (isLoaded_) {
      pc = 0;
      instructionDispatchLoop(bs);
    }
    return bs.parallelReduce();
  }

  SelectionSet SelectionEvaluator::evaluate(int frame) {
    SelectionSet bs = createSelectionSets();
    if (isLoaded_) {
      pc = 0;
      instructionDispatchLoop(bs, frame);
    }
    return bs.parallelReduce();
  }

  SelectionSet SelectionEvaluator::indexInstruction(const boost::any& value) {
    SelectionSet bs = createSelectionSets();

    if (value.type() == typeid(int)) {
      int index = boost::any_cast<int>(value);
      if (index < 0 || index >= bs.bitsets_[STUNTDOUBLE].size()) {
	invalidIndex(index);
      } else {
	bs = indexFinder.find(index);
      }
    } else if (value.type() == typeid(std::pair<int, int>)) {
      std::pair<int, int> indexRange= boost::any_cast<std::pair<int, int> >(value);
      assert(indexRange.first <= indexRange.second);
      if (indexRange.first < 0 || 
          indexRange.second >= bs.bitsets_[STUNTDOUBLE].size()) {
	invalidIndexRange(indexRange);
      }else {
	bs = indexFinder.find(indexRange.first, indexRange.second);
      }
    }

    return bs.parallelReduce();
  }

  SelectionSet SelectionEvaluator::allInstruction() {
    SelectionSet ss = createSelectionSets();

    SimInfo::MoleculeIterator mi;
    Molecule::AtomIterator ai;
    Molecule::RigidBodyIterator rbIter;
    Molecule::BondIterator bondIter;
    Molecule::BendIterator bendIter;
    Molecule::TorsionIterator torsionIter;
    Molecule::InversionIterator inversionIter;

    Molecule* mol;
    Atom* atom;
    RigidBody* rb;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Inversion* inversion;    

    // Doing the loop insures that we're actually on this processor.

    for (mol = info->beginMolecule(mi); mol != NULL; 
         mol = info->nextMolecule(mi)) {
      for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
        ss.bitsets_[STUNTDOUBLE].setBitOn(atom->getGlobalIndex());
      }     
      for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
           rb = mol->nextRigidBody(rbIter)) {
        ss.bitsets_[STUNTDOUBLE].setBitOn(rb->getGlobalIndex());
      }
      for (bond = mol->beginBond(bondIter); bond != NULL; 
           bond = mol->nextBond(bondIter)) {
        ss.bitsets_[BOND].setBitOn(bond->getGlobalIndex());
      }   
      for (bend = mol->beginBend(bendIter); bend != NULL; 
           bend = mol->nextBend(bendIter)) {
        ss.bitsets_[BEND].setBitOn(bend->getGlobalIndex());
      }   
      for (torsion = mol->beginTorsion(torsionIter); torsion != NULL; 
           torsion = mol->nextTorsion(torsionIter)) {
        ss.bitsets_[TORSION].setBitOn(torsion->getGlobalIndex());
      }   
      for (inversion = mol->beginInversion(inversionIter); inversion != NULL; 
           inversion = mol->nextInversion(inversionIter)) {
        ss.bitsets_[INVERSION].setBitOn(inversion->getGlobalIndex());
      }
      ss.bitsets_[MOLECULE].setBitOn(mol->getGlobalIndex());
    }

    return ss;
  }

  SelectionSet SelectionEvaluator::hull() {
    SelectionSet bs = createSelectionSets();
    
    bs = hullFinder.findHull();
    surfaceArea_ = hullFinder.getSurfaceArea();
    hasSurfaceArea_ = true;
    volume_ = hullFinder.getVolume();
    hasVolume_ = true;

    return bs.parallelReduce();
  }


  SelectionSet SelectionEvaluator::hull(int frame) {
    SelectionSet bs = createSelectionSets();
    
    bs = hullFinder.findHull(frame);

    return bs.parallelReduce();
  }
  
  RealType SelectionEvaluator::getCharge(Atom* atom) {
    RealType charge = 0.0;
    AtomType* atomType = atom->getAtomType();

    FixedChargeAdapter fca = FixedChargeAdapter(atomType);
    if ( fca.isFixedCharge() ) {
      charge = fca.getCharge();
    }
    
    FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
    if ( fqa.isFluctuatingCharge() ) {
      charge += atom->getFlucQPos();
    }
    return charge;
  }

  RealType SelectionEvaluator::getCharge(Atom* atom, int frame) {
    RealType charge = 0.0;
    AtomType* atomType = atom->getAtomType();    

    FixedChargeAdapter fca = FixedChargeAdapter(atomType);
    if ( fca.isFixedCharge() ) {
      charge = fca.getCharge();
    }
    
    FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
    if ( fqa.isFluctuatingCharge() ) {
      charge += atom->getFlucQPos(frame);
    }
    return charge;
  }

}
