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

#include "selection/SelectionEvaluator.hpp"
namespace oopse {


bool SelectionEvaluator::loadScript(const std::string& filename, const std::string& script) {
    this->filename = filename;
    this->script = script;
    if (! compiler.compile(filename, script)) {
        error = true;
        errorMessage = compiler.getErrorMessage();
        return false;
    }

    pc = 0;
    aatoken = compiler.getAatokenCompiled();
    linenumbers = compiler.getLineNumbers();
    lineIndices = compiler.getLineIndices();
    return true;
}

void SelectionEvaluator::clearState() {
    for (int i = scriptLevelMax; --i >= 0; )
        stack[i] = null;
    scriptLevel = 0;
    error = false;
    errorMessage = null;
}

bool SelectionEvaluator::loadScriptString(const std::string& script) {
    clearState();
    return loadScript(null, script);
}

bool SelectionEvaluator::loadScriptFile(const std::string& filename) {
    clearState();
    return loadScriptFileInternal(filename);
}


void SelectionEvaluator::instructionDispatchLoop(){

    while ( pc < aatoken.length) {
        statement = aatoken[pc++];
        statementLength = statement.length;
        Token token = statement[0];
        switch (token.tok) {
            case Token.define:
                define();
            break;
            case Token.select:
                select();
            break;
            default:
                unrecognizedCommand(token);
            return;
        }
    }
}

  void SelectionEvaluator::predefine(String script) {
    if (compiler.compile("#predefine", script)) {
      Token [][] aatoken = compiler.getAatokenCompiled();
      if (aatoken.length != 1) {
        viewer.scriptStatus("predefinition does not have exactly 1 command:"
                            + script);
        return;
      }
      Token[] statement = aatoken[0];
      if (statement.length > 2) {
        int tok = statement[1].tok;
        if (tok == Token.identifier ||
            (tok & Token.predefinedset) == Token.predefinedset) {
          String variable = (String)statement[1].value;
          variables.put(variable, statement);
        } else {
          viewer.scriptStatus("invalid variable name:" + script);
        }
      } else {
        viewer.scriptStatus("bad predefinition length:" + script);
      }
    } else {
      viewer.scriptStatus("predefined set compile error:" + script +
                          "\ncompile error:" + compiler.getErrorMessage());
    }
  }



  BitSet SelectionEvaluator::expression(Token[] code, int pcStart) throws ScriptException {
    int numberOfAtoms = viewer.getAtomCount();
    BitSet bs;
    BitSet[] stack = new BitSet[10];
    int sp = 0;

    for (int pc = pcStart; ; ++pc) {
      Token instruction = code[pc];
      if (logMessages)
        viewer.scriptStatus("instruction=" + instruction);
      switch (instruction.tok) {
      case Token.expressionBegin:
        break;
      case Token.expressionEnd:
        break expression_loop;
      case Token.all:
        bs = stack[sp++] = new BitSet(numberOfAtoms);
        for (int i = numberOfAtoms; --i >= 0; )
          bs.set(i);
        break;
      case Token.none:
        stack[sp++] = new BitSet();
        break;
      case Token.opOr:
        bs = stack[--sp];
        stack[sp-1].or(bs);
        break;
      case Token.opAnd:
        bs = stack[--sp];
        stack[sp-1].and(bs);
        break;
      case Token.opNot:
        bs = stack[sp - 1];
        notSet(bs);
        break;
      case Token.within:
        bs = stack[sp - 1];
        stack[sp - 1] = new BitSet();
        withinInstruction(instruction, bs, stack[sp - 1]);
        break;
      case Token.selected:
        stack[sp++] = copyBitSet(viewer.getSelectionSet());
        break;
      case Token.y:
      case Token.amino:
      case Token.backbone:
      case Token.solvent:
      case Token.identifier:
      case Token.sidechain:
      case Token.surface:
        stack[sp++] = lookupIdentifierValue((String)instruction.value);
        break;
      case Token.opLT:
      case Token.opLE:
      case Token.opGE:
      case Token.opGT:
      case Token.opEQ:
      case Token.opNE:
        bs = stack[sp++] = new BitSet();
        comparatorInstruction(instruction, bs);
        break;
      default:
        unrecognizedExpression();
      }
    }
    if (sp != 1)
      evalError("atom expression compiler error - stack over/underflow");
    return stack[0];
  }



  void SelectionEvaluator::comparatorInstruction(Token instruction, BitSet bs) {
    int comparator = instruction.tok;
    int property = instruction.intValue;
    float propertyValue = 0; // just for temperature
    int comparisonValue = ((Integer)instruction.value).intValue();
    int numberOfAtoms = viewer.getAtomCount();
    Frame frame = viewer.getFrame();
    for (int i = 0; i < numberOfAtoms; ++i) {
      Atom atom = frame.getAtomAt(i);
      switch (property) {
      case Token.atomno:
        propertyValue = atom.getAtomNumber();
        break;
      case Token.elemno:
        propertyValue = atom.getElementNumber();
        break;
      case Token.temperature:
        propertyValue = atom.getBfactor100();
        if (propertyValue < 0)
          continue;
        propertyValue /= 100;
        break; 
      case Token._atomID:
        propertyValue = atom.getSpecialAtomID();
        if (propertyValue < 0)
          continue;
        break;
      case Token._structure:
        propertyValue = getProteinStructureType(atom);
        if (propertyValue == -1)
          continue;
        break;
      case Token.radius:
        propertyValue = atom.getRasMolRadius();
        break;
      case Token._bondedcount:
        propertyValue = atom.getCovalentBondCount();
        break;
      case Token.model:
        propertyValue = atom.getModelTagNumber();
        break;
      default:
        unrecognizedAtomProperty(property);
      }
      boolean match = false;
      switch (comparator) {
      case Token.opLT:
        match = propertyValue < comparisonValue;
        break;
      case Token.opLE:
        match = propertyValue <= comparisonValue;
        break;
      case Token.opGE:
        match = propertyValue >= comparisonValue;
        break;
      case Token.opGT:
        match = propertyValue > comparisonValue;
        break;
      case Token.opEQ:
        match = propertyValue == comparisonValue;
        break;
      case Token.opNE:
        match = propertyValue != comparisonValue;
        break;
      }
      if (match)
        bs.set(i);
    }
  }

void SelectionEvaluator::withinInstruction(Token instruction, BitSet bs, BitSet bsResult)

    Object withinSpec = instruction.value;
    if (withinSpec instanceof Float) {
        withinDistance(((Float)withinSpec).floatValue(), bs, bsResult);
        return;
    }
    evalError("Unrecognized within parameter:" + withinSpec);
}

  void SelectionEvaluator::withinDistance(float distance, BitSet bs, BitSet bsResult) {
    Frame frame = viewer.getFrame();
    for (int i = frame.getAtomCount(); --i >= 0; ) {
      if (bs.get(i)) {
        Atom atom = frame.getAtomAt(i);
        AtomIterator iterWithin =
          frame.getWithinIterator(atom, distance);
        while (iterWithin.hasNext())
          bsResult.set(iterWithin.next().getAtomIndex());
      }
    }
  }

  void SelectionEvaluator::define() throws ScriptException {
    String variable = (String)statement[1].value;
    variables.put(variable, (expression(statement, 2)));
  }
}

  void SelectionEvaluator::predefine(Token[] statement) {
    String variable = (String)statement[1].value;
    variables.put(variable, statement);
  }

  void SelectionEvaluator::select(){
    // NOTE this is called by restrict()
    if (statementLength == 1) {
      viewer.selectAll();
      if (!viewer.getRasmolHydrogenSetting())
        viewer.excludeSelectionSet(getHydrogenSet());
      if (!viewer.getRasmolHeteroSetting())
        viewer.excludeSelectionSet(getHeteroSet());
    } else {
      viewer.setSelectionSet(expression(statement, 1));
    }
    viewer.scriptStatus("" + viewer.getSelectionCount() + " atoms selected");
  }
  
  }