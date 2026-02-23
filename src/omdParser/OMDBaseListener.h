
#include "antlr4-runtime.h"
#include "utils/StringUtils.hpp"
#include "omdParser/FilenameObserver.hpp"


// Generated from OMD.g4 by ANTLR 4.13.2

#pragma once


#include "antlr4-runtime.h"
#include "OMDListener.h"


/**
 * This class provides an empty implementation of OMDListener,
 * which can be extended to create a listener which only needs to handle a subset
 * of the available methods.
 */
class  OMDBaseListener : public OMDListener {
public:

  virtual void enterOmdfile(OMDParser::OmdfileContext * /*ctx*/) override { }
  virtual void exitOmdfile(OMDParser::OmdfileContext * /*ctx*/) override { }

  virtual void enterStatement(OMDParser::StatementContext * /*ctx*/) override { }
  virtual void exitStatement(OMDParser::StatementContext * /*ctx*/) override { }

  virtual void enterAssignment(OMDParser::AssignmentContext * /*ctx*/) override { }
  virtual void exitAssignment(OMDParser::AssignmentContext * /*ctx*/) override { }

  virtual void enterConstant(OMDParser::ConstantContext * /*ctx*/) override { }
  virtual void exitConstant(OMDParser::ConstantContext * /*ctx*/) override { }

  virtual void enterComponentblock(OMDParser::ComponentblockContext * /*ctx*/) override { }
  virtual void exitComponentblock(OMDParser::ComponentblockContext * /*ctx*/) override { }

  virtual void enterZconstraintblock(OMDParser::ZconstraintblockContext * /*ctx*/) override { }
  virtual void exitZconstraintblock(OMDParser::ZconstraintblockContext * /*ctx*/) override { }

  virtual void enterRestraintblock(OMDParser::RestraintblockContext * /*ctx*/) override { }
  virtual void exitRestraintblock(OMDParser::RestraintblockContext * /*ctx*/) override { }

  virtual void enterFlucqblock(OMDParser::FlucqblockContext * /*ctx*/) override { }
  virtual void exitFlucqblock(OMDParser::FlucqblockContext * /*ctx*/) override { }

  virtual void enterRnemdblock(OMDParser::RnemdblockContext * /*ctx*/) override { }
  virtual void exitRnemdblock(OMDParser::RnemdblockContext * /*ctx*/) override { }

  virtual void enterLightblock(OMDParser::LightblockContext * /*ctx*/) override { }
  virtual void exitLightblock(OMDParser::LightblockContext * /*ctx*/) override { }

  virtual void enterMinimizerblock(OMDParser::MinimizerblockContext * /*ctx*/) override { }
  virtual void exitMinimizerblock(OMDParser::MinimizerblockContext * /*ctx*/) override { }

  virtual void enterMoleculeblock(OMDParser::MoleculeblockContext * /*ctx*/) override { }
  virtual void exitMoleculeblock(OMDParser::MoleculeblockContext * /*ctx*/) override { }

  virtual void enterMoleculestatement(OMDParser::MoleculestatementContext * /*ctx*/) override { }
  virtual void exitMoleculestatement(OMDParser::MoleculestatementContext * /*ctx*/) override { }

  virtual void enterAtomblock(OMDParser::AtomblockContext * /*ctx*/) override { }
  virtual void exitAtomblock(OMDParser::AtomblockContext * /*ctx*/) override { }

  virtual void enterAtomstatement(OMDParser::AtomstatementContext * /*ctx*/) override { }
  virtual void exitAtomstatement(OMDParser::AtomstatementContext * /*ctx*/) override { }

  virtual void enterBondblock(OMDParser::BondblockContext * /*ctx*/) override { }
  virtual void exitBondblock(OMDParser::BondblockContext * /*ctx*/) override { }

  virtual void enterBondstatement(OMDParser::BondstatementContext * /*ctx*/) override { }
  virtual void exitBondstatement(OMDParser::BondstatementContext * /*ctx*/) override { }

  virtual void enterBendblock(OMDParser::BendblockContext * /*ctx*/) override { }
  virtual void exitBendblock(OMDParser::BendblockContext * /*ctx*/) override { }

  virtual void enterBendstatement(OMDParser::BendstatementContext * /*ctx*/) override { }
  virtual void exitBendstatement(OMDParser::BendstatementContext * /*ctx*/) override { }

  virtual void enterTorsionblock(OMDParser::TorsionblockContext * /*ctx*/) override { }
  virtual void exitTorsionblock(OMDParser::TorsionblockContext * /*ctx*/) override { }

  virtual void enterTorsionstatement(OMDParser::TorsionstatementContext * /*ctx*/) override { }
  virtual void exitTorsionstatement(OMDParser::TorsionstatementContext * /*ctx*/) override { }

  virtual void enterInversionblock(OMDParser::InversionblockContext * /*ctx*/) override { }
  virtual void exitInversionblock(OMDParser::InversionblockContext * /*ctx*/) override { }

  virtual void enterInversionstatement(OMDParser::InversionstatementContext * /*ctx*/) override { }
  virtual void exitInversionstatement(OMDParser::InversionstatementContext * /*ctx*/) override { }

  virtual void enterRigidbodyblock(OMDParser::RigidbodyblockContext * /*ctx*/) override { }
  virtual void exitRigidbodyblock(OMDParser::RigidbodyblockContext * /*ctx*/) override { }

  virtual void enterRigidbodystatement(OMDParser::RigidbodystatementContext * /*ctx*/) override { }
  virtual void exitRigidbodystatement(OMDParser::RigidbodystatementContext * /*ctx*/) override { }

  virtual void enterCutoffgroupblock(OMDParser::CutoffgroupblockContext * /*ctx*/) override { }
  virtual void exitCutoffgroupblock(OMDParser::CutoffgroupblockContext * /*ctx*/) override { }

  virtual void enterCutoffgroupstatement(OMDParser::CutoffgroupstatementContext * /*ctx*/) override { }
  virtual void exitCutoffgroupstatement(OMDParser::CutoffgroupstatementContext * /*ctx*/) override { }

  virtual void enterNodesblock(OMDParser::NodesblockContext * /*ctx*/) override { }
  virtual void exitNodesblock(OMDParser::NodesblockContext * /*ctx*/) override { }

  virtual void enterNodesstatement(OMDParser::NodesstatementContext * /*ctx*/) override { }
  virtual void exitNodesstatement(OMDParser::NodesstatementContext * /*ctx*/) override { }

  virtual void enterFragmentblock(OMDParser::FragmentblockContext * /*ctx*/) override { }
  virtual void exitFragmentblock(OMDParser::FragmentblockContext * /*ctx*/) override { }

  virtual void enterFragmentstatement(OMDParser::FragmentstatementContext * /*ctx*/) override { }
  virtual void exitFragmentstatement(OMDParser::FragmentstatementContext * /*ctx*/) override { }

  virtual void enterConstraintblock(OMDParser::ConstraintblockContext * /*ctx*/) override { }
  virtual void exitConstraintblock(OMDParser::ConstraintblockContext * /*ctx*/) override { }

  virtual void enterConstraintstatement(OMDParser::ConstraintstatementContext * /*ctx*/) override { }
  virtual void exitConstraintstatement(OMDParser::ConstraintstatementContext * /*ctx*/) override { }

  virtual void enterSequencestring(OMDParser::SequencestringContext * /*ctx*/) override { }
  virtual void exitSequencestring(OMDParser::SequencestringContext * /*ctx*/) override { }

  virtual void enterDoubleNumberTuple(OMDParser::DoubleNumberTupleContext * /*ctx*/) override { }
  virtual void exitDoubleNumberTuple(OMDParser::DoubleNumberTupleContext * /*ctx*/) override { }

  virtual void enterInttuple(OMDParser::InttupleContext * /*ctx*/) override { }
  virtual void exitInttuple(OMDParser::InttupleContext * /*ctx*/) override { }

  virtual void enterIntConst(OMDParser::IntConstContext * /*ctx*/) override { }
  virtual void exitIntConst(OMDParser::IntConstContext * /*ctx*/) override { }

  virtual void enterDoubleNumber(OMDParser::DoubleNumberContext * /*ctx*/) override { }
  virtual void exitDoubleNumber(OMDParser::DoubleNumberContext * /*ctx*/) override { }

  virtual void enterFloatConst(OMDParser::FloatConstContext * /*ctx*/) override { }
  virtual void exitFloatConst(OMDParser::FloatConstContext * /*ctx*/) override { }

  virtual void enterVectorConst(OMDParser::VectorConstContext * /*ctx*/) override { }
  virtual void exitVectorConst(OMDParser::VectorConstContext * /*ctx*/) override { }


  virtual void enterEveryRule(antlr4::ParserRuleContext * /*ctx*/) override { }
  virtual void exitEveryRule(antlr4::ParserRuleContext * /*ctx*/) override { }
  virtual void visitTerminal(antlr4::tree::TerminalNode * /*node*/) override { }
  virtual void visitErrorNode(antlr4::tree::ErrorNode * /*node*/) override { }

};

