
#include "antlr4-runtime.h"
#include "utils/StringUtils.hpp"
#include "omdParser/FilenameObserver.hpp"


// Generated from OMD.g4 by ANTLR 4.13.2

#pragma once


#include "antlr4-runtime.h"
#include "OMDParser.h"


/**
 * This interface defines an abstract listener for a parse tree produced by OMDParser.
 */
class  OMDListener : public antlr4::tree::ParseTreeListener {
public:

  virtual void enterOmdfile(OMDParser::OmdfileContext *ctx) = 0;
  virtual void exitOmdfile(OMDParser::OmdfileContext *ctx) = 0;

  virtual void enterStatement(OMDParser::StatementContext *ctx) = 0;
  virtual void exitStatement(OMDParser::StatementContext *ctx) = 0;

  virtual void enterAssignment(OMDParser::AssignmentContext *ctx) = 0;
  virtual void exitAssignment(OMDParser::AssignmentContext *ctx) = 0;

  virtual void enterConstant(OMDParser::ConstantContext *ctx) = 0;
  virtual void exitConstant(OMDParser::ConstantContext *ctx) = 0;

  virtual void enterComponentblock(OMDParser::ComponentblockContext *ctx) = 0;
  virtual void exitComponentblock(OMDParser::ComponentblockContext *ctx) = 0;

  virtual void enterZconstraintblock(OMDParser::ZconstraintblockContext *ctx) = 0;
  virtual void exitZconstraintblock(OMDParser::ZconstraintblockContext *ctx) = 0;

  virtual void enterRestraintblock(OMDParser::RestraintblockContext *ctx) = 0;
  virtual void exitRestraintblock(OMDParser::RestraintblockContext *ctx) = 0;

  virtual void enterFlucqblock(OMDParser::FlucqblockContext *ctx) = 0;
  virtual void exitFlucqblock(OMDParser::FlucqblockContext *ctx) = 0;

  virtual void enterRnemdblock(OMDParser::RnemdblockContext *ctx) = 0;
  virtual void exitRnemdblock(OMDParser::RnemdblockContext *ctx) = 0;

  virtual void enterLightblock(OMDParser::LightblockContext *ctx) = 0;
  virtual void exitLightblock(OMDParser::LightblockContext *ctx) = 0;

  virtual void enterMinimizerblock(OMDParser::MinimizerblockContext *ctx) = 0;
  virtual void exitMinimizerblock(OMDParser::MinimizerblockContext *ctx) = 0;

  virtual void enterMoleculeblock(OMDParser::MoleculeblockContext *ctx) = 0;
  virtual void exitMoleculeblock(OMDParser::MoleculeblockContext *ctx) = 0;

  virtual void enterMoleculestatement(OMDParser::MoleculestatementContext *ctx) = 0;
  virtual void exitMoleculestatement(OMDParser::MoleculestatementContext *ctx) = 0;

  virtual void enterAtomblock(OMDParser::AtomblockContext *ctx) = 0;
  virtual void exitAtomblock(OMDParser::AtomblockContext *ctx) = 0;

  virtual void enterAtomstatement(OMDParser::AtomstatementContext *ctx) = 0;
  virtual void exitAtomstatement(OMDParser::AtomstatementContext *ctx) = 0;

  virtual void enterBondblock(OMDParser::BondblockContext *ctx) = 0;
  virtual void exitBondblock(OMDParser::BondblockContext *ctx) = 0;

  virtual void enterBondstatement(OMDParser::BondstatementContext *ctx) = 0;
  virtual void exitBondstatement(OMDParser::BondstatementContext *ctx) = 0;

  virtual void enterBendblock(OMDParser::BendblockContext *ctx) = 0;
  virtual void exitBendblock(OMDParser::BendblockContext *ctx) = 0;

  virtual void enterBendstatement(OMDParser::BendstatementContext *ctx) = 0;
  virtual void exitBendstatement(OMDParser::BendstatementContext *ctx) = 0;

  virtual void enterTorsionblock(OMDParser::TorsionblockContext *ctx) = 0;
  virtual void exitTorsionblock(OMDParser::TorsionblockContext *ctx) = 0;

  virtual void enterTorsionstatement(OMDParser::TorsionstatementContext *ctx) = 0;
  virtual void exitTorsionstatement(OMDParser::TorsionstatementContext *ctx) = 0;

  virtual void enterInversionblock(OMDParser::InversionblockContext *ctx) = 0;
  virtual void exitInversionblock(OMDParser::InversionblockContext *ctx) = 0;

  virtual void enterInversionstatement(OMDParser::InversionstatementContext *ctx) = 0;
  virtual void exitInversionstatement(OMDParser::InversionstatementContext *ctx) = 0;

  virtual void enterRigidbodyblock(OMDParser::RigidbodyblockContext *ctx) = 0;
  virtual void exitRigidbodyblock(OMDParser::RigidbodyblockContext *ctx) = 0;

  virtual void enterRigidbodystatement(OMDParser::RigidbodystatementContext *ctx) = 0;
  virtual void exitRigidbodystatement(OMDParser::RigidbodystatementContext *ctx) = 0;

  virtual void enterCutoffgroupblock(OMDParser::CutoffgroupblockContext *ctx) = 0;
  virtual void exitCutoffgroupblock(OMDParser::CutoffgroupblockContext *ctx) = 0;

  virtual void enterCutoffgroupstatement(OMDParser::CutoffgroupstatementContext *ctx) = 0;
  virtual void exitCutoffgroupstatement(OMDParser::CutoffgroupstatementContext *ctx) = 0;

  virtual void enterNodesblock(OMDParser::NodesblockContext *ctx) = 0;
  virtual void exitNodesblock(OMDParser::NodesblockContext *ctx) = 0;

  virtual void enterNodesstatement(OMDParser::NodesstatementContext *ctx) = 0;
  virtual void exitNodesstatement(OMDParser::NodesstatementContext *ctx) = 0;

  virtual void enterFragmentblock(OMDParser::FragmentblockContext *ctx) = 0;
  virtual void exitFragmentblock(OMDParser::FragmentblockContext *ctx) = 0;

  virtual void enterFragmentstatement(OMDParser::FragmentstatementContext *ctx) = 0;
  virtual void exitFragmentstatement(OMDParser::FragmentstatementContext *ctx) = 0;

  virtual void enterConstraintblock(OMDParser::ConstraintblockContext *ctx) = 0;
  virtual void exitConstraintblock(OMDParser::ConstraintblockContext *ctx) = 0;

  virtual void enterConstraintstatement(OMDParser::ConstraintstatementContext *ctx) = 0;
  virtual void exitConstraintstatement(OMDParser::ConstraintstatementContext *ctx) = 0;

  virtual void enterSequencestring(OMDParser::SequencestringContext *ctx) = 0;
  virtual void exitSequencestring(OMDParser::SequencestringContext *ctx) = 0;

  virtual void enterDoubleNumberTuple(OMDParser::DoubleNumberTupleContext *ctx) = 0;
  virtual void exitDoubleNumberTuple(OMDParser::DoubleNumberTupleContext *ctx) = 0;

  virtual void enterInttuple(OMDParser::InttupleContext *ctx) = 0;
  virtual void exitInttuple(OMDParser::InttupleContext *ctx) = 0;

  virtual void enterIntConst(OMDParser::IntConstContext *ctx) = 0;
  virtual void exitIntConst(OMDParser::IntConstContext *ctx) = 0;

  virtual void enterDoubleNumber(OMDParser::DoubleNumberContext *ctx) = 0;
  virtual void exitDoubleNumber(OMDParser::DoubleNumberContext *ctx) = 0;

  virtual void enterFloatConst(OMDParser::FloatConstContext *ctx) = 0;
  virtual void exitFloatConst(OMDParser::FloatConstContext *ctx) = 0;

  virtual void enterVectorConst(OMDParser::VectorConstContext *ctx) = 0;
  virtual void exitVectorConst(OMDParser::VectorConstContext *ctx) = 0;


};

