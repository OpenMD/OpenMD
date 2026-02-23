
#include "antlr4-runtime.h"
#include "utils/StringUtils.hpp"
#include "omdParser/FilenameObserver.hpp"


// Generated from OMD.g4 by ANTLR 4.13.2

#pragma once


#include "antlr4-runtime.h"
#include "OMDVisitor.h"


/**
 * This class provides an empty implementation of OMDVisitor, which can be
 * extended to create a visitor which only needs to handle a subset of the available methods.
 */
class  OMDBaseVisitor : public OMDVisitor {
public:

  virtual std::any visitOmdfile(OMDParser::OmdfileContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitStatement(OMDParser::StatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitAssignment(OMDParser::AssignmentContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitConstant(OMDParser::ConstantContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitComponentblock(OMDParser::ComponentblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitZconstraintblock(OMDParser::ZconstraintblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitRestraintblock(OMDParser::RestraintblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitFlucqblock(OMDParser::FlucqblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitRnemdblock(OMDParser::RnemdblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitLightblock(OMDParser::LightblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitMinimizerblock(OMDParser::MinimizerblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitMoleculeblock(OMDParser::MoleculeblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitMoleculestatement(OMDParser::MoleculestatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitAtomblock(OMDParser::AtomblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitAtomstatement(OMDParser::AtomstatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitBondblock(OMDParser::BondblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitBondstatement(OMDParser::BondstatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitBendblock(OMDParser::BendblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitBendstatement(OMDParser::BendstatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitTorsionblock(OMDParser::TorsionblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitTorsionstatement(OMDParser::TorsionstatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitInversionblock(OMDParser::InversionblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitInversionstatement(OMDParser::InversionstatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitRigidbodyblock(OMDParser::RigidbodyblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitRigidbodystatement(OMDParser::RigidbodystatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitCutoffgroupblock(OMDParser::CutoffgroupblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitCutoffgroupstatement(OMDParser::CutoffgroupstatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitNodesblock(OMDParser::NodesblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitNodesstatement(OMDParser::NodesstatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitFragmentblock(OMDParser::FragmentblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitFragmentstatement(OMDParser::FragmentstatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitConstraintblock(OMDParser::ConstraintblockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitConstraintstatement(OMDParser::ConstraintstatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitSequencestring(OMDParser::SequencestringContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitDoubleNumberTuple(OMDParser::DoubleNumberTupleContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitInttuple(OMDParser::InttupleContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitIntConst(OMDParser::IntConstContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitDoubleNumber(OMDParser::DoubleNumberContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitFloatConst(OMDParser::FloatConstContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitVectorConst(OMDParser::VectorConstContext *ctx) override {
    return visitChildren(ctx);
  }


};

