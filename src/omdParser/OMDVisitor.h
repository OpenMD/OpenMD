
#include "antlr4-runtime.h"
#include "utils/StringUtils.hpp"
#include "omdParser/FilenameObserver.hpp"


// Generated from OMD.g4 by ANTLR 4.13.2

#pragma once


#include "antlr4-runtime.h"
#include "OMDParser.h"



/**
 * This class defines an abstract visitor for a parse tree
 * produced by OMDParser.
 */
class  OMDVisitor : public antlr4::tree::AbstractParseTreeVisitor {
public:

  /**
   * Visit parse trees produced by OMDParser.
   */
    virtual std::any visitOmdfile(OMDParser::OmdfileContext *context) = 0;

    virtual std::any visitStatement(OMDParser::StatementContext *context) = 0;

    virtual std::any visitAssignment(OMDParser::AssignmentContext *context) = 0;

    virtual std::any visitConstant(OMDParser::ConstantContext *context) = 0;

    virtual std::any visitComponentblock(OMDParser::ComponentblockContext *context) = 0;

    virtual std::any visitZconstraintblock(OMDParser::ZconstraintblockContext *context) = 0;

    virtual std::any visitRestraintblock(OMDParser::RestraintblockContext *context) = 0;

    virtual std::any visitFlucqblock(OMDParser::FlucqblockContext *context) = 0;

    virtual std::any visitRnemdblock(OMDParser::RnemdblockContext *context) = 0;

    virtual std::any visitLightblock(OMDParser::LightblockContext *context) = 0;

    virtual std::any visitMinimizerblock(OMDParser::MinimizerblockContext *context) = 0;

    virtual std::any visitMoleculeblock(OMDParser::MoleculeblockContext *context) = 0;

    virtual std::any visitMoleculestatement(OMDParser::MoleculestatementContext *context) = 0;

    virtual std::any visitAtomblock(OMDParser::AtomblockContext *context) = 0;

    virtual std::any visitAtomstatement(OMDParser::AtomstatementContext *context) = 0;

    virtual std::any visitBondblock(OMDParser::BondblockContext *context) = 0;

    virtual std::any visitBondstatement(OMDParser::BondstatementContext *context) = 0;

    virtual std::any visitBendblock(OMDParser::BendblockContext *context) = 0;

    virtual std::any visitBendstatement(OMDParser::BendstatementContext *context) = 0;

    virtual std::any visitTorsionblock(OMDParser::TorsionblockContext *context) = 0;

    virtual std::any visitTorsionstatement(OMDParser::TorsionstatementContext *context) = 0;

    virtual std::any visitInversionblock(OMDParser::InversionblockContext *context) = 0;

    virtual std::any visitInversionstatement(OMDParser::InversionstatementContext *context) = 0;

    virtual std::any visitRigidbodyblock(OMDParser::RigidbodyblockContext *context) = 0;

    virtual std::any visitRigidbodystatement(OMDParser::RigidbodystatementContext *context) = 0;

    virtual std::any visitCutoffgroupblock(OMDParser::CutoffgroupblockContext *context) = 0;

    virtual std::any visitCutoffgroupstatement(OMDParser::CutoffgroupstatementContext *context) = 0;

    virtual std::any visitNodesblock(OMDParser::NodesblockContext *context) = 0;

    virtual std::any visitNodesstatement(OMDParser::NodesstatementContext *context) = 0;

    virtual std::any visitFragmentblock(OMDParser::FragmentblockContext *context) = 0;

    virtual std::any visitFragmentstatement(OMDParser::FragmentstatementContext *context) = 0;

    virtual std::any visitConstraintblock(OMDParser::ConstraintblockContext *context) = 0;

    virtual std::any visitConstraintstatement(OMDParser::ConstraintstatementContext *context) = 0;

    virtual std::any visitSequencestring(OMDParser::SequencestringContext *context) = 0;

    virtual std::any visitDoubleNumberTuple(OMDParser::DoubleNumberTupleContext *context) = 0;

    virtual std::any visitInttuple(OMDParser::InttupleContext *context) = 0;

    virtual std::any visitIntConst(OMDParser::IntConstContext *context) = 0;

    virtual std::any visitDoubleNumber(OMDParser::DoubleNumberContext *context) = 0;

    virtual std::any visitFloatConst(OMDParser::FloatConstContext *context) = 0;

    virtual std::any visitVectorConst(OMDParser::VectorConstContext *context) = 0;


};

