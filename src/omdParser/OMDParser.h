
#include "antlr4-runtime.h"
#include "utils/StringUtils.hpp"
#include "omdParser/FilenameObserver.hpp"


// Generated from OMD.g4 by ANTLR 4.13.2

#pragma once


#include "antlr4-runtime.h"




class  OMDParser : public antlr4::Parser {
public:
  enum {
    TRUE = 1, FALSE = 2, COMPONENT = 3, MOLECULE = 4, ZCONSTRAINT = 5, RESTRAINT = 6, 
    ATOM = 7, BOND = 8, BEND = 9, TORSION = 10, INVERSION = 11, RIGIDBODY = 12, 
    CUTOFFGROUP = 13, CONSTRAINT = 14, DISTANCE = 15, FRAGMENT = 16, SEQUENCE = 17, 
    MEMBERS = 18, CENTER = 19, SATELLITES = 20, POSITION = 21, ORIENTATION = 22, 
    FLUCQ = 23, RNEMD = 24, LIGHT = 25, VELOCITYFIELD = 26, MINIMIZER = 27, 
    FIXED = 28, HARMONIC = 29, CUBIC = 30, QUARTIC = 31, POLYNOMIAL = 32, 
    MORSE = 33, GHOSTBEND = 34, UREYBRADLEY = 35, COSINE = 36, GHOSTTORSION = 37, 
    CHARMM = 38, OPLS = 39, TRAPPE = 40, AMBERIMPROPER = 41, IMPROPERCOSINE = 42, 
    CENTRALATOMHEIGHT = 43, DREIDING = 44, CHARGE = 45, NODES = 46, ASSIGNEQUAL = 47, 
    COLON = 48, COMMA = 49, QUESTIONMARK = 50, SEMICOLON = 51, DOT = 52, 
    LPAREN = 53, RPAREN = 54, LBRACKET = 55, RBRACKET = 56, LCURLY = 57, 
    RCURLY = 58, NUM_LONG = 59, NUM_INT = 60, NUM_FLOAT = 61, NUM_DOUBLE = 62, 
    CharLiteral = 63, StringLiteral = 64, ID = 65, Whitespace = 66, Newline = 67, 
    LineContinuation = 68, Comment = 69, CPPComment = 70, PREPROC_DIRECTIVE = 71
  };

  enum {
    RuleOmdfile = 0, RuleStatement = 1, RuleAssignment = 2, RuleConstant = 3, 
    RuleComponentblock = 4, RuleZconstraintblock = 5, RuleRestraintblock = 6, 
    RuleFlucqblock = 7, RuleRnemdblock = 8, RuleLightblock = 9, RuleVelocityfieldblock = 10, 
    RuleMinimizerblock = 11, RuleMoleculeblock = 12, RuleMoleculestatement = 13, 
    RuleAtomblock = 14, RuleAtomstatement = 15, RuleBondblock = 16, RuleBondstatement = 17, 
    RuleBendblock = 18, RuleBendstatement = 19, RuleTorsionblock = 20, RuleTorsionstatement = 21, 
    RuleInversionblock = 22, RuleInversionstatement = 23, RuleRigidbodyblock = 24, 
    RuleRigidbodystatement = 25, RuleCutoffgroupblock = 26, RuleCutoffgroupstatement = 27, 
    RuleNodesblock = 28, RuleNodesstatement = 29, RuleFragmentblock = 30, 
    RuleFragmentstatement = 31, RuleConstraintblock = 32, RuleConstraintstatement = 33, 
    RuleSequencestring = 34, RuleDoubleNumberTuple = 35, RuleInttuple = 36, 
    RuleIntConst = 37, RuleDoubleNumber = 38, RuleFloatConst = 39, RuleVectorConst = 40
  };

  explicit OMDParser(antlr4::TokenStream *input);

  OMDParser(antlr4::TokenStream *input, const antlr4::atn::ParserATNSimulatorOptions &options);

  ~OMDParser() override;

  std::string getGrammarFileName() const override;

  const antlr4::atn::ATN& getATN() const override;

  const std::vector<std::string>& getRuleNames() const override;

  const antlr4::dfa::Vocabulary& getVocabulary() const override;

  antlr4::atn::SerializedATNView getSerializedATN() const override;


      FilenameObserver* observer;

      void setObserver(FilenameObserver* osv) {
          observer = osv;
      }


  class OmdfileContext;
  class StatementContext;
  class AssignmentContext;
  class ConstantContext;
  class ComponentblockContext;
  class ZconstraintblockContext;
  class RestraintblockContext;
  class FlucqblockContext;
  class RnemdblockContext;
  class LightblockContext;
  class VelocityfieldblockContext;
  class MinimizerblockContext;
  class MoleculeblockContext;
  class MoleculestatementContext;
  class AtomblockContext;
  class AtomstatementContext;
  class BondblockContext;
  class BondstatementContext;
  class BendblockContext;
  class BendstatementContext;
  class TorsionblockContext;
  class TorsionstatementContext;
  class InversionblockContext;
  class InversionstatementContext;
  class RigidbodyblockContext;
  class RigidbodystatementContext;
  class CutoffgroupblockContext;
  class CutoffgroupstatementContext;
  class NodesblockContext;
  class NodesstatementContext;
  class FragmentblockContext;
  class FragmentstatementContext;
  class ConstraintblockContext;
  class ConstraintstatementContext;
  class SequencestringContext;
  class DoubleNumberTupleContext;
  class InttupleContext;
  class IntConstContext;
  class DoubleNumberContext;
  class FloatConstContext;
  class VectorConstContext; 

  class  OmdfileContext : public antlr4::ParserRuleContext {
  public:
    OmdfileContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *EOF();
    std::vector<StatementContext *> statement();
    StatementContext* statement(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  OmdfileContext* omdfile();

  class  StatementContext : public antlr4::ParserRuleContext {
  public:
    StatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    AssignmentContext *assignment();
    ComponentblockContext *componentblock();
    MoleculeblockContext *moleculeblock();
    FragmentblockContext *fragmentblock();
    ZconstraintblockContext *zconstraintblock();
    RestraintblockContext *restraintblock();
    FlucqblockContext *flucqblock();
    RnemdblockContext *rnemdblock();
    LightblockContext *lightblock();
    VelocityfieldblockContext *velocityfieldblock();
    MinimizerblockContext *minimizerblock();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  StatementContext* statement();

  class  AssignmentContext : public antlr4::ParserRuleContext {
  public:
    AssignmentContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *ID();
    antlr4::tree::TerminalNode *ASSIGNEQUAL();
    ConstantContext *constant();
    antlr4::tree::TerminalNode *SEMICOLON();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  AssignmentContext* assignment();

  class  ConstantContext : public antlr4::ParserRuleContext {
  public:
    ConstantContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    IntConstContext *intConst();
    FloatConstContext *floatConst();
    VectorConstContext *vectorConst();
    antlr4::tree::TerminalNode *TRUE();
    antlr4::tree::TerminalNode *FALSE();
    antlr4::tree::TerminalNode *ID();
    antlr4::tree::TerminalNode *StringLiteral();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ConstantContext* constant();

  class  ComponentblockContext : public antlr4::ParserRuleContext {
  public:
    ComponentblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *COMPONENT();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    std::vector<AssignmentContext *> assignment();
    AssignmentContext* assignment(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ComponentblockContext* componentblock();

  class  ZconstraintblockContext : public antlr4::ParserRuleContext {
  public:
    ZconstraintblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *ZCONSTRAINT();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    std::vector<AssignmentContext *> assignment();
    AssignmentContext* assignment(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ZconstraintblockContext* zconstraintblock();

  class  RestraintblockContext : public antlr4::ParserRuleContext {
  public:
    RestraintblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *RESTRAINT();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    std::vector<AssignmentContext *> assignment();
    AssignmentContext* assignment(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  RestraintblockContext* restraintblock();

  class  FlucqblockContext : public antlr4::ParserRuleContext {
  public:
    FlucqblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *FLUCQ();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    std::vector<AssignmentContext *> assignment();
    AssignmentContext* assignment(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  FlucqblockContext* flucqblock();

  class  RnemdblockContext : public antlr4::ParserRuleContext {
  public:
    RnemdblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *RNEMD();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    std::vector<AssignmentContext *> assignment();
    AssignmentContext* assignment(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  RnemdblockContext* rnemdblock();

  class  LightblockContext : public antlr4::ParserRuleContext {
  public:
    LightblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *LIGHT();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    std::vector<AssignmentContext *> assignment();
    AssignmentContext* assignment(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  LightblockContext* lightblock();

  class  VelocityfieldblockContext : public antlr4::ParserRuleContext {
  public:
    VelocityfieldblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *VELOCITYFIELD();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    std::vector<AssignmentContext *> assignment();
    AssignmentContext* assignment(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  VelocityfieldblockContext* velocityfieldblock();

  class  MinimizerblockContext : public antlr4::ParserRuleContext {
  public:
    MinimizerblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *MINIMIZER();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    std::vector<AssignmentContext *> assignment();
    AssignmentContext* assignment(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  MinimizerblockContext* minimizerblock();

  class  MoleculeblockContext : public antlr4::ParserRuleContext {
  public:
    MoleculeblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *MOLECULE();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    std::vector<MoleculestatementContext *> moleculestatement();
    MoleculestatementContext* moleculestatement(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  MoleculeblockContext* moleculeblock();

  class  MoleculestatementContext : public antlr4::ParserRuleContext {
  public:
    MoleculestatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    AssignmentContext *assignment();
    AtomblockContext *atomblock();
    BondblockContext *bondblock();
    BendblockContext *bendblock();
    TorsionblockContext *torsionblock();
    InversionblockContext *inversionblock();
    RigidbodyblockContext *rigidbodyblock();
    CutoffgroupblockContext *cutoffgroupblock();
    ConstraintblockContext *constraintblock();
    SequencestringContext *sequencestring();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  MoleculestatementContext* moleculestatement();

  class  AtomblockContext : public antlr4::ParserRuleContext {
  public:
    AtomblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *ATOM();
    antlr4::tree::TerminalNode *LBRACKET();
    IntConstContext *intConst();
    antlr4::tree::TerminalNode *RBRACKET();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    std::vector<AtomstatementContext *> atomstatement();
    AtomstatementContext* atomstatement(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  AtomblockContext* atomblock();

  class  AtomstatementContext : public antlr4::ParserRuleContext {
  public:
    AtomstatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    AssignmentContext *assignment();
    antlr4::tree::TerminalNode *POSITION();
    antlr4::tree::TerminalNode *LPAREN();
    DoubleNumberTupleContext *doubleNumberTuple();
    antlr4::tree::TerminalNode *RPAREN();
    antlr4::tree::TerminalNode *SEMICOLON();
    antlr4::tree::TerminalNode *ORIENTATION();
    antlr4::tree::TerminalNode *CHARGE();
    FloatConstContext *floatConst();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  AtomstatementContext* atomstatement();

  class  BondblockContext : public antlr4::ParserRuleContext {
  public:
    BondblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *BOND();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    antlr4::tree::TerminalNode *LBRACKET();
    IntConstContext *intConst();
    antlr4::tree::TerminalNode *RBRACKET();
    std::vector<BondstatementContext *> bondstatement();
    BondstatementContext* bondstatement(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  BondblockContext* bondblock();

  class  BondstatementContext : public antlr4::ParserRuleContext {
  public:
    BondstatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    AssignmentContext *assignment();
    antlr4::tree::TerminalNode *MEMBERS();
    antlr4::tree::TerminalNode *LPAREN();
    InttupleContext *inttuple();
    antlr4::tree::TerminalNode *RPAREN();
    antlr4::tree::TerminalNode *SEMICOLON();
    antlr4::tree::TerminalNode *FIXED();
    FloatConstContext *floatConst();
    antlr4::tree::TerminalNode *HARMONIC();
    DoubleNumberTupleContext *doubleNumberTuple();
    antlr4::tree::TerminalNode *CUBIC();
    antlr4::tree::TerminalNode *QUARTIC();
    antlr4::tree::TerminalNode *POLYNOMIAL();
    antlr4::tree::TerminalNode *MORSE();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  BondstatementContext* bondstatement();

  class  BendblockContext : public antlr4::ParserRuleContext {
  public:
    BendblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *BEND();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    antlr4::tree::TerminalNode *LBRACKET();
    IntConstContext *intConst();
    antlr4::tree::TerminalNode *RBRACKET();
    std::vector<BendstatementContext *> bendstatement();
    BendstatementContext* bendstatement(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  BendblockContext* bendblock();

  class  BendstatementContext : public antlr4::ParserRuleContext {
  public:
    BendstatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    AssignmentContext *assignment();
    antlr4::tree::TerminalNode *MEMBERS();
    antlr4::tree::TerminalNode *LPAREN();
    InttupleContext *inttuple();
    antlr4::tree::TerminalNode *RPAREN();
    antlr4::tree::TerminalNode *SEMICOLON();
    antlr4::tree::TerminalNode *HARMONIC();
    DoubleNumberTupleContext *doubleNumberTuple();
    antlr4::tree::TerminalNode *GHOSTBEND();
    antlr4::tree::TerminalNode *UREYBRADLEY();
    antlr4::tree::TerminalNode *CUBIC();
    antlr4::tree::TerminalNode *QUARTIC();
    antlr4::tree::TerminalNode *POLYNOMIAL();
    antlr4::tree::TerminalNode *COSINE();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  BendstatementContext* bendstatement();

  class  TorsionblockContext : public antlr4::ParserRuleContext {
  public:
    TorsionblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *TORSION();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    antlr4::tree::TerminalNode *LBRACKET();
    IntConstContext *intConst();
    antlr4::tree::TerminalNode *RBRACKET();
    std::vector<TorsionstatementContext *> torsionstatement();
    TorsionstatementContext* torsionstatement(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  TorsionblockContext* torsionblock();

  class  TorsionstatementContext : public antlr4::ParserRuleContext {
  public:
    TorsionstatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    AssignmentContext *assignment();
    antlr4::tree::TerminalNode *MEMBERS();
    antlr4::tree::TerminalNode *LPAREN();
    InttupleContext *inttuple();
    antlr4::tree::TerminalNode *RPAREN();
    antlr4::tree::TerminalNode *SEMICOLON();
    antlr4::tree::TerminalNode *GHOSTTORSION();
    DoubleNumberTupleContext *doubleNumberTuple();
    antlr4::tree::TerminalNode *CUBIC();
    antlr4::tree::TerminalNode *QUARTIC();
    antlr4::tree::TerminalNode *POLYNOMIAL();
    antlr4::tree::TerminalNode *CHARMM();
    antlr4::tree::TerminalNode *OPLS();
    antlr4::tree::TerminalNode *TRAPPE();
    antlr4::tree::TerminalNode *HARMONIC();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  TorsionstatementContext* torsionstatement();

  class  InversionblockContext : public antlr4::ParserRuleContext {
  public:
    InversionblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *INVERSION();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    antlr4::tree::TerminalNode *LBRACKET();
    IntConstContext *intConst();
    antlr4::tree::TerminalNode *RBRACKET();
    std::vector<InversionstatementContext *> inversionstatement();
    InversionstatementContext* inversionstatement(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  InversionblockContext* inversionblock();

  class  InversionstatementContext : public antlr4::ParserRuleContext {
  public:
    InversionstatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    AssignmentContext *assignment();
    antlr4::tree::TerminalNode *CENTER();
    antlr4::tree::TerminalNode *LPAREN();
    IntConstContext *intConst();
    antlr4::tree::TerminalNode *RPAREN();
    antlr4::tree::TerminalNode *SEMICOLON();
    antlr4::tree::TerminalNode *SATELLITES();
    InttupleContext *inttuple();
    antlr4::tree::TerminalNode *AMBERIMPROPER();
    DoubleNumberTupleContext *doubleNumberTuple();
    antlr4::tree::TerminalNode *IMPROPERCOSINE();
    antlr4::tree::TerminalNode *HARMONIC();
    antlr4::tree::TerminalNode *CENTRALATOMHEIGHT();
    antlr4::tree::TerminalNode *DREIDING();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  InversionstatementContext* inversionstatement();

  class  RigidbodyblockContext : public antlr4::ParserRuleContext {
  public:
    RigidbodyblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *RIGIDBODY();
    antlr4::tree::TerminalNode *LBRACKET();
    IntConstContext *intConst();
    antlr4::tree::TerminalNode *RBRACKET();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    std::vector<RigidbodystatementContext *> rigidbodystatement();
    RigidbodystatementContext* rigidbodystatement(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  RigidbodyblockContext* rigidbodyblock();

  class  RigidbodystatementContext : public antlr4::ParserRuleContext {
  public:
    RigidbodystatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    AssignmentContext *assignment();
    antlr4::tree::TerminalNode *MEMBERS();
    antlr4::tree::TerminalNode *LPAREN();
    InttupleContext *inttuple();
    antlr4::tree::TerminalNode *RPAREN();
    antlr4::tree::TerminalNode *SEMICOLON();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  RigidbodystatementContext* rigidbodystatement();

  class  CutoffgroupblockContext : public antlr4::ParserRuleContext {
  public:
    CutoffgroupblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *CUTOFFGROUP();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    antlr4::tree::TerminalNode *LBRACKET();
    IntConstContext *intConst();
    antlr4::tree::TerminalNode *RBRACKET();
    std::vector<CutoffgroupstatementContext *> cutoffgroupstatement();
    CutoffgroupstatementContext* cutoffgroupstatement(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  CutoffgroupblockContext* cutoffgroupblock();

  class  CutoffgroupstatementContext : public antlr4::ParserRuleContext {
  public:
    CutoffgroupstatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    AssignmentContext *assignment();
    antlr4::tree::TerminalNode *MEMBERS();
    antlr4::tree::TerminalNode *LPAREN();
    InttupleContext *inttuple();
    antlr4::tree::TerminalNode *RPAREN();
    antlr4::tree::TerminalNode *SEMICOLON();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  CutoffgroupstatementContext* cutoffgroupstatement();

  class  NodesblockContext : public antlr4::ParserRuleContext {
  public:
    NodesblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *NODES();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    antlr4::tree::TerminalNode *LBRACKET();
    IntConstContext *intConst();
    antlr4::tree::TerminalNode *RBRACKET();
    std::vector<NodesstatementContext *> nodesstatement();
    NodesstatementContext* nodesstatement(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  NodesblockContext* nodesblock();

  class  NodesstatementContext : public antlr4::ParserRuleContext {
  public:
    NodesstatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    AssignmentContext *assignment();
    antlr4::tree::TerminalNode *MEMBERS();
    antlr4::tree::TerminalNode *LPAREN();
    InttupleContext *inttuple();
    antlr4::tree::TerminalNode *RPAREN();
    antlr4::tree::TerminalNode *SEMICOLON();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  NodesstatementContext* nodesstatement();

  class  FragmentblockContext : public antlr4::ParserRuleContext {
  public:
    FragmentblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *FRAGMENT();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    std::vector<FragmentstatementContext *> fragmentstatement();
    FragmentstatementContext* fragmentstatement(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  FragmentblockContext* fragmentblock();

  class  FragmentstatementContext : public antlr4::ParserRuleContext {
  public:
    FragmentstatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    AssignmentContext *assignment();
    AtomblockContext *atomblock();
    BondblockContext *bondblock();
    BendblockContext *bendblock();
    TorsionblockContext *torsionblock();
    InversionblockContext *inversionblock();
    RigidbodyblockContext *rigidbodyblock();
    CutoffgroupblockContext *cutoffgroupblock();
    ConstraintblockContext *constraintblock();
    NodesblockContext *nodesblock();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  FragmentstatementContext* fragmentstatement();

  class  ConstraintblockContext : public antlr4::ParserRuleContext {
  public:
    ConstraintblockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *CONSTRAINT();
    antlr4::tree::TerminalNode *LCURLY();
    antlr4::tree::TerminalNode *RCURLY();
    antlr4::tree::TerminalNode *LBRACKET();
    IntConstContext *intConst();
    antlr4::tree::TerminalNode *RBRACKET();
    std::vector<ConstraintstatementContext *> constraintstatement();
    ConstraintstatementContext* constraintstatement(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ConstraintblockContext* constraintblock();

  class  ConstraintstatementContext : public antlr4::ParserRuleContext {
  public:
    ConstraintstatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    AssignmentContext *assignment();
    antlr4::tree::TerminalNode *MEMBERS();
    antlr4::tree::TerminalNode *LPAREN();
    InttupleContext *inttuple();
    antlr4::tree::TerminalNode *RPAREN();
    antlr4::tree::TerminalNode *SEMICOLON();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ConstraintstatementContext* constraintstatement();

  class  SequencestringContext : public antlr4::ParserRuleContext {
  public:
    SequencestringContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *SEQUENCE();
    antlr4::tree::TerminalNode *ASSIGNEQUAL();
    ConstantContext *constant();
    antlr4::tree::TerminalNode *SEMICOLON();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  SequencestringContext* sequencestring();

  class  DoubleNumberTupleContext : public antlr4::ParserRuleContext {
  public:
    DoubleNumberTupleContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<DoubleNumberContext *> doubleNumber();
    DoubleNumberContext* doubleNumber(size_t i);
    std::vector<antlr4::tree::TerminalNode *> COMMA();
    antlr4::tree::TerminalNode* COMMA(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  DoubleNumberTupleContext* doubleNumberTuple();

  class  InttupleContext : public antlr4::ParserRuleContext {
  public:
    InttupleContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<IntConstContext *> intConst();
    IntConstContext* intConst(size_t i);
    std::vector<antlr4::tree::TerminalNode *> COMMA();
    antlr4::tree::TerminalNode* COMMA(size_t i);

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  InttupleContext* inttuple();

  class  IntConstContext : public antlr4::ParserRuleContext {
  public:
    IntConstContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *NUM_INT();
    antlr4::tree::TerminalNode *NUM_LONG();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  IntConstContext* intConst();

  class  DoubleNumberContext : public antlr4::ParserRuleContext {
  public:
    DoubleNumberContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    IntConstContext *intConst();
    FloatConstContext *floatConst();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  DoubleNumberContext* doubleNumber();

  class  FloatConstContext : public antlr4::ParserRuleContext {
  public:
    FloatConstContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *NUM_FLOAT();
    antlr4::tree::TerminalNode *NUM_DOUBLE();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  FloatConstContext* floatConst();

  class  VectorConstContext : public antlr4::ParserRuleContext {
  public:
    VectorConstContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *LPAREN();
    DoubleNumberTupleContext *doubleNumberTuple();
    antlr4::tree::TerminalNode *RPAREN();

    virtual void enterRule(antlr4::tree::ParseTreeListener *listener) override;
    virtual void exitRule(antlr4::tree::ParseTreeListener *listener) override;

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  VectorConstContext* vectorConst();


  // By default the static state used to implement the parser is lazily initialized during the first
  // call to the constructor. You can call this function if you wish to initialize the static state
  // ahead of time.
  static void initialize();

private:
};

