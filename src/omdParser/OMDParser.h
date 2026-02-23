
#include "antlr4-runtime.h"
#include "utils/StringUtils.hpp"
#include "omdParser/FilenameObserver.hpp"


// Generated from OMD.g4 by ANTLR 4.13.2

#pragma once


#include "antlr4-runtime.h"




class  OMDParser : public antlr4::Parser {
public:
  enum {
    COMPONENT = 1, MOLECULE = 2, ZCONSTRAINT = 3, RESTRAINT = 4, ATOM = 5, 
    BOND = 6, BEND = 7, TORSION = 8, INVERSION = 9, RIGIDBODY = 10, CUTOFFGROUP = 11, 
    CONSTRAINT = 12, DISTANCE = 13, FRAGMENT = 14, SEQUENCE = 15, MEMBERS = 16, 
    CENTER = 17, SATELLITES = 18, POSITION = 19, ORIENTATION = 20, FLUCQ = 21, 
    RNEMD = 22, LIGHT = 23, MINIMIZER = 24, FIXED = 25, HARMONIC = 26, CUBIC = 27, 
    QUARTIC = 28, POLYNOMIAL = 29, MORSE = 30, GHOSTBEND = 31, UREYBRADLEY = 32, 
    COSINE = 33, GHOSTTORSION = 34, CHARMM = 35, OPLS = 36, TRAPPE = 37, 
    AMBERIMPROPER = 38, IMPROPERCOSINE = 39, CENTRALATOMHEIGHT = 40, DREIDING = 41, 
    CHARGE = 42, NODES = 43, ASSIGNEQUAL = 44, COLON = 45, COMMA = 46, QUESTIONMARK = 47, 
    SEMICOLON = 48, DOT = 49, LPAREN = 50, RPAREN = 51, LBRACKET = 52, RBRACKET = 53, 
    LCURLY = 54, RCURLY = 55, ID = 56, NUM_LONG = 57, NUM_INT = 58, NUM_FLOAT = 59, 
    NUM_DOUBLE = 60, CharLiteral = 61, StringLiteral = 62, Whitespace = 63, 
    Newline = 64, LineContinuation = 65, Comment = 66, CPPComment = 67, 
    PREPROC_DIRECTIVE = 68
  };

  enum {
    RuleOmdfile = 0, RuleStatement = 1, RuleAssignment = 2, RuleConstant = 3, 
    RuleComponentblock = 4, RuleZconstraintblock = 5, RuleRestraintblock = 6, 
    RuleFlucqblock = 7, RuleRnemdblock = 8, RuleLightblock = 9, RuleMinimizerblock = 10, 
    RuleMoleculeblock = 11, RuleMoleculestatement = 12, RuleAtomblock = 13, 
    RuleAtomstatement = 14, RuleBondblock = 15, RuleBondstatement = 16, 
    RuleBendblock = 17, RuleBendstatement = 18, RuleTorsionblock = 19, RuleTorsionstatement = 20, 
    RuleInversionblock = 21, RuleInversionstatement = 22, RuleRigidbodyblock = 23, 
    RuleRigidbodystatement = 24, RuleCutoffgroupblock = 25, RuleCutoffgroupstatement = 26, 
    RuleNodesblock = 27, RuleNodesstatement = 28, RuleFragmentblock = 29, 
    RuleFragmentstatement = 30, RuleConstraintblock = 31, RuleConstraintstatement = 32, 
    RuleSequencestring = 33, RuleDoubleNumberTuple = 34, RuleInttuple = 35, 
    RuleIntConst = 36, RuleDoubleNumber = 37, RuleFloatConst = 38, RuleVectorConst = 39
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

