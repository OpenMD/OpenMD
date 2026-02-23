
// ANTLR v4 Visitor for OpenMD Tree Parsing
// Converted from ANTLR v2 MDTreeParser.g
// This uses the Visitor pattern instead of tree walking

#ifndef OMDTREEVISITOR_HPP
#define OMDTREEVISITOR_HPP

#include <stack>
#include <vector>
#include "io/Globals.hpp"
#include "utils/StringUtils.hpp"
#include "OMDBaseVisitor.h"

using namespace std;
using namespace OpenMD;

class OMDTreeVisitor : public OMDBaseVisitor {
private:
  Globals* currConf;
  stack<DataHolder*> blockStack;

public:
  OMDTreeVisitor() : currConf(nullptr) {}
    
  Globals* walkTree(antlr4::tree::ParseTree* tree) {
    currConf = new Globals;
    blockStack.push(currConf);
    visit(tree);
    return currConf;
  }

  // Visit mdfile - top level
  virtual antlrcpp::Any visitOmdfile(OMDParser::OmdfileContext *ctx) override {
    visitChildren(ctx);
    blockStack.top()->validate();
    blockStack.pop();
    return nullptr;
  }

  // Visit assignment
  virtual antlrcpp::Any visitAssignment(OMDParser::AssignmentContext *ctx) override {
    std::string id = ctx->ID()->getText();
    visitConstant(ctx->constant(), id);
    return nullptr;
  }

  // Helper method to handle constants
  void visitConstant(OMDParser::ConstantContext *ctx, const std::string& id) {
    if (ctx->intConst()) {
      int ival = getIntConst(ctx->intConst());
      blockStack.top()->assign(id, ival);
    }
    else if (ctx->floatConst()) {
      RealType dval = getFloatConst(ctx->floatConst());
      blockStack.top()->assign(id, dval);
    }
    else if (ctx->ID()) {
      std::string str = ctx->ID()->getText();
      blockStack.top()->assign(id, str);
    }
    else if (ctx->StringLiteral()) {
      std::string s = ctx->StringLiteral()->getText();
      s = s.substr(1, s.length()-2);
      blockStack.top()->assign(id, s);
    }
    else if (ctx->vectorConst()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(
			      	ctx->vectorConst()->doubleNumberTuple());
      blockStack.top()->assign(id, dvec);
    }
  }

  // Component block
  virtual antlrcpp::Any visitComponentblock(OMDParser::ComponentblockContext *ctx) override {
    Component* currComponent = new Component();
    blockStack.push(currComponent);
        
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    currConf->addComponent(currComponent);
    return nullptr;
  }

  // ZConstraint block
  virtual antlrcpp::Any visitZconstraintblock(OMDParser::ZconstraintblockContext *ctx) override {
    ZConsStamp* currZConsStamp = new ZConsStamp();
    blockStack.push(currZConsStamp);
        
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    currConf->addZConsStamp(currZConsStamp);
    return nullptr;
  }

  // Restraint block
  virtual antlrcpp::Any visitRestraintblock(OMDParser::RestraintblockContext *ctx) override {
    RestraintStamp* currRestraintStamp = new RestraintStamp();
    blockStack.push(currRestraintStamp);
        
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    currConf->addRestraintStamp(currRestraintStamp);
    return nullptr;
  }

  // FlucQ block
  virtual antlrcpp::Any visitFlucqblock(OMDParser::FlucqblockContext *ctx) override {
    FluctuatingChargeParameters* flucQpars = new FluctuatingChargeParameters();
    blockStack.push(flucQpars);
        
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    currConf->addFluctuatingChargeParameters(flucQpars);
    return nullptr;
  }

  // RNEMD block
  virtual antlrcpp::Any visitRnemdblock(OMDParser::RnemdblockContext *ctx) override {
    RNEMD::RNEMDParameters* rnemdPars = new RNEMD::RNEMDParameters();
    blockStack.push(rnemdPars);
        
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    currConf->addRNEMDParameters(rnemdPars);
    return nullptr;
  }

  // Light block
  virtual antlrcpp::Any visitLightblock(OMDParser::LightblockContext *ctx) override {
    Perturbations::LightParameters* lightPars = new Perturbations::LightParameters();
    blockStack.push(lightPars);
        
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    currConf->addLightParameters(lightPars);
    return nullptr;
  }

  // Minimizer block
  virtual antlrcpp::Any visitMinimizerblock(OMDParser::MinimizerblockContext *ctx) override {
    MinimizerParameters* minimizerPars = new MinimizerParameters();
    blockStack.push(minimizerPars);
        
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    currConf->addMinimizerParameters(minimizerPars);
    return nullptr;
  }

  // Molecule block
  virtual antlrcpp::Any visitMoleculeblock(OMDParser::MoleculeblockContext *ctx) override {
    MoleculeStamp* currMoleculeStamp = new MoleculeStamp();
    blockStack.push(currMoleculeStamp);
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    currConf->addMoleculeStamp(currMoleculeStamp);
    return nullptr;
  }

  // Sequence string
  virtual antlrcpp::Any visitSequencestring(OMDParser::SequencestringContext *ctx) override {
    std::string id = "sequence";
    visitConstant(ctx->constant(), id);
    return nullptr;
  }

  // Fragment block
  virtual antlrcpp::Any visitFragmentblock(OMDParser::FragmentblockContext *ctx) override {
    FragmentStamp* currFragmentStamp = new FragmentStamp();
    blockStack.push(currFragmentStamp);
        
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    currConf->addFragmentStamp(currFragmentStamp);
    return nullptr;
  }

  // Atom block
  virtual antlrcpp::Any visitAtomblock(OMDParser::AtomblockContext *ctx) override {
    int index = getIntConst(ctx->intConst());
    AtomStamp* currAtomStamp = new AtomStamp(index);
    blockStack.push(currAtomStamp);
        
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
    currMoleculeStamp->addAtomStamp(currAtomStamp);
    return nullptr;
  }

  // Atom statement
  virtual antlrcpp::Any visitAtomstatement(OMDParser::AtomstatementContext *ctx) override {
    AtomStamp* currAtomStamp = static_cast<AtomStamp*>(blockStack.top());
        
    if (ctx->assignment()) {
      visitAssignment(ctx->assignment());
    }
    else if (ctx->POSITION()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currAtomStamp->setPosition(dvec);
    }
    else if (ctx->ORIENTATION()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currAtomStamp->setOrientation(dvec);
    }
    else if (ctx->CHARGE()) {
      RealType rval = getFloatConst(ctx->floatConst());
      currAtomStamp->overrideCharge(rval);
    }
    return nullptr;
  }

  // Bond block
  virtual antlrcpp::Any visitBondblock(OMDParser::BondblockContext *ctx) override {
    BondStamp* currBondStamp = new BondStamp();
    blockStack.push(currBondStamp);
        
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
    currMoleculeStamp->addBondStamp(currBondStamp);
    return nullptr;
  }

  // Bond statement
  virtual antlrcpp::Any visitBondstatement(OMDParser::BondstatementContext *ctx) override {
    BondStamp* currBondStamp = static_cast<BondStamp*>(blockStack.top());
        
    if (ctx->assignment()) {
      visitAssignment(ctx->assignment());
    }
    else if (ctx->MEMBERS()) {
      std::vector<int> ivec = getIntTuple(ctx->inttuple());
      currBondStamp->setMembers(ivec);
    }
    else if (ctx->FIXED()) {
      RealType rval = getFloatConst(ctx->floatConst());
      currBondStamp->overrideType("Fixed", rval);
    }
    else if (ctx->HARMONIC()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currBondStamp->overrideType("Harmonic", dvec);
    }
    else if (ctx->CUBIC()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currBondStamp->overrideType("Cubic", dvec);
    }
    else if (ctx->QUARTIC()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currBondStamp->overrideType("Quartic", dvec);
    }
    else if (ctx->POLYNOMIAL()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currBondStamp->overrideType("Polynomial", dvec);
    }
    else if (ctx->MORSE()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currBondStamp->overrideType("Morse", dvec);
    }
    return nullptr;
  }

  // Bend block
  virtual antlrcpp::Any visitBendblock(OMDParser::BendblockContext *ctx) override {
    BendStamp* currBendStamp = new BendStamp();
    blockStack.push(currBendStamp);
        
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
    currMoleculeStamp->addBendStamp(currBendStamp);
    return nullptr;
  }

  // Bend statement
  virtual antlrcpp::Any visitBendstatement(OMDParser::BendstatementContext *ctx) override {
    BendStamp* currBendStamp = static_cast<BendStamp*>(blockStack.top());
        
    if (ctx->assignment()) {
      visitAssignment(ctx->assignment());
    }
    else if (ctx->MEMBERS()) {
      std::vector<int> ivec = getIntTuple(ctx->inttuple());
      currBendStamp->setMembers(ivec);
    }
    else if (ctx->HARMONIC()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currBendStamp->overrideType("Harmonic", dvec);
    }
    else if (ctx->GHOSTBEND()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currBendStamp->overrideType("GhostBend", dvec);
    }
    else if (ctx->UREYBRADLEY()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currBendStamp->overrideType("UreyBradley", dvec);
    }
    else if (ctx->CUBIC()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currBendStamp->overrideType("Cubic", dvec);
    }
    else if (ctx->QUARTIC()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currBendStamp->overrideType("Quartic", dvec);
    }
    else if (ctx->POLYNOMIAL()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currBendStamp->overrideType("Polynomial", dvec);
    }
    else if (ctx->COSINE()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currBendStamp->overrideType("Cosine", dvec);
    }
    return nullptr;
  }

  // Torsion block
  virtual antlrcpp::Any visitTorsionblock(OMDParser::TorsionblockContext *ctx) override {
    TorsionStamp* currTorsionStamp = new TorsionStamp();
    blockStack.push(currTorsionStamp);
        
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
    currMoleculeStamp->addTorsionStamp(currTorsionStamp);
    return nullptr;
  }

  // Torsion statement
  virtual antlrcpp::Any visitTorsionstatement(OMDParser::TorsionstatementContext *ctx) override {
    TorsionStamp* currTorsionStamp = static_cast<TorsionStamp*>(blockStack.top());
        
    if (ctx->assignment()) {
      visitAssignment(ctx->assignment());
    }
    else if (ctx->MEMBERS()) {
      std::vector<int> ivec = getIntTuple(ctx->inttuple());
      currTorsionStamp->setMembers(ivec);
    }
    else if (ctx->GHOSTTORSION()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currTorsionStamp->overrideType("GhostTorsion", dvec);
    }
    else if (ctx->CUBIC()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currTorsionStamp->overrideType("Cubic", dvec);
    }
    else if (ctx->QUARTIC()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currTorsionStamp->overrideType("Quartic", dvec);
    }
    else if (ctx->POLYNOMIAL()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currTorsionStamp->overrideType("Polynomial", dvec);
    }
    else if (ctx->CHARMM()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currTorsionStamp->overrideType("Charmm", dvec);
    }
    else if (ctx->OPLS()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currTorsionStamp->overrideType("Opls", dvec);
    }
    else if (ctx->TRAPPE()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currTorsionStamp->overrideType("Trappe", dvec);
    }
    else if (ctx->HARMONIC()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currTorsionStamp->overrideType("Harmonic", dvec);
    }
    return nullptr;
  }

  // Inversion block
  virtual antlrcpp::Any visitInversionblock(OMDParser::InversionblockContext *ctx) override {
    InversionStamp* currInversionStamp = new InversionStamp();
    blockStack.push(currInversionStamp);
        
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
    currMoleculeStamp->addInversionStamp(currInversionStamp);
    return nullptr;
  }

  // Inversion statement
  virtual antlrcpp::Any visitInversionstatement(OMDParser::InversionstatementContext *ctx) override {
    InversionStamp* currInversionStamp = static_cast<InversionStamp*>(blockStack.top());
        
    if (ctx->assignment()) {
      visitAssignment(ctx->assignment());
    }
    else if (ctx->CENTER()) {
      int icent = getIntConst(ctx->intConst());
      currInversionStamp->setCenter(icent);
    }
    else if (ctx->SATELLITES()) {
      std::vector<int> ivec = getIntTuple(ctx->inttuple());
      currInversionStamp->setSatellites(ivec);
    }
    else if (ctx->AMBERIMPROPER()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currInversionStamp->overrideType("AmberImproper", dvec);
    }
    else if (ctx->IMPROPERCOSINE()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currInversionStamp->overrideType("ImproperCosine", dvec);
    }
    else if (ctx->HARMONIC()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currInversionStamp->overrideType("Harmonic", dvec);
    }
    else if (ctx->CENTRALATOMHEIGHT()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currInversionStamp->overrideType("CentralAtomHeight", dvec);
    }
    else if (ctx->DREIDING()) {
      std::vector<RealType> dvec = getDoubleNumberTuple(ctx->doubleNumberTuple());
      currInversionStamp->overrideType("Dreiding", dvec);
    }
    return nullptr;
  }

  // Rigid body block
  virtual antlrcpp::Any visitRigidbodyblock(OMDParser::RigidbodyblockContext *ctx) override {
    int index = getIntConst(ctx->intConst());
    RigidBodyStamp* currRigidBodyStamp = new RigidBodyStamp(index);
    blockStack.push(currRigidBodyStamp);
        
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
    currMoleculeStamp->addRigidBodyStamp(currRigidBodyStamp);
    return nullptr;
  }

  // Rigid body statement
  virtual antlrcpp::Any visitRigidbodystatement(OMDParser::RigidbodystatementContext *ctx) override {
    RigidBodyStamp* currRigidBodyStamp = static_cast<RigidBodyStamp*>(blockStack.top());
        
    if (ctx->assignment()) {
      visitAssignment(ctx->assignment());
    }
    else if (ctx->MEMBERS()) {
      std::vector<int> ivec = getIntTuple(ctx->inttuple());
      currRigidBodyStamp->setMembers(ivec);
    }
    return nullptr;
  }

  // Cutoff group block
  virtual antlrcpp::Any visitCutoffgroupblock(OMDParser::CutoffgroupblockContext *ctx) override {
    CutoffGroupStamp* currCutoffGroupStamp = new CutoffGroupStamp();
    blockStack.push(currCutoffGroupStamp);
        
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
    currMoleculeStamp->addCutoffGroupStamp(currCutoffGroupStamp);
    return nullptr;
  }

  // Cutoff group statement
  virtual antlrcpp::Any visitCutoffgroupstatement(OMDParser::CutoffgroupstatementContext *ctx) override {
    CutoffGroupStamp* currCutoffGroupStamp = static_cast<CutoffGroupStamp*>(blockStack.top());
        
    if (ctx->assignment()) {
      visitAssignment(ctx->assignment());
    }
    else if (ctx->MEMBERS()) {
      std::vector<int> ivec = getIntTuple(ctx->inttuple());
      currCutoffGroupStamp->setMembers(ivec);
    }
    return nullptr;
  }

  // Constraint block
  virtual antlrcpp::Any visitConstraintblock(OMDParser::ConstraintblockContext *ctx) override {
    ConstraintStamp* currConstraintStamp = new ConstraintStamp();
    blockStack.push(currConstraintStamp);
        
    visitChildren(ctx);
        
    blockStack.pop();
    MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
    currMoleculeStamp->addConstraintStamp(currConstraintStamp);
    return nullptr;
  }

  // Constraint statement
  virtual antlrcpp::Any visitConstraintstatement(OMDParser::ConstraintstatementContext *ctx) override {
    ConstraintStamp* currConstraintStamp = static_cast<ConstraintStamp*>(blockStack.top());
        
    if (ctx->assignment()) {
      visitAssignment(ctx->assignment());
    }
    else if (ctx->MEMBERS()) {
      std::vector<int> ivec = getIntTuple(ctx->inttuple());
      currConstraintStamp->setMembers(ivec);
    }
    return nullptr;
  }

  // Nodes block
  virtual antlrcpp::Any visitNodesblock(OMDParser::NodesblockContext *ctx) override {
    NodesStamp* currNodesStamp = new NodesStamp();
    blockStack.push(currNodesStamp);
        
    visitChildren(ctx);
        
    blockStack.top()->validate();
    blockStack.pop();
    FragmentStamp* currFragmentStamp = static_cast<FragmentStamp*>(blockStack.top());
    currFragmentStamp->addNodesStamp(currNodesStamp);
    return nullptr;
  }

  // Nodes statement
  virtual antlrcpp::Any visitNodesstatement(OMDParser::NodesstatementContext *ctx) override {
    NodesStamp* currNodesStamp = static_cast<NodesStamp*>(blockStack.top());
        
    if (ctx->assignment()) {
      visitAssignment(ctx->assignment());
    }
    else if (ctx->MEMBERS()) {
      std::vector<int> ivec = getIntTuple(ctx->inttuple());
      currNodesStamp->setMembers(ivec);
    }
    return nullptr;
  }

private:
  // Helper methods to extract values from parse tree
    
  int getIntConst(OMDParser::IntConstContext *ctx) {
    if (ctx->NUM_INT()) {
      return lexi_cast<int>(ctx->NUM_INT()->getText());
    } else if (ctx->NUM_LONG()) {
      return lexi_cast<int>(ctx->NUM_LONG()->getText());
    }
    return 0;
  }
    
  RealType getFloatConst(OMDParser::FloatConstContext *ctx) {
    if (ctx->NUM_FLOAT()) {
      return lexi_cast<RealType>(ctx->NUM_FLOAT()->getText());
    } else if (ctx->NUM_DOUBLE()) {
      return lexi_cast<RealType>(ctx->NUM_DOUBLE()->getText());
    }
    return 0.0;
  }
    
  RealType getDoubleNumber(OMDParser::DoubleNumberContext *ctx) {
    if (ctx->intConst()) {
      return static_cast<RealType>(getIntConst(ctx->intConst()));
    } else if (ctx->floatConst()) {
      return getFloatConst(ctx->floatConst());
    }
    return 0.0;
  }
    
  std::vector<RealType> getDoubleNumberTuple(OMDParser::DoubleNumberTupleContext *ctx) {
    std::vector<RealType> dvec;
    for (auto dn : ctx->doubleNumber()) {
      dvec.push_back(getDoubleNumber(dn));
    }
    return dvec;
  }
    
  std::vector<int> getIntTuple(OMDParser::InttupleContext *ctx) {
    std::vector<int> ivec;
    for (auto ic : ctx->intConst()) {
      ivec.push_back(getIntConst(ic));
    }
    return ivec;
  }
};

#endif // OMDTREEVISITOR_HPP
