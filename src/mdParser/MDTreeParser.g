header
{
#include <stack>
#include "io/Globals.hpp"
#include "utils/StringUtils.hpp"
using namespace std;
using namespace OpenMD;
}
options
  {
  language = "Cpp";
  }
                     
class MDTreeParser extends TreeParser;

options
{
        k = 1;
        importVocab = MD;
}
{
  public:
    Globals* walkTree(ANTLR_USE_NAMESPACE(antlr)RefAST tree)
    {
      currConf = new Globals;
      blockStack.push(currConf);
      mdfile(tree);
      return currConf;
    }
  private:
    Globals* currConf;
    stack<DataHolder*> blockStack;    
}
mdfile  : (statement)* {blockStack.top()->validate(); blockStack.pop();}
        ;

statement : assignment
    | componentblock
    | moleculeblock
    | zconstraintblock
    | restraintblock
    | flucqblock
    | rnemdblock
    | minimizerblock
    ;

assignment  : #(ASSIGNEQUAL id:ID constant[#id]) //{blockStack.top()->assign(#ID->getText(),);}
            ;
            
constant [ANTLR_USE_NAMESPACE(antlr)RefAST id]
{
    int ival;
    RealType dval, x, y, z;
    Vector3d dvec;
}    
    : ival=intConst {blockStack.top()->assign(id->getText(), ival);}
    | dval=floatConst {blockStack.top()->assign(id->getText(), dval);}
    | #(LPAREN x=doubleNumber COMMA y=doubleNumber COMMA z=doubleNumber RPAREN) 
        {   dvec.x() = x; 
            dvec.y() = y; 
            dvec.z() = z; 
            blockStack.top()->assign(id->getText(), dvec);
        }
    | str1:ID {blockStack.top()->assign(id->getText(), str1->getText());}
    | str2:StringLiteral {std::string s =  str2->getText();
            s = s.substr(1, s.length()-2);
            blockStack.top()->assign(id->getText(),s);
        }
    ;
            

componentblock  : #(COMPONENT  {Component* currComponet = new Component(); blockStack.push(currComponet);}
                      (assignment)* 
                       ENDBLOCK ) {blockStack.top()->validate();blockStack.pop(); currConf->addComponent(currComponet);}
                ;
    
zconstraintblock  : #(ZCONSTRAINT {ZConsStamp* currZConsStamp = new ZConsStamp(); blockStack.push(currZConsStamp);}
                        (assignment)* 
                         ENDBLOCK ) {blockStack.top()->validate();blockStack.pop(); currConf->addZConsStamp(currZConsStamp);}
                  ;

restraintblock  : #(RESTRAINT {RestraintStamp* currRestraintStamp = new RestraintStamp(); blockStack.push(currRestraintStamp);}
                        (assignment)* 
                         ENDBLOCK ) {blockStack.top()->validate();blockStack.pop(); currConf->addRestraintStamp(currRestraintStamp);}
                  ;
  
flucqblock  : #(FLUCQ  {FluctuatingChargeParameters* flucQpars = new FluctuatingChargeParameters(); blockStack.push(flucQpars);}
                      (assignment)* 
                       ENDBLOCK ) {blockStack.top()->validate();blockStack.pop(); currConf->addFluctuatingChargeParameters(flucQpars);}
                ;

rnemdblock  : #(RNEMD  {RNEMDParameters* rnemdPars = new RNEMDParameters(); blockStack.push(rnemdPars);}
                      (assignment)* 
                       ENDBLOCK ) {blockStack.top()->validate();blockStack.pop(); currConf->addRNEMDParameters(rnemdPars);}
                ;

minimizerblock  : #(MINIMIZER  {MinimizerParameters* minimizerPars = new MinimizerParameters(); blockStack.push(minimizerPars);}
                      (assignment)* 
                       ENDBLOCK ) {blockStack.top()->validate();blockStack.pop(); currConf->addMinimizerParameters(minimizerPars);}
                ;


moleculeblock : #(MOLECULE {MoleculeStamp* currMoleculeStamp = new MoleculeStamp(); blockStack.push(currMoleculeStamp);}
                    (moleculestatement)* 
                     ENDBLOCK ) {blockStack.top()->validate(); blockStack.pop(); currConf->addMoleculeStamp(currMoleculeStamp);}
              ;

moleculestatement : assignment
                  | atomblock
                  | bondblock
                  | bendblock
                  | torsionblock
                  | inversionblock
                  | rigidbodyblock
                  | cutoffgroupblock
                  | fragmentblock
                  ;

atomblock 
{
  int index;
}
          : #(ATOM index=intConst {AtomStamp* currAtomStamp = new AtomStamp(index); blockStack.push(currAtomStamp);} 
                  (atomstatement)* 
                  ENDBLOCK ) {
                                blockStack.top()->validate();
                                blockStack.pop(); 
                                MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
                                currMoleculeStamp->addAtomStamp(currAtomStamp); 
                             }
          ;

atomstatement 
{
vector<RealType> dvec;
AtomStamp* currAtomStamp =  static_cast<AtomStamp*>(blockStack.top());

}
              : assignment
              | #(POSITION dvec=doubleNumberTuple) {currAtomStamp->setPosition(dvec);}
              | #(ORIENTATION dvec=doubleNumberTuple) {currAtomStamp->setOrientation(dvec);}
              ;

                      
bondblock : #(BOND {BondStamp* currBondStamp = new BondStamp(); blockStack.push(currBondStamp);}
                (bondstatement)* 
                 ENDBLOCK )  {
                                blockStack.pop(); 
                                MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
                                currMoleculeStamp->addBondStamp(currBondStamp); 
                             }
          ;

bondstatement 
{
  vector<int> ivec;
  BondStamp* currBondStamp = static_cast<BondStamp*>(blockStack.top());
}
              : assignment
              | #(MEMBERS ivec=inttuple) {currBondStamp->setMembers(ivec);}
              ;

bendblock : #(BEND {BendStamp* currBendStamp = new BendStamp(); blockStack.push(currBendStamp);}
                  (bendstatement)* 
                   ENDBLOCK)  {
                                blockStack.top()->validate();
                                blockStack.pop(); 
                                MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
                                currMoleculeStamp->addBendStamp(currBendStamp); 
                             }
          ;

bendstatement 
{
  vector<int> ivec;
  BendStamp* currBendStamp = static_cast<BendStamp*>(blockStack.top());
}
              : assignment
              | #(MEMBERS ivec=inttuple) {currBendStamp->setMembers(ivec);}
              ;

torsionblock  : #(TORSION {TorsionStamp* currTorsionStamp = new TorsionStamp(); blockStack.push(currTorsionStamp);}
                   (torsionstatement)*
                    ENDBLOCK )  {
                                  blockStack.top()->validate();
                                  blockStack.pop(); 
                                  MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
                                  currMoleculeStamp->addTorsionStamp(currTorsionStamp); 
                                }
          ;

torsionstatement
{
  vector<int> ivec;
  TorsionStamp* currTorsionStamp = static_cast<TorsionStamp*>(blockStack.top());
}  
              : assignment
              | #(MEMBERS ivec=inttuple) {currTorsionStamp->setMembers(ivec);}
              ;

inversionblock  : #(INVERSION {InversionStamp* currInversionStamp = new InversionStamp(); blockStack.push(currInversionStamp);}
                   (inversionstatement)*
                    ENDBLOCK )  {
                                  blockStack.top()->validate();
                                  blockStack.pop(); 
                                  MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
                                  currMoleculeStamp->addInversionStamp(currInversionStamp); 
                                }
          ;

inversionstatement
{
  int icent;
  InversionStamp* currInversionStamp = static_cast<InversionStamp*>(blockStack.top());
}  
              : assignment
              | #(CENTER icent=intConst) {currInversionStamp->setCenter(icent);}
              ;

rigidbodyblock
{
int index;
}
              : #(RIGIDBODY index=intConst {RigidBodyStamp* currRigidBodyStamp = new RigidBodyStamp(index); blockStack.push(currRigidBodyStamp);}
                        (rigidbodystatement)* 
                        ENDBLOCK )  {
                                      blockStack.top()->validate();
                                      blockStack.pop(); 
                                      MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
                                      currMoleculeStamp->addRigidBodyStamp(currRigidBodyStamp); 
                                    }
                ;

rigidbodystatement  
{
  vector<int> ivec;
  RigidBodyStamp* currRigidBodyStamp = static_cast<RigidBodyStamp*>(blockStack.top());
}
              : assignment
              | #(MEMBERS ivec=inttuple) {currRigidBodyStamp->setMembers(ivec);}
              ;

cutoffgroupblock  : #(CUTOFFGROUP  {CutoffGroupStamp* currCutoffGroupStamp = new CutoffGroupStamp(); blockStack.push(currCutoffGroupStamp);}
                          (cutoffgroupstatement)* 
                           ENDBLOCK ) {
                                        blockStack.top()->validate();
                                        blockStack.pop(); 
                                        MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
                                        currMoleculeStamp->addCutoffGroupStamp(currCutoffGroupStamp); 
                                      }
                  ;

cutoffgroupstatement
{
  vector<int> ivec;
  CutoffGroupStamp* currCutoffGroupStamp = static_cast<CutoffGroupStamp*>(blockStack.top());
}  
              : assignment
              | #(MEMBERS ivec=inttuple) {currCutoffGroupStamp->setMembers(ivec);}             
              ;

fragmentblock {int ival;}
               : #(FRAGMENT ival=intConst {FragmentStamp* currFragmentStamp = new FragmentStamp(ival); blockStack.push(currFragmentStamp);}
                      (fragmentstatement)* 
                      ENDBLOCK) {
                                  blockStack.top()->validate();
                                  blockStack.pop(); 
                                  MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
                                  currMoleculeStamp->addFragmentStamp(currFragmentStamp); 
                                }
              ;

fragmentstatement : assignment
    ;


              
doubleNumberTuple returns [vector<RealType> dvec]
{
    RealType dval;
}   
    : (dval=doubleNumber {dvec.push_back(dval);})+  
    ;


inttuple returns [vector<int> ivec]
{
    int ival;
}
    : (ival=intConst {ivec.push_back(ival);})+ 
    ;

protected
intConst returns [int ival]
    : i1:NUM_INT {ival = lexi_cast<int>(i1->getText());} 
    | i2:NUM_LONG {ival = lexi_cast<int>(i2->getText());}
    ;

protected
doubleNumber returns [RealType dval]
    : ic:intConst {dval = lexi_cast<RealType>(ic->getText());}
    | fc:floatConst {dval = lexi_cast<RealType>(fc->getText());}         
    ;
              
protected
floatConst returns [RealType dval]
    : d1:NUM_FLOAT {dval = lexi_cast<RealType>(d1->getText());}  
    | d2:NUM_DOUBLE {dval = lexi_cast<RealType>(d2->getText());} 
    ;
  
