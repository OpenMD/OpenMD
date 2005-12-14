header
{
#include <stack>
#include "io/Globals.hpp"
#include "utils/StringUtils.hpp"
using namespace std;
using namespace oopse;
}
options
  {
  language = "Cpp";
  }
                     
class MDTreeParser extends TreeParser;

options
{
        k = 3;
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
          ;


assignment  : #(ASSIGNEQUAL id:ID constant[#id]) //{blockStack.top()->assign(#ID->getText(),);}
            ;
            
constant [ANTLR_USE_NAMESPACE(antlr)RefAST id]    
            : signedIntOrFloat[#id]
            | str1:ID {blockStack.top()->assign(id->getText(), str1->getText());}
            | str2:StringLiteral { std::string s =  str2->getText();
                                   s = s.substr(1, s.length()-2);
                                   blockStack.top()->assign(id->getText(),s);
                                 }
            ;
            
signedIntOrFloat [ANTLR_USE_NAMESPACE(antlr)RefAST id]  
{
  int ival;
  double dval;
}
              : #(MINUS (icMinus:intConst {ival = lexi_cast<int>(icMinus->getText()); ival = -ival; blockStack.top()->assign(id->getText(), ival);}
                | fcMinus:floatConst) {dval = lexi_cast<double>(fcMinus->getText());dval = -dval;  blockStack.top()->assign(id->getText(), dval);}
                ) 
              | (ic:intConst {ival = lexi_cast<int>(ic->getText()); blockStack.top()->assign(id->getText(), ival);}
                | fc:floatConst {dval = lexi_cast<double>(fc->getText());  blockStack.top()->assign(id->getText(), dval);} 
                )               
              ;

componentblock  : #(COMPONENT  {Component* currComponet = new Component(); blockStack.push(currComponet);}
                      (assignment)* 
                       ENDBLOCK ) {blockStack.top()->validate();blockStack.pop(); currConf->addComponent(currComponet);}
                ;
    
zconstraintblock  : #(ZCONSTRAINT {ZConsStamp* currZConsStamp = new ZConsStamp(); blockStack.push(currZConsStamp);}
                        (assignment)* 
                         ENDBLOCK ) {blockStack.top()->validate();blockStack.pop(); currConf->addZConsStamp(currZConsStamp);}
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
vector<double> dvec;
AtomStamp* currAtomStamp =  static_cast<AtomStamp*>(blockStack.top());

}
              : assignment
              | #(POSITION dvec=signedNumberTuple) {currAtomStamp->setPosition(dvec);}
              | #(ORIENTATION dvec=signedNumberTuple) {currAtomStamp->setOrientation(dvec);}
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


              
signedNumberTuple   returns [vector<double> dvec]
{
  double dval;
}
              : (dval=signedNumber {dvec.push_back(dval);})+  
              ;
                          
inttuple  returns [vector<int> ivec]
{
  int ival;
}
              : (ival=intConst {ivec.push_back(ival);})+ 
              ;

protected
intConst returns [int ival]
        : oival:OCTALINT {ival = lexi_cast<int>(oival->getText());} 
        | dival:DECIMALINT {ival = lexi_cast<int>(dival->getText());}
        | hival:HEXADECIMALINT {ival = lexi_cast<int>(hival->getText());}
        ;

protected
signedNumber  returns [double dval]
              : 
                ic:intConst {dval = lexi_cast<double>(ic->getText());}
                | fc:floatConst {dval = lexi_cast<double>(fc->getText());} 
                               
              ;
              
protected
floatConst returns [double dval]
        : d1:FLOATONE {dval = lexi_cast<double>(d1->getText());}  
        | d2:FLOATTWO {dval = lexi_cast<double>(d2->getText());} 
        ;
        
