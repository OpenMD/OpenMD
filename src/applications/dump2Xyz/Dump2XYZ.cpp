#include <iostream>
#include <fstream>
#include <string>

#include "brains/SimSetup.hpp"
#include "applications/dump2Xyz/Dump2XYZCmd.h"
#include "visitors/AtomVisitor.hpp"
#include "visitors/CompositeVisitor.hpp"
#include "visitors/RigidBodyVisitor.hpp"
#include "visitors/OtherVisitor.hpp"
#include "visitors/ZconsVisitor.hpp"

using namespace std;

int main(int argc, char* argv[]){
  gengetopt_args_info args_info;
  string dumpFileName;
  string mdFileName;
  char inFileName[2002];
  string xyzFileName;
  SimInfo* info;
  SimSetup startMe;
  DumpReader* dumpReader;
  ofstream xyzStream;
  int nframes;
  Molecule* mol;  
  vector<StuntDouble*> integrableObjects;
  vector<StuntDouble*>::iterator iter;
  vector<RigidBody*> myRigidBodies;
  vector<RigidBody*>::iterator rbIter;
  
  CompositeVisitor* compositeVisitor;
  SSDAtomVisitor* ssdVisitor;
  DefaultAtomVisitor* defaultAtomVisitor;
  LipidHeadVisitor* lipidVisitor;
  RBCOMVisitor* rbCOMVisitor;
  ReplicateVisitor* replicateVisitor;
  WrappingVisitor* wrappingVisitor;
  IgnoreVisitor* ignoreVisitor;
  XYZVisitor* xyzVisitor;
  ZConsVisitor* zconsVisitor;
  PrepareVisitor* prepareVisitor;
  WaterTypeVisitor* waterTypeVisitor;
  
  //parse the command line option
    if (cmdline_parser (argc, argv, &args_info) != 0)
      exit(1) ;


  //get the dumpfile name and meta-data file name
  if (args_info.input_given){
    dumpFileName = args_info.input_arg;
  }
  else{
    cerr << "Does not have input file name" << endl;
    exit(1);
  }
  mdFileName = dumpFileName;
  mdFileName = mdFileName.substr(0, mdFileName.rfind(".")) + ".md";

  if (args_info.output_given){
    xyzFileName = args_info.output_arg;
  }
  else{
    xyzFileName = dumpFileName;
    xyzFileName = xyzFileName.substr(0, xyzFileName.rfind(".")) + ".xyz";
  }

  //parse md file and set up the system
  info = new SimInfo();
  startMe.setSimInfo(info );

  strcpy(inFileName, mdFileName.c_str() );
  startMe.parseFile( inFileName );

  startMe.createSim();

  
  //creat visitor list
  compositeVisitor = new CompositeVisitor();

  //creat ignore visitor
  if(args_info.ignore_given ||args_info.water_flag){
    
    ignoreVisitor = new IgnoreVisitor();

    for(int i = 0; i < args_info.ignore_given; i++)
      ignoreVisitor->addIgnoreType(args_info.ignore_arg[i]);

    //ignore water 
    if(args_info.water_flag){
      ignoreVisitor->addIgnoreType("SSD");
      ignoreVisitor->addIgnoreType("SSD1");
      ignoreVisitor->addIgnoreType("SSD_E");
      ignoreVisitor->addIgnoreType("SSD_RF");
      ignoreVisitor->addIgnoreType("TIP3P_RB_0");
      ignoreVisitor->addIgnoreType("TIP4P_RB_0");
      ignoreVisitor->addIgnoreType("TIP5P_RB_0");
      ignoreVisitor->addIgnoreType("SPCE_RB_0");      
      ignoreVisitor->addIgnoreType("DPD_RB_0");
    }

     compositeVisitor->addVisitor(ignoreVisitor, 1000);
  }

  //creat RigidBody Visitor
  if(args_info.rigidbody_flag){
    rbCOMVisitor = new RBCOMVisitor(info);
    compositeVisitor->addVisitor(rbCOMVisitor, 900);
  }

  //compositeVisitor->addVisitor(lipidVisitor, 900);

  //creat SSD atom visitor
  ssdVisitor = new SSDAtomVisitor(info);
  compositeVisitor->addVisitor(ssdVisitor, 800);

  //creat default atom visitor
  defaultAtomVisitor = new DefaultAtomVisitor(info);
  compositeVisitor->addVisitor(defaultAtomVisitor, 700);

  //creat waterType visitor
  if(args_info.watertype_flag){
    waterTypeVisitor = new WaterTypeVisitor;
    compositeVisitor->addVisitor(waterTypeVisitor, 600);
  }

  //create ZconsVisitor
  if(args_info.zconstraint_flag){
    
    zconsVisitor = new ZConsVisitor(info);

    if(zconsVisitor->haveZconsMol())
      compositeVisitor->addVisitor(zconsVisitor, 500);
    else
      delete zconsVisitor;
  }

  //creat wrapping visitor

  if(args_info.periodicBox_flag){
    wrappingVisitor = new WrappingVisitor(info);
    compositeVisitor->addVisitor(wrappingVisitor, 400);
  }
  
  //creat replicate visitor
  if(args_info.repeatX_given > 0 || args_info.repeatY_given > 0 ||args_info.repeatY_given > 0){
    IntVec3 replicateOpt(args_info.repeatX_arg, args_info.repeatY_arg, args_info.repeatZ_arg);
    replicateVisitor = new ReplicateVisitor(info, replicateOpt);
    compositeVisitor->addVisitor(replicateVisitor, 300);
  }

  //creat xyzVisitor
  xyzVisitor = new XYZVisitor(info);
  compositeVisitor->addVisitor(xyzVisitor, 200);

  cout << compositeVisitor->toString();

  //creat prepareVisitor
  prepareVisitor = new PrepareVisitor();

  //open dump file
  dumpReader = new DumpReader(dumpFileName.c_str());
  nframes = dumpReader->getNframes();

  xyzStream .open(xyzFileName.c_str());
  
  for (int i = 0; i < nframes; i += args_info.frame_arg){
    dumpReader->readFrame(info, i);

    mol = info->molecules;

    //update atoms of rigidbody
    for(int j = 0; j < info->n_mol; j++){
      myRigidBodies = mol[j].getMyRigidBodies();

      for(rbIter = myRigidBodies.begin(); rbIter != myRigidBodies.end(); ++rbIter)
        (*rbIter)->updateAtoms();
    }    

    
    //prepare visit
    for(int j = 0; j < info->n_mol; j++){
      integrableObjects = mol[j].getIntegrableObjects();

      for(iter = integrableObjects.begin(); iter != integrableObjects.end(); ++iter)
        (*iter)->accept(prepareVisitor);
    }    

    //update visitor
    compositeVisitor->update();
    
    //visit stuntdouble
    for(int j = 0; j < info->n_mol; j++){
      integrableObjects = mol[j].getIntegrableObjects();

      for(iter = integrableObjects.begin(); iter != integrableObjects.end(); ++iter)
        (*iter)->accept(compositeVisitor);
    }

    xyzVisitor->writeFrame(xyzStream);
    xyzVisitor->clear();
    
  }//end for (int i = 0; i < nframes; i += args_info.frame_arg)

  xyzStream.close();


  delete compositeVisitor;
  delete info;


}
