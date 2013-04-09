/*
 * Copyright (c) 2013 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *                                                                      
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <algorithm>

#include "config.h"
#include "icosahedralBuilderCmd.h"
#include "utils/MoLocator.hpp"
#include "utils/Tuple.hpp"
#include "brains/Register.hpp"
#include "brains/SimInfo.hpp"
#include "brains/SimCreator.hpp"
#include "io/DumpWriter.hpp"
#include "math/Vector3.hpp"
#include "math/SquareMatrix3.hpp"
#include "utils/StringUtils.hpp"

using namespace OpenMD;
using namespace std;

//
// Create Mackay icosaheron structure.
//
// Heavily modified from a code created by: Yanting Wang    07/21/2003
//

vector<Vector3d> Points;
vector<std::pair<int, int> > Edges;
vector<tuple3<int, int, int> > Facets;
vector<Vector3d> Basis; // Basis vectors of the edges

//
// function np
//
// Calculate number of particles on the nth layer.
//

int np( int n ) {
  if( n<0 ) return -1;
  else if( n==0 ) return 1;
  else if( n==1 ) return 12;
  else if( n==2 ) return 42;
  else {
    int count = 0;    
    count += 12;   // edge particles
    count += (n-1)*30;   // side particles
    for( int i = 1; i <= n-2; i++ ) count += i*20;   // body particles
    return count;
  }
}

//
// function init
//
// Initialize some constants.
//

void init() {
  
  Basis.clear();
  Edges.clear();
  Facets.clear();
  Points.clear();
  
  //
  // Initialize Basis vectors.
  //
  const RealType HT = ( sqrt(5.0) + 1.0 ) / 4.0;   // half Tau
  
  Basis.push_back( Vector3d( HT, 0.0, 0.5 ));
  Basis.push_back( Vector3d( HT, 0.0, -0.5 ));
  Basis.push_back( Vector3d( 0.5, HT, 0.0 ));
  Basis.push_back( Vector3d( -0.5, HT, 0.0 ));
  Basis.push_back( Vector3d( 0.0, 0.5, HT ));
  Basis.push_back( Vector3d( 0.0, -0.5, HT ));
  Basis.push_back( Vector3d( 0.5, -HT, 0.0 ));
  Basis.push_back( Vector3d( 0.0, 0.5, -HT ));
  Basis.push_back( Vector3d( -HT, 0.0, 0.5 ));
  Basis.push_back( Vector3d( 0.0, -0.5, -HT ));
  Basis.push_back( Vector3d( -HT, 0.0, -0.5 ));
  Basis.push_back( Vector3d( -0.5, -HT, 0.0 ));
  
  //
  // Initialize 30 edges
  //
  
  Edges.push_back(std::make_pair(0, 1));
  Edges.push_back(std::make_pair(0, 2));
  Edges.push_back(std::make_pair(0, 4));
  Edges.push_back(std::make_pair(0, 5));
  Edges.push_back(std::make_pair(0, 6));
  
  Edges.push_back(std::make_pair(10, 3));
  Edges.push_back(std::make_pair(10, 7));
  Edges.push_back(std::make_pair(10, 8));
  Edges.push_back(std::make_pair(10, 9));
  Edges.push_back(std::make_pair(10, 11));
  
  Edges.push_back(std::make_pair(1, 2));
  Edges.push_back(std::make_pair(1, 6));
  Edges.push_back(std::make_pair(1, 7));
  Edges.push_back(std::make_pair(1, 9));
  
  Edges.push_back(std::make_pair(8, 3));
  Edges.push_back(std::make_pair(8, 4));
  Edges.push_back(std::make_pair(8, 5));
  Edges.push_back(std::make_pair(8, 11));
  
  Edges.push_back(std::make_pair(2, 3));
  Edges.push_back(std::make_pair(2, 4));
  Edges.push_back(std::make_pair(2, 7));
  
  Edges.push_back(std::make_pair(11, 5));
  Edges.push_back(std::make_pair(11, 6));
  Edges.push_back(std::make_pair(11, 9));
  
  Edges.push_back(std::make_pair(6, 5));
  Edges.push_back(std::make_pair(6, 9));
  
  Edges.push_back(std::make_pair(3, 4));
  Edges.push_back(std::make_pair(3, 7));
  
  Edges.push_back(std::make_pair(7, 9));
  
  Edges.push_back(std::make_pair(5, 4));
  
  //
  // Initialize 20 facets
  //
  
  Facets.push_back(make_tuple3(0, 1, 2));
  Facets.push_back(make_tuple3(0, 2, 4));
  Facets.push_back(make_tuple3(0, 4, 5));
  Facets.push_back(make_tuple3(0, 5, 6));
  Facets.push_back(make_tuple3(0, 1, 6));
  
  Facets.push_back(make_tuple3(10, 3, 7));
  Facets.push_back(make_tuple3(10, 3, 8));
  Facets.push_back(make_tuple3(10, 8, 11));
  Facets.push_back(make_tuple3(10, 9, 11));
  Facets.push_back(make_tuple3(10, 7, 9));
                                    
  Facets.push_back(make_tuple3(1, 2, 7));
  Facets.push_back(make_tuple3(1, 7, 9));
  Facets.push_back(make_tuple3(1, 6, 9));
  
  Facets.push_back(make_tuple3(8, 5, 11));
  Facets.push_back(make_tuple3(8, 4, 5));
  Facets.push_back(make_tuple3(8, 3, 4));
  
  Facets.push_back(make_tuple3(2, 3, 7));
  Facets.push_back(make_tuple3(2, 3, 4));
  
  Facets.push_back(make_tuple3(11, 5, 6));
  Facets.push_back(make_tuple3(11, 6, 9));
}

//
// function ih
//
// Create nth layer particles.
// The distance between nearest neighbors has the unit length of 1.

vector<Vector3d> ih( int n ) {

  if( n < 0 ) return Points; 
  
  if( n==0 ) {
    // center particle only
    Points.push_back(Vector3d( 0.0, 0.0, 0.0 ));
    return Points;
  } 
  
  //
  // Generate edge particles
  //
  for( vector<Vector3d>::iterator i = Basis.begin(); i != Basis.end(); ++i ) {
    
    Points.push_back( (*i) * RealType(n) );
  }
  
  //
  // Generate side particles
  //
  
  if( n<2 ) return Points;
  
  for( vector<pair<int,int> >::iterator i=Edges.begin(); 
       i != Edges.end(); ++i ) {
    
    Vector3d e1 = Basis[ (*i).first  ] * RealType(n);
    Vector3d e2 = Basis[ (*i).second ] * RealType(n);
    
    for( int j = 1; j <= n-1; j++ ) {
      Points.push_back( e1 + (e2-e1) * RealType(j) / RealType(n));
    }      
  }
  
  //
  // Generate body particles
  //
  
  if( n<3 ) return Points;
  
  for( vector<tuple3<int,int,int> >::iterator i = Facets.begin();
       i != Facets.end(); ++i) {
    
    Vector3d e1 = Basis[ (*i).first  ] * RealType(n);
    Vector3d e2 = Basis[ (*i).second ] * RealType(n);
    Vector3d e3 = Basis[ (*i).third  ] * RealType(n);
    
    for( int j=1; j<=n-2; j++ ) {
      
      Vector3d v1 = e1 + (e2-e1) * RealType(j+1) / RealType(n);
      Vector3d v2 = e1 + (e3-e1) * RealType(j+1) / RealType(n);
      
      for( int k=1; k<=j; k++ ) {
        Points.push_back(v1 + (v2-v1) * RealType(k) / RealType(j+1));
      }
    }
  }
  return Points;
}


void createMdFile(const std::string&oldMdFileName, 
                  const std::string&newMdFileName,
                  int nMol) {
  ifstream oldMdFile;
  ofstream newMdFile;
  const int MAXLEN = 65535;
  char buffer[MAXLEN];
  
  //create new .md file based on old .md file
  oldMdFile.open(oldMdFileName.c_str());
  newMdFile.open(newMdFileName.c_str());
  oldMdFile.getline(buffer, MAXLEN);

  unsigned int i = 0;
  while (!oldMdFile.eof()) {

    //correct molecule number
    if (strstr(buffer, "nMol") != NULL) {
      sprintf(buffer, "\tnMol = %i;", nMol);
      newMdFile << buffer << std::endl;
    } else {
      newMdFile << buffer << std::endl;
    }
    
    oldMdFile.getline(buffer, MAXLEN);
  }
  
  oldMdFile.close();
  newMdFile.close();
}

int main(int argc, char *argv []) {
  
  gengetopt_args_info args_info;
  std::string inputFileName;
  std::string outputFileName;
  
  MoLocator* locator;
  RealType latticeConstant;
  int nShells;
 
  DumpWriter *writer;

  init();
  
  // Parse Command Line Arguments
  if (cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);
         
  /* get input file name */
  if (args_info.inputs_num)
    inputFileName = args_info.inputs[0];
  else {
    sprintf(painCave.errMsg, "No input .md file name was specified "
            "on the command line");
    painCave.isFatal = 1;
    cmdline_parser_print_help();
    simError();
  }

  if (args_info.shells_given) {   
    nShells = args_info.shells_arg;    
    if( nShells < 0 ) {
      sprintf(painCave.errMsg, "icosahedralBuilder:  The number of shells\n"
              "\tmust be greater than or equal to zero.");
      painCave.isFatal = 1;
      cmdline_parser_print_help();
      simError();
    }
  } else {
    sprintf(painCave.errMsg, "icosahedralBuilder:  The number of shells\n"
            "\tis required to build a Mackay Icosahedron.");
    painCave.isFatal = 1;
    cmdline_parser_print_help();
    simError();
  }

  if (args_info.latticeConstant_given) {    
    latticeConstant = args_info.latticeConstant_arg;
  } else  {
    
    int count=0;
    for( int i = 0; i <= nShells; i++ ) count += np( i );
    sprintf(painCave.errMsg, "icosahedralBuilder:  No lattice constant\n"
            "\tgiven.  Total number of atoms in a Mackay Icosahedron with\n"
            "\t%d shells is %d.", nShells, count);
    painCave.isFatal = 1;
    cmdline_parser_print_help();
    simError();
  }

  /* parse md file and set up the system */
  SimCreator oldCreator;
  SimInfo* oldInfo = oldCreator.createSim(inputFileName, false);

  Globals* simParams = oldInfo->getSimParams();

  //generate the coordinates
  for( int i = 0; i <= nShells; i++ ) ih( i );

  outputFileName = args_info.output_arg;
   
  //creat new .md file on fly which corrects the number of molecule    

  createMdFile(inputFileName, outputFileName, Points.size());
  
  if (oldInfo != NULL)
    delete oldInfo;
  
  SimCreator newCreator;
  SimInfo* NewInfo = newCreator.createSim(outputFileName, false);
    
  // Place molecules
  Molecule* mol;
  SimInfo::MoleculeIterator mi;
  mol = NewInfo->beginMolecule(mi);

  int l = 0;

  locator = new MoLocator(NewInfo->getMoleculeStamp(0), 
                          NewInfo->getForceField());
    
  Vector3d boxMax;
  Vector3d boxMin;

  for (int n = 0; n < Points.size(); n++) {
    mol = NewInfo->getMoleculeByGlobalIndex(l);
    Vector3d location = Points[n] * latticeConstant;
    Vector3d orientation = Vector3d(0, 0, 1.0);

    if (n == 0) {
      boxMax = location;
      boxMin = location;
    } else {
      for (int i = 0; i < 3; i++) {
        boxMax[i] = max(boxMax[i], location[i]);
        boxMin[i] = min(boxMin[i], location[i]);
      }
    }
    
    locator->placeMol(location, orientation, mol);
    l++;
  }

  Mat3x3d boundingBox;
  boundingBox(0,0) = 10.0*(boxMax[0] - boxMin[0]);
  boundingBox(1,1) = 10.0*(boxMax[1] - boxMin[1]);
  boundingBox(2,2) = 10.0*(boxMax[2] - boxMin[2]);
 
  //set Hmat
  NewInfo->getSnapshotManager()->getCurrentSnapshot()->setHmat( boundingBox );
  
  //create dumpwriter and write out the coordinates
  writer = new DumpWriter(NewInfo, outputFileName);
  
  if (writer == NULL) {
    sprintf(painCave.errMsg, "Error in creating dumpwriter object ");
    painCave.isFatal = 1;
    simError();
  }
  
  writer->writeDump();

  // deleting the writer will put the closing at the end of the dump file

  delete writer;

  // cleanup a by calling sim error.....
  sprintf(painCave.errMsg, "A new OpenMD file called \"%s\" has been "
          "generated.\n", outputFileName.c_str());
  painCave.isFatal = 0;
  painCave.severity = OPENMD_INFO;
  simError();
  return 0;
}
