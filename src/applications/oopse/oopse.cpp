#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>

#ifdef IS_MPI
#include <mpi.h>
#endif

#ifdef PROFILE
#include "profiling/mdProfile.hpp"
#endif // PROFILE

#include "utils/simError.h"
#include "brains/SimSetup.hpp"
#include "brains/SimInfo.hpp"
#include "primitives/Atom.hpp"
#include "integrators/Integrator.hpp"
#include "brains/Thermo.hpp"
#include "io/ReadWrite.hpp"
#include "minimizers/OOPSEMinimizer.hpp"

char* program_name;
using namespace std;

int main(int argc,char* argv[]){
  
  char* in_name;
  SimSetup* startMe;
  SimInfo* entry_plug;
   
  // first things first, all of the initializations

#ifdef IS_MPI
  MPI_Init( &argc, &argv ); // the MPI communicators
#endif
   
  initSimError();           // the error handler
  srand48( 1337 );          // the random number generator.

#ifdef PROFILE
  initProfile();
#endif //profile
  
  // check command line arguments, and set the input file
  
  program_name = argv[0]; // save the program name in case we need it
  
#ifdef IS_MPI
  if( worldRank == 0 ){
#endif
    std::cerr << 
      "  +----------------------------------------------------------------------+\n" <<
      "  |    ____  ____  ____  _____ ______  The OpenSource, Object-oriented   |\n" <<
      "  |   / __ \\/ __ \\/ __ \\/ ___// ____/  Parallel Simulation Engine.       |\n" <<
      "  |  / / / / / / / /_/ /\\__ \\/ __/                                       |\n" <<
      "  | / /_/ / /_/ / ____/___/ / /___     Copyright 2004 by the             |\n" <<
      "  | \\____/\\____/_/    /____/_____/     University of Notre Dame.         |\n" <<
      "  |                                    http://www.oopse.org              |\n" <<
      "  |                                                                      |\n" <<
      "  | OOPSE is an OpenScience project.  All source code is available for   |\n" <<
      "  | any use subject to only one condition:                               |\n" <<
      "  |                                                                      |\n" <<
      "  | Any published work resulting from the use of this code must cite the |\n" <<
      "  | following paper:       M. A. Meineke, C. F. Vardeman II, T. Lin,     |\n" <<
      "  |                        C. J. Fennell, and J. D. Gezelter,            |\n" <<
      "  |                        J. Comp. Chem. XX, XXXX (2004).               |\n" << 
      "  +----------------------------------------------------------------------+\n" <<
      "\n";
    
    if( argc < 2 ){
      strcpy( painCave.errMsg, "Error, a meta-data file is needed to run.\n" );
      painCave.isFatal = 1;
      simError();
    }
#ifdef IS_MPI
  }
#endif
  
  in_name = argv[1];

#ifdef IS_MPI
  strcpy( checkPointMsg, "Successful number of arguments" );
  MPIcheckPoint();
#endif
    
  // create the simulation objects, and get the show on the road

  entry_plug = new SimInfo();
  startMe = new SimSetup();

  startMe->setSimInfo( entry_plug );

#ifdef PROFILE
  startProfile( pro1 );
#endif //profile

  startMe->parseFile( in_name );

#ifdef PROFILE
  endProfile( pro1 );
  
  startProfile( pro2 );
#endif //profile

  startMe->createSim();
  delete startMe;

#ifdef PROFILE
  endProfile( pro2 );
  
  startProfile( pro3 );
#endif //profile

  if (!entry_plug->has_minimizer)
    entry_plug->the_integrator->integrate();
  else
    entry_plug->the_minimizer->minimize();
#ifdef PROFILE
  endProfile( pro3 );
 
  writeProfiles();
#endif //profile

#ifdef IS_MPI
  strcpy( checkPointMsg, "Oh what a lovely Tea Party!" );
  MPIcheckPoint();
  
  MPI_Finalize();	 
#endif

  return 0 ;
}
