 /*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
 
#include <sys/time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef IS_MPI
#include <mpi.h>

#include "brains/mpiSimulation.hpp"
#endif //is_mpi

#include "config.h"
#include "utils/simError.h"
#include "profiling/mdProfile.hpp"

namespace mdProfileSpace {

  class ProfileString{
  public:
    char myName[MAX_PROFILE_NAMELENGTH];
  };

  ProfileString theNames[N_PROFILES];
  
  struct timeval startTime[N_PROFILES];
  struct timeval endTime[N_PROFILES];
    
  double accumTime[N_PROFILES];
  
#ifdef IS_MPI
  double globalTime[N_PROFILES];
#endif //is_mpi

  
}

extern "C"{
  
  void F90_FUNC(gettimes, GETTIMES)(double* forceTime, 
				    double* commTime);
}


using namespace mdProfileSpace;


void initProfile( void ){

  int i;

  for( i=0;i<N_PROFILES;i++ ){
    
    accumTime[i] = 0.0;

#ifdef IS_MPI
    globalTime[i] = 0.0;
#endif //is_mpi
  }

  strncpy( theNames[pro1].myName, "Integrator->integrateStep()", MAX_PROFILE_NAMELENGTH );
  strncpy( theNames[pro2].myName, "Integrator->writes and stats", MAX_PROFILE_NAMELENGTH );
  strncpy( theNames[pro3].myName, "Integrator->preMove", MAX_PROFILE_NAMELENGTH );
  strncpy( theNames[pro4].myName, "Integrator->moveA", MAX_PROFILE_NAMELENGTH );
  strncpy( theNames[pro5].myName, "Integrator->CalcForce", MAX_PROFILE_NAMELENGTH );
  strncpy( theNames[pro6].myName, "Integrator->moveB", MAX_PROFILE_NAMELENGTH );
  strncpy( theNames[pro7].myName, "shortRange force calc", MAX_PROFILE_NAMELENGTH );
  strncpy( theNames[pro8].myName, "fortran force calc", MAX_PROFILE_NAMELENGTH );
}


void startProfile( proNames theProfile ){
  struct timezone tz;

  gettimeofday( &startTime[theProfile], &tz );
}

void endProfile( proNames theProfile ){
  struct timezone tz;
  double startVal, endVal;

  gettimeofday( &endTime[theProfile], &tz );

  startVal = (double)startTime[theProfile].tv_sec 
    + (double)startTime[theProfile].tv_usec / 1000000.0;

  endVal = (double)endTime[theProfile].tv_sec 
    + (double)endTime[theProfile].tv_usec / 1000000.0;
  
  accumTime[theProfile] += endVal - startVal;
}


void writeProfiles( void ){
 
  int i;
  double totalTime;
  double percentTime[N_PROFILES];
  int days, hours, minutes, secs, msecs;
  double donkey;
  
  double forceTime, commTime;
  
#ifdef IS_MPI
  int j;

  MPI_Status istatus;    

  double nodeTime, nodeForceTime, nodeCommTime;
  double nodeAccum[N_PROFILES];
  double nodePercent[N_PROFILES];

  double globalTime, globalForceTime, globalCommTime;
  double globalAccum[N_PROFILES];
  double globalPercent[N_PROFILES];
#endif // is_mpi


#ifndef IS_MPI // single processor version 

  totalTime = 0.0;
  for(i=0;i<N_PROFILES;i++)
    totalTime += accumTime[i];

  for(i=0;i<N_PROFILES;i++)
    percentTime[i] = accumTime[i] / totalTime;

  fprintf(stdout,
	  "  Time Spent      Percent Time                        Name\n"
	  "--------------  ----------------   -----------------------------------------\n"
	  );

  for(i=0;i<N_PROFILES;i++){
    fprintf(stdout,
	    " %12G    %14G     %40s\n",
	    accumTime[i],
	    percentTime[i],
	    theNames[i].myName );
  }

  days = (int)floor( totalTime / 86400 );
  donkey = totalTime - 86400 * days;

  hours = (int)floor( donkey / 3600 );
  donkey -= hours * 3600;

  minutes = (int)floor( donkey / 60 );
  donkey -= minutes * 60;

  secs = (int)donkey;
  msecs = (int)( (donkey - secs) * 1000 );

  F90_FUNC(gettimes, GETTIMES)(&forceTime, &commTime);

  fprintf( stdout,
	   "----------------------------------------------------------------------------\n"
	   "  Total Time = %03d:%02d:%02d:%02d.%03d ( %G sec )\n"
	   "\n"
	   "  From Fortran: forceTime = %G secs; communicationTime = %G secs.\n",
	   days,
	   hours,
	   minutes,
	   secs,
	   msecs,
	   totalTime,
	   forceTime,
	   commTime);

#else // the parrallel version

  if( worldRank == 0 ){
    
    double *nodeTots = new double[mpiSim->getNProcessors()];
    double *nodePercentTots = new double[mpiSim->getNProcessors()];
    
    totalTime = 0.0;
    for(i=0;i<N_PROFILES;i++)
      totalTime += accumTime[i];
    
    for(i=0;i<N_PROFILES;i++)
      percentTime[i] = accumTime[i] / totalTime;
    
    fprintf(stdout,
	    "\n"
	    "----------------------------------------------------------------------------\n"
	    "  Output from Node %d:   \n"
	    "\n"
	    "  Time Spent      Percent Time                        Name\n"
	    "--------------  ----------------   -----------------------------------------\n",
	    worldRank);
    
    for(i=0;i<N_PROFILES;i++){
      fprintf(stdout,
	      " %12G    %14G     %40s\n",
	      accumTime[i],
	      percentTime[i],
	      theNames[i].myName );
    }
    
    days = (int)floor( totalTime / 86400 );
    donkey = totalTime - 86400 * days;
    
    hours = (int)floor( donkey / 3600 );
    donkey -= hours * 3600;
    
    minutes = (int)floor( donkey / 60 );
    donkey -= minutes * 60;
    
    secs = (int)donkey;
    msecs = (int)( (donkey - secs) * 1000 );
    
    F90_FUNC(gettimes, GETTIMES)(&forceTime, &commTime);

    fprintf( stdout,
	     "----------------------------------------------------------------------------\n"
	     "  Total Time = %03d:%02d:%02d:%02d.%03d ( %G sec )\n"
	     "\n"
	     "  From Fortran: forceTime = %G secs; communicationTime = %G secs.\n",
	     days,
	     hours,
	     minutes,
	     secs,
	     msecs,
	     totalTime,
	     forceTime,
	     commTime);
    
    // now the rest of the nodes
    
    nodeTots[0] = totalTime;
    
    globalTime = totalTime;
    globalForceTime = forceTime;
    globalCommTime = commTime;
    for(i=0;i<N_PROFILES;i++)
      globalAccum[i] = accumTime[i];
    
    
    for(j=1;j<mpiSim->getNProcessors();j++){
      
      nodeTime = 0.0;
     
      MPI_Recv(nodeAccum, N_PROFILES, MPI_DOUBLE, j,
	       1, MPI_COMM_WORLD, &istatus );

      MPI_Recv(&nodeForceTime, 1, MPI_DOUBLE, j,
	       1, MPI_COMM_WORLD, &istatus );
      MPI_Recv(&nodeCommTime, 1, MPI_DOUBLE, j,
	       1, MPI_COMM_WORLD, &istatus );

      for(i=0;i<N_PROFILES;i++){
	nodeTime += nodeAccum[i];
      }
      
      for(i=0;i<N_PROFILES;i++)
	nodePercent[i] = nodeAccum[i] / nodeTime;
      
      fprintf(stdout,
	      "\n"
	      "----------------------------------------------------------------------------\n"
	      "  Output from Node %d:   \n"
	      "\n"
	      "  Time Spent      Percent Time                        Name\n"
	      "--------------  ----------------   -----------------------------------------\n",
	      j);
      
      for(i=0;i<N_PROFILES;i++){
	fprintf(stdout,
		" %12G    %14G     %40s\n",
		nodeAccum[i],
		nodePercent[i],
		theNames[i].myName );
      }
      
      days = (int)floor( nodeTime / 86400 );
      donkey = nodeTime - 86400 * days;
      
      hours = (int)floor( donkey / 3600 );
      donkey -= hours * 3600;
      
      minutes = (int)floor( donkey / 60 );
      donkey -= minutes * 60;
      
      secs = (int)donkey;
      msecs = (int)( (donkey - secs) * 1000 );
      
      fprintf( stdout,
	       "----------------------------------------------------------------------------\n"
	       "  Total Time = %03d:%02d:%02d:%02d.%03d ( %G sec )\n"
	       "\n"
	       "  From Fortran: forceTime = %G secs; communicationTime = %G secs.\n",
	       days,
	       hours,
	       minutes,
	       secs,
	       msecs,
	       nodeTime,
	       nodeForceTime,
	       nodeCommTime);
      
      for(i=0;i<N_PROFILES;i++)
	globalAccum[i] += nodeAccum[i];
      
      globalTime += nodeTime;
      globalForceTime += nodeForceTime;
      globalCommTime += nodeCommTime;
      nodeTots[j] = nodeTime;
    }

    // print out the totals
    
    for(j=0;j<mpiSim->getNProcessors();j++)
      nodePercentTots[j] = nodeTots[j] / globalTime;
    
    for(i=0;i<N_PROFILES;i++)
      globalPercent[i] = globalAccum[i] / globalTime;

    fprintf(stdout,
	    "\n"
	    "----------------------------------------------------------------------------\n"
	    "  Total Across Nodes\n"
	    "\n"
	    "  Time Spent      Percent Time                        Name\n"
	    "--------------  ----------------   -----------------------------------------\n",
	    j);
    
    for(i=0;i<N_PROFILES;i++){
      fprintf(stdout,
	      " %12G    %14G     %40s\n",
	      globalAccum[i],
	      globalPercent[i],
	      theNames[i].myName );
    }
    fprintf(stdout,
	    "\n"
	    "\n" );
    
    for(j=0;j<mpiSim->getNProcessors();j++){
      
      fprintf(stdout,
	      " %12G    %14G     node %d\n",
	      nodeTots[j],
	      nodePercentTots[j],
	      j );
    }
    
    days = (int)floor( globalTime / 86400 );
    donkey = nodeTime - 86400 * days;
    

    hours = (int)floor( donkey / 3600 );
    donkey -= hours * 3600;
    
    minutes = (int)floor( donkey / 60 );
    donkey -= minutes * 60;
    
    secs = (int)donkey;
    msecs = (int)( (donkey - secs) * 1000 );
    
    fprintf( stdout,
	     "----------------------------------------------------------------------------\n"
	     "  Total Time = %03d:%02d:%02d:%02d.%03d ( %G sec )\n"
	     "\n"
	     "  From Fortran: forceTime = %G secs; communicationTime = %G secs.\n",
	     days,
	     hours,
	     minutes,
	     secs,
	     msecs,
	     globalTime,
	     globalForceTime,
	     globalCommTime);
  }

  else{

    for(j=1;j<mpiSim->getNProcessors();j++){
      
      if( worldRank == j ){
	
	F90_FUNC(gettimes, GETTIMES)(&forceTime, &commTime);

	MPI_Send( accumTime, N_PROFILES, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD );
	MPI_Send( &forceTime, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD );
	MPI_Send( &commTime, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD );
      }
    }
  }
    
#endif // is_mpi
  

}
