#define _LARGEFILE_SOURCE64
#define _FILE_OFFSET_BITS 64

#include <string.h>
#include <iostream>
#include <fstream>

#include "io/ReadWrite.hpp"
#include "utils/simError.h"


StatWriter::StatWriter( SimInfo* the_entry_plug ){

  entry_plug = the_entry_plug;

#ifdef IS_MPI
  if(worldRank == 0 ){
#endif // is_mpi

    outName =  entry_plug->statusName;

    //std::cerr << "Opening " << outName << " for stat\n";

    outFile.open(outName.c_str(), ios::out | ios::trunc );
    
    if( !outFile ){
      
      sprintf( painCave.errMsg,
	       "Could not open \"%s\" for stat output.\n",
	       outName.c_str());
      painCave.isFatal = 1;
      simError();
    }

    //outFile.setf( ios::scientific );
    outFile << "#time(fs)\tE_tot\tV\tKE\tT(K)\tP(atm)\tVol(A^3)\tH_conserved";

    if (entry_plug->useSolidThermInt || entry_plug->useLiquidThermInt)
      outFile << "\tV_raw";
    
    if (entry_plug->useSolidThermInt) 
      outFile << "\tV_harm";

    outFile << "\n";

    
#ifdef IS_MPI
  }

  sprintf( checkPointMsg,
	   "Sucessfully opened output file for stating.\n");
  MPIcheckPoint();
#endif // is_mpi

  tStats = new Thermo( entry_plug );
}

StatWriter::~StatWriter( ){

#ifdef IS_MPI
  if(worldRank == 0 ){
#endif // is_mpi

    outFile.close();
    delete tStats;

#ifdef IS_MPI
  }
#endif // is_mpi
}

void StatWriter::writeStat( double currentTime ){

  double totE, potE, kinE, temp, press, vol;
  double conservedQuantity;

  totE = tStats->getTotalE();
  potE = tStats->getPotential();
  kinE = tStats->getKinetic();
  temp = tStats->getTemperature();
  press = tStats->getPressure();
  vol = tStats->getVolume();
  conservedQuantity = entry_plug->the_integrator->getConservedQuantity();

#ifdef IS_MPI
  if(worldRank == 0 ){
#endif // is_mpi

    outFile.precision(8);
    outFile
      << currentTime << "\t"
      << totE << "\t"
      << potE << "\t"
      << kinE << "\t"
      << temp << "\t"
      << press << "\t"
      << vol << "\t"
      << conservedQuantity;

    if (entry_plug->useSolidThermInt || entry_plug->useLiquidThermInt)
      outFile << "\t" << entry_plug->vRaw;
    
    if (entry_plug->useSolidThermInt) 
      outFile << "\t" << entry_plug->vHarm;

    outFile << "\n";

    outFile.flush();

#ifdef IS_MPI
  }
#endif // is_mpi
}

