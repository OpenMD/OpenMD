#define _LARGEFILE_SOURCE64
#define _FILE_OFFSET_BITS 64

#include <string.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <utility>

#ifdef IS_MPI
#include <mpi.h>
#include "mpiSimulation.hpp"

namespace dWrite{
  void DieDieDie( void );
}

using namespace dWrite;
#endif //is_mpi

#include "ReadWrite.hpp"
#include "simError.h"

DumpWriter::DumpWriter( SimInfo* the_entry_plug ){

  entry_plug = the_entry_plug;

#ifdef IS_MPI
  if(worldRank == 0 ){
#endif // is_mpi

    dumpFile.open(entry_plug->sampleName.c_str(), ios::out | ios::trunc );

    if( !dumpFile ){

      sprintf( painCave.errMsg,
	       "Could not open \"%s\" for dump output.\n",
	       entry_plug->sampleName.c_str());
      painCave.isFatal = 1;
      simError();
    }

#ifdef IS_MPI
  }

  //sort the local atoms by global index
  sortByGlobalIndex();
  
  sprintf( checkPointMsg,
	   "Sucessfully opened output file for dumping.\n");
  MPIcheckPoint();
#endif // is_mpi
}

DumpWriter::~DumpWriter( ){

#ifdef IS_MPI
  if(worldRank == 0 ){
#endif // is_mpi

    dumpFile.close();

#ifdef IS_MPI
  }
#endif // is_mpi
}

#ifdef IS_MPI

/**
 * A hook function to load balancing
 */

void DumpWriter::update(){
  sortByGlobalIndex();		
}
  
/**
 * Auxiliary sorting function
 */
 
bool indexSortingCriterion(const pair<int, int>& p1, const pair<int, int>& p2){
  return p1.second < p2.second;
}

/**
 * Sorting the local index by global index
 */
 
void DumpWriter::sortByGlobalIndex(){
  Molecule* mols = entry_plug->molecules;  
  indexArray.clear();
  
  for(int i = 0; i < entry_plug->n_mol;i++) 
    indexArray.push_back(make_pair(i, mols[i].getGlobalIndex()));
  
  sort(indexArray.begin(), indexArray.end(), indexSortingCriterion);	
}

#endif

void DumpWriter::writeDump(double currentTime){

  ofstream finalOut;
  vector<ofstream*> fileStreams;

#ifdef IS_MPI
  if(worldRank == 0 ){
#endif    
    finalOut.open( entry_plug->finalName.c_str(), ios::out | ios::trunc );
    if( !finalOut ){
      sprintf( painCave.errMsg,
	       "Could not open \"%s\" for final dump output.\n",
	       entry_plug->finalName.c_str() );
      painCave.isFatal = 1;
      simError();
    }
#ifdef IS_MPI
  }
#endif // is_mpi

  fileStreams.push_back(&finalOut); 
  fileStreams.push_back(&dumpFile);

  writeFrame(fileStreams, currentTime);

#ifdef IS_MPI
  finalOut.close();
#endif
  	
}

void DumpWriter::writeFinal(double currentTime){

  ofstream finalOut;
  vector<ofstream*> fileStreams;

#ifdef IS_MPI
  if(worldRank == 0 ){
#endif // is_mpi

    finalOut.open( entry_plug->finalName.c_str(), ios::out | ios::trunc );

    if( !finalOut ){
      sprintf( painCave.errMsg,
	       "Could not open \"%s\" for final dump output.\n",
	       entry_plug->finalName.c_str() );
      painCave.isFatal = 1;
      simError();
    }

#ifdef IS_MPI
  }
#endif // is_mpi
  
  fileStreams.push_back(&finalOut);  
  writeFrame(fileStreams, currentTime);

#ifdef IS_MPI
  finalOut.close();
#endif
  
}

void DumpWriter::writeFrame( vector<ofstream*>& outFile, double currentTime ){

  const int BUFFERSIZE = 2000;
  const int MINIBUFFERSIZE = 100;

  char tempBuffer[BUFFERSIZE];  
  char writeLine[BUFFERSIZE];

  int i;
  unsigned int k;

#ifdef IS_MPI
  
  /*********************************************************************
   * Documentation?  You want DOCUMENTATION?
   * 
   * Why all the potatoes below?  
   *
   * To make a long story short, the original version of DumpWriter
   * worked in the most inefficient way possible.  Node 0 would 
   * poke each of the node for an individual atom's formatted data 
   * as node 0 worked its way down the global index. This was particularly 
   * inefficient since the method blocked all processors at every atom 
   * (and did it twice!).
   *
   * An intermediate version of DumpWriter could be described from Node
   * zero's perspective as follows:
   * 
   *  1) Have 100 of your friends stand in a circle.
   *  2) When you say go, have all of them start tossing potatoes at
   *     you (one at a time).
   *  3) Catch the potatoes.
   *
   * It was an improvement, but MPI has buffers and caches that could 
   * best be described in this analogy as "potato nets", so there's no 
   * need to block the processors atom-by-atom.
   * 
   * This new and improved DumpWriter works in an even more efficient 
   * way:
   * 
   *  1) Have 100 of your friend stand in a circle.
   *  2) When you say go, have them start tossing 5-pound bags of 
   *     potatoes at you.
   *  3) Once you've caught a friend's bag of potatoes,
   *     toss them a spud to let them know they can toss another bag.
   *
   * How's THAT for documentation?
   *
   *********************************************************************/

  int *potatoes;
  int myPotato;

  int nProc;
  int j, which_node, done, which_atom, local_index, currentIndex;
  double atomData[13];
  int isDirectional;
  char* atomTypeString;
  char MPIatomTypeString[MINIBUFFERSIZE];
  int nObjects;
  int msgLen; // the length of message actually recieved at master nodes
#endif //is_mpi

  double q[4], ji[3];
  DirectionalAtom* dAtom;
  double pos[3], vel[3];
  int nTotObjects;
  StuntDouble* sd;
  char* molName;
  vector<StuntDouble*> integrableObjects;
  vector<StuntDouble*>::iterator iter;
  nTotObjects = entry_plug->getTotIntegrableObjects();
#ifndef IS_MPI
  
  for(k = 0; k < outFile.size(); k++){
    *outFile[k] << nTotObjects << "\n";

    *outFile[k] << currentTime << ";\t"
               << entry_plug->Hmat[0][0] << "\t"
	             << entry_plug->Hmat[1][0] << "\t"
	             << entry_plug->Hmat[2][0] << ";\t"
               
               << entry_plug->Hmat[0][1] << "\t"
	             << entry_plug->Hmat[1][1] << "\t"
	             << entry_plug->Hmat[2][1] << ";\t"

	             << entry_plug->Hmat[0][2] << "\t"
	             << entry_plug->Hmat[1][2] << "\t"
	             << entry_plug->Hmat[2][2] << ";";

    //write out additional parameters, such as chi and eta
    *outFile[k] << entry_plug->the_integrator->getAdditionalParameters() << endl;
  }
  
  for( i=0; i< entry_plug->n_mol; i++ ){

    integrableObjects = entry_plug->molecules[i].getIntegrableObjects();
    molName = (entry_plug->compStamps[entry_plug->molecules[i].getStampID()])->getID();
    
    for( iter = integrableObjects.begin();iter !=  integrableObjects.end(); ++iter){
      sd = *iter;
      sd->getPos(pos);
      sd->getVel(vel);

      sprintf( tempBuffer,
  	     "%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",
  	     sd->getType(),
  	     pos[0],
  	     pos[1],
  	     pos[2],
  	     vel[0],
  	     vel[1],
  	     vel[2]);
      strcpy( writeLine, tempBuffer );

      if( sd->isDirectional() ){

        sd->getQ( q );
        sd->getJ( ji );

        sprintf( tempBuffer,
  	       "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
  	       q[0],
  	       q[1],
  	       q[2],
  	       q[3],
                 ji[0],
                 ji[1],
                 ji[2]);
        strcat( writeLine, tempBuffer );
      }
      else
        strcat( writeLine, "0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n" );
    
      for(k = 0; k < outFile.size(); k++)
        *outFile[k] << writeLine;      
    }

}

#else // is_mpi

  /* code to find maximum tag value */
  
  int *tagub, flag, MAXTAG;
  MPI_Attr_get(MPI_COMM_WORLD, MPI_TAG_UB, &tagub, &flag);
  if (flag) {
    MAXTAG = *tagub;
  } else {
    MAXTAG = 32767;
  }  

  int haveError;

  MPI_Status istatus;
  int nCurObj;
  int *MolToProcMap = mpiSim->getMolToProcMap();

  // write out header and node 0's coordinates

  if( worldRank == 0 ){

    // Node 0 needs a list of the magic potatoes for each processor;

    nProc = mpiSim->getNProcessors();
    potatoes = new int[nProc];

    //write out the comment lines
    for (i = 0; i < nProc; i++) 
      potatoes[i] = 0;
    
      for(k = 0; k < outFile.size(); k++){
        *outFile[k] << nTotObjects << "\n";

        *outFile[k] << currentTime << ";\t"
	                 << entry_plug->Hmat[0][0] << "\t"
	                 << entry_plug->Hmat[1][0] << "\t"
	                 << entry_plug->Hmat[2][0] << ";\t"

	                 << entry_plug->Hmat[0][1] << "\t"
	                 << entry_plug->Hmat[1][1] << "\t"
	                 << entry_plug->Hmat[2][1] << ";\t"

	                 << entry_plug->Hmat[0][2] << "\t"
	                 << entry_plug->Hmat[1][2] << "\t"
	                 << entry_plug->Hmat[2][2] << ";";
  
        *outFile[k] << entry_plug->the_integrator->getAdditionalParameters() << endl;
    }

    currentIndex = 0;

    for (i = 0 ; i < mpiSim->getNMolGlobal(); i++ ) {
      
      // Get the Node number which has this atom;
      
      which_node = MolToProcMap[i];
      
      if (which_node != 0) {
        
        if (potatoes[which_node] + 1 >= MAXTAG) {
          // The potato was going to exceed the maximum value, 
          // so wrap this processor potato back to 0:         

          potatoes[which_node] = 0;          
          MPI_Send(&potatoes[which_node], 1, MPI_INT, which_node, 0, MPI_COMM_WORLD);
          
        }

        myPotato = potatoes[which_node];        

        //recieve the number of integrableObject in current molecule
        MPI_Recv(&nCurObj, 1, MPI_INT, which_node,
		 myPotato, MPI_COMM_WORLD, &istatus);
        myPotato++;
        
        for(int l = 0; l < nCurObj; l++){

          if (potatoes[which_node] + 2 >= MAXTAG) {
            // The potato was going to exceed the maximum value, 
            // so wrap this processor potato back to 0:         

            potatoes[which_node] = 0;          
            MPI_Send(&potatoes[which_node], 1, MPI_INT, which_node, 0, MPI_COMM_WORLD);
            
          }

          MPI_Recv(MPIatomTypeString, MINIBUFFERSIZE, MPI_CHAR, which_node,
          myPotato, MPI_COMM_WORLD, &istatus);

          atomTypeString = MPIatomTypeString;

          myPotato++;

          MPI_Recv(atomData, 13, MPI_DOUBLE, which_node, myPotato, MPI_COMM_WORLD, &istatus);
          myPotato++;

          MPI_Get_count(&istatus, MPI_DOUBLE, &msgLen);

          if(msgLen  == 13) 
            isDirectional = 1;
          else
            isDirectional = 0;
          
          // If we've survived to here, format the line:
            
          if (!isDirectional) {
    	
            sprintf( writeLine,
    		 "%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",
    		 atomTypeString,
    		 atomData[0],
    		 atomData[1],
    		 atomData[2],
    		 atomData[3],
    		 atomData[4],
    		 atomData[5]);
    	
           strcat( writeLine, "0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n" );
    	
          } 
          else {
    	
          	sprintf( writeLine,
          		 "%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
          		 atomTypeString,
          		 atomData[0],
          		 atomData[1],
          		 atomData[2],
          		 atomData[3],
          		 atomData[4],
          		 atomData[5],
          		 atomData[6],
          		 atomData[7],
          		 atomData[8],
          		 atomData[9],
          		 atomData[10],
          		 atomData[11],
          		 atomData[12]);
            
          }
          
          for(k = 0; k < outFile.size(); k++)
            *outFile[k] << writeLine;            

        }// end for(int l =0)
        potatoes[which_node] = myPotato;

      }
      else {
        
       	haveError = 0;
        
  	    local_index = indexArray[currentIndex].first;        

        integrableObjects = (entry_plug->molecules[local_index]).getIntegrableObjects(); 

        for(iter= integrableObjects.begin(); iter != integrableObjects.end(); ++iter){    
  	        sd = *iter;
    	    atomTypeString = sd->getType();
    	    
    	    sd->getPos(pos);
    	    sd->getVel(vel);          
    	  
            atomData[0] = pos[0];
            atomData[1] = pos[1];
            atomData[2] = pos[2];

            atomData[3] = vel[0];
            atomData[4] = vel[1];
            atomData[5] = vel[2];
              
            isDirectional = 0;

            if( sd->isDirectional() ){

              isDirectional = 1;
                
              sd->getQ( q );
              sd->getJ( ji );

              for (int j = 0; j < 6 ; j++)
                atomData[j] = atomData[j];            
              
              atomData[6] = q[0];
              atomData[7] = q[1];
              atomData[8] = q[2];
              atomData[9] = q[3];
              
              atomData[10] = ji[0];
              atomData[11] = ji[1];
              atomData[12] = ji[2];
            }
            
            // If we've survived to here, format the line:
            
            if (!isDirectional) {
      	
              sprintf( writeLine,
      		 "%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",
      		 atomTypeString,
      		 atomData[0],
      		 atomData[1],
      		 atomData[2],
      		 atomData[3],
      		 atomData[4],
      		 atomData[5]);
      	
             strcat( writeLine, "0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n" );
      	
            } 
            else {
      	
            	sprintf( writeLine,
            		 "%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
            		 atomTypeString,
            		 atomData[0],
            		 atomData[1],
            		 atomData[2],
            		 atomData[3],
            		 atomData[4],
            		 atomData[5],
            		 atomData[6],
            		 atomData[7],
            		 atomData[8],
            		 atomData[9],
            		 atomData[10],
            		 atomData[11],
            		 atomData[12]);
              
            }
            
            for(k = 0; k < outFile.size(); k++)
              *outFile[k] << writeLine;
            
            
        }//end for(iter = integrableObject.begin())
        
      currentIndex++;
      }

    }//end for(i = 0; i < mpiSim->getNmol())
    
    for(k = 0; k < outFile.size(); k++)
      outFile[k]->flush();
    
    sprintf( checkPointMsg,
             "Sucessfully took a dump.\n");
    
    MPIcheckPoint();        
    
    delete[] potatoes;
    
  } else {

    // worldRank != 0, so I'm a remote node.  

    // Set my magic potato to 0:

    myPotato = 0;
    currentIndex = 0;
    
    for (i = 0 ; i < mpiSim->getNMolGlobal(); i++ ) {
      
      // Am I the node which has this integrableObject?
      
      if (MolToProcMap[i] == worldRank) {


        if (myPotato + 1 >= MAXTAG) {
	  
          // The potato was going to exceed the maximum value, 
          // so wrap this processor potato back to 0 (and block until
          // node 0 says we can go:
	  
          MPI_Recv(&myPotato, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &istatus);
          
        }

          local_index = indexArray[currentIndex].first;        
          integrableObjects = entry_plug->molecules[local_index].getIntegrableObjects();
          
          nCurObj = integrableObjects.size();
                      
          MPI_Send(&nCurObj, 1, MPI_INT, 0,
  		             myPotato, MPI_COMM_WORLD);
          myPotato++;

          for( iter = integrableObjects.begin(); iter  != integrableObjects.end(); iter++){

            if (myPotato + 2 >= MAXTAG) {
    	  
              // The potato was going to exceed the maximum value, 
              // so wrap this processor potato back to 0 (and block until
              // node 0 says we can go:
    	  
              MPI_Recv(&myPotato, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &istatus);
              
            }
            
            sd = *iter;
            
            atomTypeString = sd->getType();

            sd->getPos(pos);
            sd->getVel(vel);

            atomData[0] = pos[0];
            atomData[1] = pos[1];
            atomData[2] = pos[2];

            atomData[3] = vel[0];
            atomData[4] = vel[1];
            atomData[5] = vel[2];
              
            isDirectional = 0;

            if( sd->isDirectional() ){

                isDirectional = 1;
                
                sd->getQ( q );
                sd->getJ( ji );
                
                
                atomData[6] = q[0];
                atomData[7] = q[1];
                atomData[8] = q[2];
                atomData[9] = q[3];
      
                atomData[10] = ji[0];
                atomData[11] = ji[1];
                atomData[12] = ji[2];
              }

             
            strncpy(MPIatomTypeString, atomTypeString, MINIBUFFERSIZE);

            // null terminate the string before sending (just in case):
            MPIatomTypeString[MINIBUFFERSIZE-1] = '\0';

            MPI_Send(MPIatomTypeString, MINIBUFFERSIZE, MPI_CHAR, 0,
    		             myPotato, MPI_COMM_WORLD);
            
            myPotato++;
            
            if (isDirectional) {

              MPI_Send(atomData, 13, MPI_DOUBLE, 0,
                       myPotato, MPI_COMM_WORLD);
              
            } else {

              MPI_Send(atomData, 6, MPI_DOUBLE, 0,
                       myPotato, MPI_COMM_WORLD);
            }

            myPotato++;  

          }

          currentIndex++;    
          
        }
      
      }

    sprintf( checkPointMsg,
             "Sucessfully took a dump.\n");
    MPIcheckPoint();                
    
    }


  
#endif // is_mpi
}

#ifdef IS_MPI

// a couple of functions to let us escape the write loop

void dWrite::DieDieDie( void ){

  MPI_Finalize();
  exit (0);
}

#endif //is_mpi
