#define _LARGEFILE_SOURCE64
#define _FILE_OFFSET_BITS 64

#include <string.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <utility>

#ifdef IS_MPI
#include <mpi.h>
#include "brains/mpiSimulation.hpp"
#endif // is_mpi

#include "io/ReadWrite.hpp"
#include "utils/simError.h"

RestraintWriter::RestraintWriter( SimInfo* the_entry_plug ){

  entry_plug = the_entry_plug;

#ifdef IS_MPI
  //sort the local atoms by global index
  sortByGlobalIndex();
#endif // is_mpi

}

RestraintWriter::~RestraintWriter( ){}

#ifdef IS_MPI

/**
 * A hook function to load balancing
 */

void RestraintWriter::update(){
  sortByGlobalIndex();		
}
  
/**
 * Auxiliary sorting function
 */
 
bool indexSortingCriterion2(const pair<int, int>& p1, const pair<int, int>& p2){
  return p1.second < p2.second;
}

/**
 * Sorting the local index by global index
 */
 
void RestraintWriter::sortByGlobalIndex(){
  Molecule* mols = entry_plug->molecules;  
  indexArray.clear();
  
  for(int i = 0; i < entry_plug->n_mol;i++) 
    indexArray.push_back(make_pair(i, mols[i].getGlobalIndex()));
  
  sort(indexArray.begin(), indexArray.end(), indexSortingCriterion2);	
}

#endif

void RestraintWriter::writeZangle(double currentTime){

  ofstream angleOut;
  vector<ofstream*> fileStreams;

#ifdef IS_MPI
  if(worldRank == 0 ){
#endif    
    angleOut.open( entry_plug->zAngleName.c_str(), ios::out | ios::trunc );
    if( !angleOut ){
      sprintf( painCave.errMsg,
	       "Could not open \"%s\" for zAngle output.\n",
	       entry_plug->zAngleName.c_str() );
      painCave.isFatal = 1;
      simError();
    }
#ifdef IS_MPI
  }
#endif // is_mpi

  fileStreams.push_back(&angleOut); 
  writeFrame(fileStreams, currentTime);

#ifdef IS_MPI
  angleOut.close();
#endif
  
}

void RestraintWriter::writeZangle(double currentTime, const char *in_name){

  ofstream finalOut;
  vector<ofstream*> fileStreams;

#ifdef IS_MPI
  if(worldRank == 0 ){
#endif    
    finalOut.open( in_name, ios::out | ios::trunc );
    if( !finalOut ){
      sprintf( painCave.errMsg,
	       "Could not open \"%s\" for zAngle output.\n",
	       in_name );
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

void RestraintWriter::writeFrame( vector<ofstream*>& outFile, double currentTime ){

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
   * To make a long story short, the original version of RestraintWriter
   * worked in the most inefficient way possible.  Node 0 would 
   * poke each of the node for an individual atom's formatted data 
   * as node 0 worked its way down the global index. This was particularly 
   * inefficient since the method blocked all processors at every atom 
   * (and did it twice!).
   *
   * An intermediate version of RestraintWriter could be described from Node
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
   * This new and improved RestraintWriter works in an even more efficient 
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
  double atomData;
  int isDirectional;
  char* atomTypeString;
  char MPIatomTypeString[MINIBUFFERSIZE];
  int nObjects;
  int msgLen; // the length of message actually recieved at master nodes
#endif //is_mpi

  double angle;
  DirectionalAtom* dAtom;
  int nTotObjects;
  StuntDouble* sd;
  char* molName;
  vector<StuntDouble*> integrableObjects;
  vector<StuntDouble*>::iterator iter;
  nTotObjects = entry_plug->getTotIntegrableObjects();
#ifndef IS_MPI
  
  for(k = 0; k < outFile.size(); k++)
    *outFile[k] << currentTime << " : omega values at this time\n";

  for( i=0; i<nTotObjects; i++ ){
    
    integrableObjects = entry_plug->molecules[i].getIntegrableObjects();
    
    for( iter = integrableObjects.begin();iter !=  integrableObjects.end(); ++iter){
      sd = *iter;
      
      sprintf( tempBuffer,
	       "%14.10lf\n",
	       sd->getZangle());
      strcpy( writeLine, tempBuffer );

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

    for(k = 0; k < outFile.size(); k++)
      *outFile[k] << currentTime << " fs: omega values at this time\n";
    
    currentIndex = 0;
    
    for (i = 0 ; i < mpiSim->getNMolGlobal(); i++ ) {
      
      // Get the Node number which has this atom;
      
      which_node = MolToProcMap[i];
      
      if (which_node != 0) {
	
        if (potatoes[which_node] + 1 >= MAXTAG) {
          // The potato was going to exceed the maximum value, 
          // so wrap this processor potato back to 0:         
	  
          potatoes[which_node] = 0;          
          MPI_Send(&potatoes[which_node], 1, MPI_INT, which_node, 0, 
		   MPI_COMM_WORLD);
          
        }
	
        myPotato = potatoes[which_node];        
	
        //recieve the number of integrableObject in current molecule
        MPI_Recv(&nCurObj, 1, MPI_INT, which_node,
		 myPotato, MPI_COMM_WORLD, &istatus);
        myPotato++;
        
        for(int l = 0; l < nCurObj; l++){
	  
          if (potatoes[which_node] + 1 >= MAXTAG) {
            // The potato was going to exceed the maximum value, 
            // so wrap this processor potato back to 0:         
	    
            potatoes[which_node] = 0;          
            MPI_Send(&potatoes[which_node], 1, MPI_INT, which_node, 0, 
		     MPI_COMM_WORLD);
            
          }
	  	  
          MPI_Recv(&atomData, 1, MPI_DOUBLE, which_node, myPotato, MPI_COMM_WORLD, &istatus);
          myPotato++;
	  
          // If we've survived to here, format the line:
	  sprintf( writeLine,
		   "%14.10lf\n",
		   atomData);

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
	  atomData = sd->getZangle();
	  
	  // If we've survived to here, format the line:
	  
	  sprintf( writeLine,
		   "%14.10lf\n",
		   atomData);
	  
	  for(k = 0; k < outFile.size(); k++)
	    *outFile[k] << writeLine;

	}//end for(iter = integrableObject.begin())
	
	currentIndex++;
      }

    }//end for(i = 0; i < mpiSim->getNmol())

    for(k = 0; k < outFile.size(); k++)
      outFile[k]->flush();
    
    sprintf( checkPointMsg,
	     "Successfully printed a zAngle.\n");
  
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
	
	for( iter = integrableObjects.begin(); 
	     iter != integrableObjects.end(); iter++ ){
	  
	  if (myPotato + 1 >= MAXTAG) {
	    
	    // The potato was going to exceed the maximum value, 
	    // so wrap this processor potato back to 0 (and block until
	    // node 0 says we can go:
	    
	    MPI_Recv(&myPotato, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &istatus);
	  }
          
	  sd = *iter;
          
	  atomData = sd->getZangle();
	  
	  MPI_Send(&atomData, 1, MPI_DOUBLE, 0,
		   myPotato, MPI_COMM_WORLD);
	  
	  myPotato++;  
	}
	
	currentIndex++;    
      }
    }
    
    sprintf( checkPointMsg,
             "Successfully dropped a zAngle.\n");
    MPIcheckPoint();                
 }
#endif // is_mpi
}

