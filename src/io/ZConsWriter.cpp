#include <algorithm>
#include <iostream>
#include <vector>
//#include <pair>
#include "ZConsWriter.hpp"
#include "simError.h"

using namespace std;

ZConsWriter::ZConsWriter(const char* filename, vector<ZConsParaItem>* thePara)
{
  //use master - slave mode, only master node writes to disk
#ifdef IS_MPI
  if(worldRank == 0){
#endif

   output.open(filename);
   
   if(!output){
     sprintf( painCave.errMsg,
              "Could not open %s for z constrain output \n",
         filename);
     painCave.isFatal = 1;
     simError();
   }
   output << "#number of z constrain molecules" << endl;
   output << "#global Index of molecule\tzPos" << endl;
   output << "#every frame will contain below data" <<endl;
   output << "#time(fs)" << endl;
   output << "#number of fixed z-constrain molecules" << endl;
   output << "#global Index of molecule\tzconstrain force\tcurrentZPos" << endl;

   parameters = thePara;
   writeZPos();

#ifdef IS_MPI
  }
#endif  

}

ZConsWriter::~ZConsWriter()
{

#ifdef IS_MPI
  if(worldRank == 0 ){
#endif  
  output.close();  
#ifdef IS_MPI  
  }
#endif
}

/**
 *
 */
void ZConsWriter::writeFZ(double time, int num, int* index, double* fz, double* curZPos, double* zpos){

#ifndef IS_MPI
  output << time << endl;
  output << num << endl;
  
  for(int i = 0; i < num; i++)
    output << index[i] <<"\t" << fz[i] << "\t" << curZPos[i] << "\t" << zpos[i] <<endl;

#else
  int totalNum;
  MPI_Allreduce(&num, &totalNum, 1, MPI_INT,MPI_SUM, MPI_COMM_WORLD); 
  
  if(worldRank == 0){
    output << time << endl;
    output << totalNum << endl;
  }
  
  int whichNode;
  enum CommType { RequesPosAndForce, EndOfRequest} status;
  double pos;
  double force;
  double zconsPos;
  int localIndex;
  MPI_Status ierr;
  int tag = 0;
  
  if(worldRank == 0){
    
    int globalIndexOfCurMol;
    int *MolToProcMap;
    MolToProcMap = mpiSim->getMolToProcMap();
    
    for(int i = 0; i < (int)(parameters->size()); i++){
      
      globalIndexOfCurMol = (*parameters)[i].zconsIndex;
      whichNode = MolToProcMap[globalIndexOfCurMol];
      
      if(whichNode == 0){
        
       for(int j = 0; j < num; j++)
        if(index[j] == globalIndexOfCurMol){
          localIndex = j;
          break;
        }

      force = fz[localIndex];
      pos = curZPos[localIndex];
      
      }
      else{
        status = RequesPosAndForce;
        MPI_Send(&status, 1, MPI_INT, whichNode, tag, MPI_COMM_WORLD);
        MPI_Send(&globalIndexOfCurMol, 1, MPI_INT, whichNode, tag, MPI_COMM_WORLD);
        MPI_Recv(&force, 1, MPI_DOUBLE, whichNode, tag, MPI_COMM_WORLD, &ierr);
        MPI_Recv(&pos, 1, MPI_DOUBLE, whichNode, tag, MPI_COMM_WORLD, &ierr);
        MPI_Recv(&zconsPos, 1, MPI_DOUBLE, whichNode, tag, MPI_COMM_WORLD, &ierr);
      }

     output << globalIndexOfCurMol << "\t" << force << "\t" << pos << "\t"<<  zconsPos << endl;
              
    } //End of Request Loop
    
    //Send ending request message to slave nodes    
    status = EndOfRequest;
    for(int i =1; i < mpiSim->getNProcessors(); i++)
      MPI_Send(&status, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
     
  }
  else{
  
    int whichMol;
    bool done = false;

    while (!done){  
      
      MPI_Recv(&status, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &ierr);
    
      switch (status){
          
         case RequesPosAndForce : 
          
           MPI_Recv(&whichMol, 1, MPI_INT, 0, tag, MPI_COMM_WORLD,&ierr);
    
           for(int i = 0; i < num; i++)
           if(index[i] == whichMol){
             localIndex = i;
             break;
           }
    
           MPI_Send(&fz[localIndex], 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);    
           MPI_Send(&curZPos[localIndex], 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);    
           MPI_Send(&zpos[localIndex], 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);     
           break;
       
        case EndOfRequest :
         
         done = true;
         break;
      }
      
    }
          
  }

#endif

}

/*
 *
 */
void ZConsWriter::writeZPos(){

#ifdef IS_MPI
  if(worldRank == 0){
#endif
    
    output << parameters->size() << endl;     
    
    for(int i =0 ; i < (int)(parameters->size()); i++)
      output << (*parameters)[i].zconsIndex << "\t" <<  (*parameters)[i].zPos << endl;

#ifdef IS_MPI
  }
#endif
}
