#include "Integrator.hpp"
#include "simError.h"
#include <math.h>

const double INFINITE_TIME = 10e30;
template<typename T> ZConstraint<T>::ZConstraint(SimInfo* theInfo,
                                                 ForceFields* the_ff): T(theInfo, the_ff),
                                                                       fzOut(NULL),
                                                                       curZconsTime(0),
                                                                       forcePolicy(NULL),
                                                                       usingSMD(false),
                                                                       hasZConsGap(false){
  //get properties from SimInfo
  GenericData* data;
  ZConsParaData* zConsParaData;
  DoubleData* sampleTime;
  DoubleData* tolerance;
  DoubleData* gap;
  DoubleData* fixtime;
  StringData* policy;
  StringData* filename; 
  IntData* smdFlag;
  double COM[3];

  //by default, the direction of constraint is z
  // 0 --> x
  // 1 --> y
  // 2 --> z
  whichDirection = 2;

  //estimate the force constant of harmonical potential
  double Kb = 1.986E-3 ; //in kcal/K

  double halfOfLargestBox = max(info->boxL[0], max(info->boxL[1], info->boxL[2])) /
                            2;
  zForceConst = Kb * info->target_temp / (halfOfLargestBox * halfOfLargestBox);

  //creat force Subtraction policy
  data = info->getProperty(ZCONSFORCEPOLICY_ID);
  if (!data){
    sprintf(painCave.errMsg,
            "ZConstraint Warning: User does not set force Subtraction policy, "
            "PolicyByMass is used\n");
    painCave.isFatal = 0;
    simError();      

    forcePolicy = (ForceSubtractionPolicy *) new PolicyByMass(this);
  }
  else{
    policy = dynamic_cast<StringData*>(data);

    if (!policy){
      sprintf(painCave.errMsg,
              "ZConstraint Error: Convertion from GenericData to StringData failure, "
              "PolicyByMass is used\n");
      painCave.isFatal = 0;
      simError();      

      forcePolicy = (ForceSubtractionPolicy *) new PolicyByMass(this);
    }
    else{
      if (policy->getData() == "BYNUMBER")
        forcePolicy = (ForceSubtractionPolicy *) new PolicyByNumber(this);
      else if (policy->getData() == "BYMASS")
        forcePolicy = (ForceSubtractionPolicy *) new PolicyByMass(this);
      else{
        sprintf(painCave.errMsg,
                "ZConstraint Warning: unknown force Subtraction policy, "
                "PolicyByMass is used\n");
        painCave.isFatal = 0;
        simError();      
        forcePolicy = (ForceSubtractionPolicy *) new PolicyByMass(this);
      }
    }
  }


  //retrieve sample time of z-contraint 
  data = info->getProperty(ZCONSTIME_ID);

  if (!data){
    sprintf(painCave.errMsg,
            "ZConstraint error: If you use an ZConstraint\n"
            " , you must set sample time.\n");
    painCave.isFatal = 1;
    simError();
  }
  else{
    sampleTime = dynamic_cast<DoubleData*>(data);

    if (!sampleTime){
      sprintf(painCave.errMsg,
              "ZConstraint error: Can not get property from SimInfo\n");
      painCave.isFatal = 1;
      simError();
    }
    else{
      this->zconsTime = sampleTime->getData();
    }
  }

  //retrieve output filename of z force
  data = info->getProperty(ZCONSFILENAME_ID);
  if (!data){
    sprintf(painCave.errMsg,
            "ZConstraint error: If you use an ZConstraint\n"
            " , you must set output filename of z-force.\n");
    painCave.isFatal = 1;
    simError();
  }
  else{
    filename = dynamic_cast<StringData*>(data);

    if (!filename){
      sprintf(painCave.errMsg,
              "ZConstraint error: Can not get property from SimInfo\n");
      painCave.isFatal = 1;
      simError();
    }
    else{
      this->zconsOutput = filename->getData();
    }
  }

  //retrieve tolerance for z-constraint molecuels
  data = info->getProperty(ZCONSTOL_ID);

  if (!data){
    sprintf(painCave.errMsg, "ZConstraint error: can not get tolerance \n");
    painCave.isFatal = 1;
    simError();
  }
  else{
    tolerance = dynamic_cast<DoubleData*>(data);

    if (!tolerance){
      sprintf(painCave.errMsg,
              "ZConstraint error: Can not get property from SimInfo\n");
      painCave.isFatal = 1;
      simError();
    }
    else{
      this->zconsTol = tolerance->getData();
    }
  }

  //quick hack here
  data = info->getProperty(ZCONSGAP_ID);

  if (data){
    gap = dynamic_cast<DoubleData*>(data);

    if (!gap){
      sprintf(painCave.errMsg,
              "ZConstraint error: Can not get property from SimInfo\n");
      painCave.isFatal = 1;
      simError();
    }
    else{
      this->hasZConsGap = true;
      this->zconsGap = gap->getData();
    }
  }



  data = info->getProperty(ZCONSFIXTIME_ID);

  if (data){
    fixtime = dynamic_cast<DoubleData*>(data);
    if (!fixtime){
      sprintf(painCave.errMsg,
              "ZConstraint error: Can not get zconsFixTime from SimInfo\n");
      painCave.isFatal = 1;
      simError();
    }
    else{
      this->zconsFixTime = fixtime->getData();
    }
  }
  else if(hasZConsGap){
      sprintf(painCave.errMsg,
              "ZConstraint error: must set fixtime if already set zconsGap\n");
      painCave.isFatal = 1;
      simError();
  }



  data = info->getProperty(ZCONSUSINGSMD_ID);

  if (data){
    smdFlag = dynamic_cast<IntData*>(data);

    if (!smdFlag){
      sprintf(painCave.errMsg,
              "ZConstraint error: Can not get property from SimInfo\n");
      painCave.isFatal = 1;
      simError();
    }
    else{
      this->usingSMD= smdFlag->getData() ? true : false;
    }

  }



  //retrieve index of z-constraint molecules
  data = info->getProperty(ZCONSPARADATA_ID);
  if (!data){
    sprintf(painCave.errMsg,
            "ZConstraint error: If you use an ZConstraint\n"
            " , you must set index of z-constraint molecules.\n");
    painCave.isFatal = 1;
    simError();
  }
  else{
    zConsParaData = dynamic_cast<ZConsParaData*>(data);

    if (!zConsParaData){
      sprintf(painCave.errMsg,
              "ZConstraint error: Can not get parameters of zconstraint method from SimInfo\n");
      painCave.isFatal = 1;
      simError();
    }
    else{
      parameters = zConsParaData->getData();

      //check the range of zconsIndex
      //and the minimum value of index is the first one (we already sorted the data)
      //the maximum value of index is the last one

      int maxIndex;
      int minIndex;
      int totalNumMol;

      minIndex = (*parameters)[0].zconsIndex;
      if (minIndex < 0){
        sprintf(painCave.errMsg, "ZConstraint error: index is out of range\n");
        painCave.isFatal = 1;
        simError();
      }

      maxIndex = (*parameters)[parameters->size() - 1].zconsIndex;

#ifndef IS_MPI
      totalNumMol = nMols;
#else
      totalNumMol = mpiSim->getNMolGlobal();   
#endif      

      if (maxIndex > totalNumMol - 1){
        sprintf(painCave.errMsg, "ZConstraint error: index is out of range\n");
        painCave.isFatal = 1;
        simError();
      }

      //if user does not specify the zpos for the zconstraint molecule
      //its initial z coordinate  will be used as default
      for (int i = 0; i < (int) (parameters->size()); i++){
        if (!(*parameters)[i].havingZPos){
#ifndef IS_MPI
          for (int j = 0; j < nMols; j++){
            if (molecules[j].getGlobalIndex() == (*parameters)[i].zconsIndex){
              molecules[j].getCOM(COM);
              break;
            }
          }
#else
          //query which processor current zconstraint molecule belongs to
          int* MolToProcMap;
          int whichNode;

          MolToProcMap = mpiSim->getMolToProcMap();
          whichNode = MolToProcMap[(*parameters)[i].zconsIndex];

          //broadcast the zpos of current z-contraint molecule
          //the node which contain this 

          if (worldRank == whichNode){
            for (int j = 0; j < nMols; j++)
              if (molecules[j].getGlobalIndex() == (*parameters)[i].zconsIndex){
                molecules[j].getCOM(COM);
                break;
              }
          }

          MPI_Bcast(&COM[whichDirection], 1, MPI_DOUBLE, whichNode,
                    MPI_COMM_WORLD);        
#endif

          (*parameters)[i].zPos = COM[whichDirection];

          sprintf(painCave.errMsg,
                  "ZConstraint warning: Does not specify zpos for z-constraint molecule "
                  "initial z coornidate will be used \n");
          painCave.isFatal = 0;
          simError();
        }
      }
    }//end if (!zConsParaData)

  }//end  if (!data)

  //  
#ifdef IS_MPI
  update(); 
#else  
  int searchResult;

  for (int i = 0; i < nMols; i++){
    searchResult = isZConstraintMol(&molecules[i]); 

    if (searchResult > -1){
      zconsMols.push_back(&molecules[i]);      
      massOfZConsMols.push_back(molecules[i].getTotalMass());  

      zPos.push_back((*parameters)[searchResult].zPos);
      kz.push_back((*parameters)[searchResult]. kRatio * zForceConst);
      
      if(usingSMD)
        cantVel.push_back((*parameters)[searchResult].cantVel);

    }
    else{
      unconsMols.push_back(&molecules[i]);
      massOfUnconsMols.push_back(molecules[i].getTotalMass());
    }
  }

  fz.resize(zconsMols.size());
  curZPos.resize(zconsMols.size());
  indexOfZConsMols.resize(zconsMols.size());  

  //determine the states of z-constraint molecules
  for (size_t i = 0; i < zconsMols.size(); i++){
    indexOfZConsMols[i] = zconsMols[i]->getGlobalIndex();

    zconsMols[i]->getCOM(COM);
    
    if (fabs(zPos[i] - COM[whichDirection]) < zconsTol){
      states.push_back(zcsFixed);

      if (hasZConsGap)
        endFixTime.push_back(info->getTime() + zconsFixTime);
    }
    else{
      states.push_back(zcsMoving);

      if (hasZConsGap)
        endFixTime.push_back(INFINITE_TIME);
    }

    if(usingSMD)
      cantPos.push_back(COM[whichDirection]);    
  }

  if(usingSMD)
    prevCantPos = cantPos;
#endif 

  
  //get total masss of unconstraint molecules
  double totalMassOfUncons_local;
  totalMassOfUncons_local = 0;

  for (size_t i = 0; i < unconsMols.size(); i++)
    totalMassOfUncons_local += unconsMols[i]->getTotalMass();

#ifndef IS_MPI
  totalMassOfUncons = totalMassOfUncons_local;
#else
  MPI_Allreduce(&totalMassOfUncons_local, &totalMassOfUncons, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);  
#endif

  //get total number of unconstrained atoms
  int nUnconsAtoms_local;
  nUnconsAtoms_local = 0;
  for (int i = 0; i < (int) (unconsMols.size()); i++)
    nUnconsAtoms_local += unconsMols[i]->getNAtoms();

#ifndef IS_MPI
  totNumOfUnconsAtoms = nUnconsAtoms_local;
#else
  MPI_Allreduce(&nUnconsAtoms_local, &totNumOfUnconsAtoms, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);  
#endif  

  forcePolicy->update();
}

template<typename T> ZConstraint<T>::~ZConstraint(){

  if (fzOut){
    delete fzOut;
  }

  if (forcePolicy){
    delete forcePolicy;
  }
}


/**
 * 
 */

#ifdef IS_MPI
template<typename T> void ZConstraint<T>::update(){
  double COM[3];
  int index;

  zconsMols.clear();
  massOfZConsMols.clear();
  zPos.clear();
  kz.clear();
  cantPos.clear();
  cantVel.clear();

  unconsMols.clear();
  massOfUnconsMols.clear();


  //creat zconsMol and unconsMol lists
  for (int i = 0; i < nMols; i++){
    index = isZConstraintMol(&molecules[i]); 

    if (index > -1){
      zconsMols.push_back(&molecules[i]);      
      zPos.push_back((*parameters)[index].zPos);
      kz.push_back((*parameters)[index].kRatio * zForceConst);
      massOfZConsMols.push_back(molecules[i].getTotalMass());
      
      if(usingSMD)
        cantVel.push_back((*parameters)[index].cantVel);

    }
    else{
      unconsMols.push_back(&molecules[i]);
      massOfUnconsMols.push_back(molecules[i].getTotalMass());
    }
  }

  fz.resize(zconsMols.size());
  curZPos.resize(zconsMols.size());
  indexOfZConsMols.resize(zconsMols.size());  
 
  for (size_t i = 0; i < zconsMols.size(); i++){
    indexOfZConsMols[i] = zconsMols[i]->getGlobalIndex();
  }
    
  //determine the states of z-constraint molecules
  for (int i = 0; i < (int) (zconsMols.size()); i++){

    zconsMols[i]->getCOM(COM);
    
    if (fabs(zPos[i] - COM[whichDirection]) < zconsTol){
      states.push_back(zcsFixed);

      if (hasZConsGap)
        endFixTime.push_back(info->getTime() + zconsFixTime);
    }
    else{
      states.push_back(zcsMoving);

      if (hasZConsGap)
        endFixTime.push_back(INFINITE_TIME);
    }

    if(usingSMD)
      cantPos.push_back(COM[whichDirection]);        
  } 

  if(usingSMD)
  prevCantPos = cantPos;

  forcePolicy->update();
}

#endif

/**
 *  Function Name: isZConstraintMol
 *  Parameter
 *    Molecule* mol
 *  Return value:
 *    -1, if the molecule is not z-constraint molecule, 
 *    other non-negative values, its index in indexOfAllZConsMols vector
 */

template<typename T> int ZConstraint<T>::isZConstraintMol(Molecule* mol){
  int index;
  int low;
  int high;
  int mid;

  index = mol->getGlobalIndex();

  low = 0;
  high = parameters->size() - 1;

  //Binary Search (we have sorted the array)  
  while (low <= high){
    mid = (low + high) / 2;
    if ((*parameters)[mid].zconsIndex == index)
      return mid;
    else if ((*parameters)[mid].zconsIndex > index)
      high = mid - 1;
    else
      low = mid + 1;
  }

  return -1;
}

template<typename T> void ZConstraint<T>::integrate(){
  // creat zconsWriter  
  fzOut = new ZConsWriter(zconsOutput.c_str(), parameters);   

  if (!fzOut){
    sprintf(painCave.errMsg, "Memory allocation failure in class Zconstraint\n");
    painCave.isFatal = 1;
    simError();
  }

  //zero out the velocities of center of mass of unconstrained molecules 
  //and the velocities of center of mass of every single z-constrained molecueles
  zeroOutVel();

  curZconsTime = zconsTime + info->getTime();

  T::integrate();
}

template<typename T> void ZConstraint<T>::calcForce(int calcPot, int calcStress){
  double zsys;
  double COM[3];
  double force[3];
  double zSysCOMVel;

  T::calcForce(calcPot, calcStress);


  if (hasZConsGap){
    updateZPos();
  }

  if (checkZConsState()){
    zeroOutVel();    
    forcePolicy->update();
  }  

  zsys = calcZSys();
  zSysCOMVel = calcSysCOMVel();
#ifdef IS_MPI
  if (worldRank == 0){
#endif

#ifdef IS_MPI
  }
#endif

  //do zconstraint force; 
  if (haveFixedZMols()){
    this->doZconstraintForce();
  }

  //use external force to move the molecules to the specified positions
  if (haveMovingZMols()){
    if (usingSMD)
      this->doHarmonic(cantPos);
    else
      this->doHarmonic(zPos);      
  }

  //write out forces and current positions of z-constraint molecules
  if (info->getTime() >= curZconsTime){
    for (int i = 0; i < (int) (zconsMols.size()); i++){
      zconsMols[i]->getCOM(COM);
      curZPos[i] = COM[whichDirection];

      //if the z-constraint molecule is still moving, just record its force
      if (states[i] == zcsMoving){
        fz[i] = 0;
        Atom** movingZAtoms;
        movingZAtoms = zconsMols[i]->getMyAtoms();
        for (int j = 0; j < zconsMols[i]->getNAtoms(); j++){
          movingZAtoms[j]->getFrc(force);
          fz[i] += force[whichDirection];
        }
      }
    }
    fzOut->writeFZ(info->getTime(), zconsMols.size(), &indexOfZConsMols[0], &fz[0],
                   &curZPos[0], &zPos[0]);
    curZconsTime += zconsTime;
  }

  zSysCOMVel = calcSysCOMVel();  
#ifdef IS_MPI
  if (worldRank == 0){
#endif
#ifdef IS_MPI
  }
#endif
}


template<typename T> double ZConstraint<T>::calcZSys(){
  //calculate reference z coordinate for z-constraint molecules
  double totalMass_local;
  double totalMass;
  double totalMZ_local;
  double totalMZ;
  double massOfCurMol;
  double COM[3];

  totalMass_local = 0;
  totalMZ_local = 0;

  for (int i = 0; i < nMols; i++){
    massOfCurMol = molecules[i].getTotalMass(); 
    molecules[i].getCOM(COM); 

    totalMass_local += massOfCurMol;
    totalMZ_local += massOfCurMol * COM[whichDirection];
  }


#ifdef IS_MPI  
  MPI_Allreduce(&totalMass_local, &totalMass, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&totalMZ_local, &totalMZ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  
#else
  totalMass = totalMass_local;
  totalMZ = totalMZ_local;
#endif  

  double zsys;
  zsys = totalMZ / totalMass;

  return zsys;
}

template<typename T> void ZConstraint<T>::thermalize(void){
  T::thermalize();
  zeroOutVel();
}

template<typename T> void ZConstraint<T>::zeroOutVel(){
  Atom** fixedZAtoms;  
  double COMvel[3];
  double vel[3];
  double zSysCOMVel;

  //zero out the velocities of center of mass of fixed z-constrained molecules

  for (int i = 0; i < (int) (zconsMols.size()); i++){
    if (states[i] == zcsFixed){
      zconsMols[i]->getCOMvel(COMvel);      
      //cout << "before resetting " << indexOfZConsMols[i] <<"'s vz is " << COMvel[whichDirection] << endl;

      fixedZAtoms = zconsMols[i]->getMyAtoms(); 

      for (int j = 0; j < zconsMols[i]->getNAtoms(); j++){
        fixedZAtoms[j]->getVel(vel);
        vel[whichDirection] -= COMvel[whichDirection];
        fixedZAtoms[j]->setVel(vel);
      }

      zconsMols[i]->getCOMvel(COMvel);
    }
  }

  zSysCOMVel = calcSysCOMVel();
#ifdef IS_MPI
  if (worldRank == 0){
#endif
#ifdef IS_MPI
  }
#endif

  // calculate the vz of center of mass of unconstrained molecules and moving z-constrained molecules
  double MVzOfMovingMols_local;
  double MVzOfMovingMols;
  double totalMassOfMovingZMols_local;
  double totalMassOfMovingZMols;

  MVzOfMovingMols_local = 0;
  totalMassOfMovingZMols_local = 0;

  for (int i = 0; i < (int) (unconsMols.size()); i++){
    unconsMols[i]->getCOMvel(COMvel);
    MVzOfMovingMols_local += massOfUnconsMols[i] * COMvel[whichDirection];
  } 

  for (int i = 0; i < (int) (zconsMols.size()); i++){
    if (states[i] == zcsMoving){
      zconsMols[i]->getCOMvel(COMvel);
      MVzOfMovingMols_local += massOfZConsMols[i] * COMvel[whichDirection];   
      totalMassOfMovingZMols_local += massOfZConsMols[i];
    }
  }

#ifndef IS_MPI
  MVzOfMovingMols = MVzOfMovingMols_local;
  totalMassOfMovingZMols = totalMassOfMovingZMols_local;
#else
  MPI_Allreduce(&MVzOfMovingMols_local, &MVzOfMovingMols, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&totalMassOfMovingZMols_local, &totalMassOfMovingZMols, 1,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  
#endif

  double vzOfMovingMols;
  vzOfMovingMols = MVzOfMovingMols /
                   (totalMassOfUncons + totalMassOfMovingZMols);

  //modify the velocites of unconstrained molecules  
  Atom** unconsAtoms;
  for (int i = 0; i < (int) (unconsMols.size()); i++){
    unconsAtoms = unconsMols[i]->getMyAtoms();
    for (int j = 0; j < unconsMols[i]->getNAtoms(); j++){
      unconsAtoms[j]->getVel(vel);
      vel[whichDirection] -= vzOfMovingMols;
      unconsAtoms[j]->setVel(vel);
    }
  }  

  //modify the velocities of moving z-constrained molecuels
  Atom** movingZAtoms;
  for (int i = 0; i < (int) (zconsMols.size()); i++){
    if (states[i] == zcsMoving){
      movingZAtoms = zconsMols[i]->getMyAtoms();
      for (int j = 0; j < zconsMols[i]->getNAtoms(); j++){
        movingZAtoms[j]->getVel(vel);
        vel[whichDirection] -= vzOfMovingMols;
        movingZAtoms[j]->setVel(vel);
      }
    }
  }


  zSysCOMVel = calcSysCOMVel();
#ifdef IS_MPI
  if (worldRank == 0){
#endif
#ifdef IS_MPI
  }
#endif
}


template<typename T> void ZConstraint<T>::doZconstraintForce(){
  Atom** zconsAtoms;
  double totalFZ; 
  double totalFZ_local;
  double COM[3];
  double force[3];

  //constrain the molecules which do not reach the specified positions  

  //Zero Out the force of z-contrained molecules    
  totalFZ_local = 0;

  //calculate the total z-contrained force of fixed z-contrained molecules

  for (int i = 0; i < (int) (zconsMols.size()); i++){
    if (states[i] == zcsFixed){
      zconsMols[i]->getCOM(COM);
      zconsAtoms = zconsMols[i]->getMyAtoms();  

      fz[i] = 0;      
      for (int j = 0; j < zconsMols[i]->getNAtoms(); j++){
        zconsAtoms[j]->getFrc(force);
        fz[i] += force[whichDirection];
      } 
      totalFZ_local += fz[i];

    }
  } 

  //calculate total z-constraint force
#ifdef IS_MPI
  MPI_Allreduce(&totalFZ_local, &totalFZ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  totalFZ = totalFZ_local;
#endif


  // apply negative to fixed z-constrained molecues;
  force[0] = 0;
  force[1] = 0;
  force[2] = 0;

  for (int i = 0; i < (int) (zconsMols.size()); i++){
    if (states[i] == zcsFixed){
      int nAtomOfCurZConsMol = zconsMols[i]->getNAtoms();
      zconsAtoms = zconsMols[i]->getMyAtoms();  

      for (int j = 0; j < nAtomOfCurZConsMol; j++){
        //force[whichDirection] = -fz[i]/ nAtomOfCurZConsMol;
        force[whichDirection] = -forcePolicy->getZFOfFixedZMols(zconsMols[i],
                                                                zconsAtoms[j],
                                                                fz[i]);
        zconsAtoms[j]->addFrc(force);
      }
    }
  } 

  force[0] = 0;
  force[1] = 0;
  force[2] = 0;

  //modify the forces of unconstrained molecules
  for (int i = 0; i < (int) (unconsMols.size()); i++){
    Atom** unconsAtoms = unconsMols[i]->getMyAtoms();

    for (int j = 0; j < unconsMols[i]->getNAtoms(); j++){
      //force[whichDirection] = totalFZ / (totNumOfUnconsAtoms + nMovingZAtoms);
      force[whichDirection] = forcePolicy->getZFOfMovingMols(unconsAtoms[j],
                                                             totalFZ);
      unconsAtoms[j]->addFrc(force);
    }
  }      

  //modify the forces of moving z-constrained molecules
  for (int i = 0; i < (int) (zconsMols.size()); i++){
    if (states[i] == zcsMoving){
      Atom** movingZAtoms = zconsMols[i]->getMyAtoms();     

      for (int j = 0; j < zconsMols[i]->getNAtoms(); j++){
        //force[whichDirection] = totalFZ / (totNumOfUnconsAtoms + nMovingZAtoms);
        force[whichDirection] = forcePolicy->getZFOfMovingMols(movingZAtoms[j],
                                                               totalFZ);
        movingZAtoms[j]->addFrc(force);
      }
    }
  }
}


template<typename T> void ZConstraint<T>::doHarmonic(vector<double>& resPos){
  double force[3];
  double harmonicU;
  double harmonicF;
  double COM[3];
  double diff;
  double totalFZ_local;
  double totalFZ;

  force[0] = 0;
  force[1] = 0;
  force[2] = 0;

  totalFZ_local = 0;

  for (int i = 0; i < (int) (zconsMols.size()); i++){
    if (states[i] == zcsMoving){
      zconsMols[i]->getCOM(COM);

      diff = COM[whichDirection] - resPos[i];

      harmonicU = 0.5 * kz[i] * diff * diff;  
      info->lrPot += harmonicU;

      harmonicF = -kz[i] * diff;
      totalFZ_local += harmonicF;

      //adjust force

      Atom** movingZAtoms = zconsMols[i]->getMyAtoms();     

      for (int j = 0; j < zconsMols[i]->getNAtoms(); j++){
        force[whichDirection] = forcePolicy->getHFOfFixedZMols(zconsMols[i],
                                                               movingZAtoms[j],
                                                               harmonicF);
        movingZAtoms[j]->addFrc(force);
      }
    }
  }

#ifndef IS_MPI
  totalFZ = totalFZ_local;
#else
  MPI_Allreduce(&totalFZ_local, &totalFZ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  
#endif

  force[0] = 0;
  force[1] = 0;
  force[2] = 0;

  //modify the forces of unconstrained molecules
  for (int i = 0; i < (int) (unconsMols.size()); i++){
    Atom** unconsAtoms = unconsMols[i]->getMyAtoms();

    for (int j = 0; j < unconsMols[i]->getNAtoms(); j++){
      //force[whichDirection] = - totalFZ /totNumOfUnconsAtoms;
      force[whichDirection] = -forcePolicy->getHFOfUnconsMols(unconsAtoms[j],
                                                              totalFZ);
      unconsAtoms[j]->addFrc(force);
    }
  }   

}

template<typename T> bool ZConstraint<T>::checkZConsState(){
  double COM[3];
  double diff;

  int changed_local;
  int changed;

  changed_local = 0;

  for (int i = 0; i < (int) (zconsMols.size()); i++){
    zconsMols[i]->getCOM(COM);
    diff = fabs(COM[whichDirection] - zPos[i]);  
    if (diff <= zconsTol && states[i] == zcsMoving){
      states[i] = zcsFixed;
      changed_local = 1;

      if(usingSMD)
        prevCantPos = cantPos;

      if (hasZConsGap)
        endFixTime[i] = info->getTime() + zconsFixTime;
    }
    else if (diff > zconsTol && states[i] == zcsFixed){
      states[i] = zcsMoving;
      changed_local = 1;   

      if(usingSMD)
         cantPos = prevCantPos;
      
      if (hasZConsGap)
        endFixTime[i] = INFINITE_TIME;
    }
  }

#ifndef IS_MPI
  changed = changed_local; 
#else
  MPI_Allreduce(&changed_local, &changed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  return (changed > 0);
}

template<typename T> bool ZConstraint<T>::haveFixedZMols(){
  int havingFixed_local;
  int havingFixed;

  havingFixed_local = 0;

  for (int i = 0; i < (int) (zconsMols.size()); i++)
    if (states[i] == zcsFixed){
      havingFixed_local = 1;
      break;
    }

#ifndef IS_MPI
  havingFixed = havingFixed_local;
#else
  MPI_Allreduce(&havingFixed_local, &havingFixed, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
#endif

  return (havingFixed > 0);
}


template<typename T> bool ZConstraint<T>::haveMovingZMols(){
  int havingMoving_local;
  int havingMoving;

  havingMoving_local = 0;

  for (int i = 0; i < (int) (zconsMols.size()); i++)
    if (states[i] == zcsMoving){
      havingMoving_local = 1;
      break;
    }

#ifndef IS_MPI
  havingMoving = havingMoving_local;
#else
  MPI_Allreduce(&havingMoving_local, &havingMoving, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
#endif

  return (havingMoving > 0);
}


template<typename T> double ZConstraint<T>::calcMovingMolsCOMVel(){
  double MVzOfMovingMols_local;
  double MVzOfMovingMols;
  double totalMassOfMovingZMols_local;
  double totalMassOfMovingZMols;
  double COMvel[3];

  MVzOfMovingMols_local = 0;
  totalMassOfMovingZMols_local = 0;

  for (int i = 0; i < unconsMols.size(); i++){
    unconsMols[i]->getCOMvel(COMvel);
    MVzOfMovingMols_local += massOfUnconsMols[i] * COMvel[whichDirection];
  } 

  for (int i = 0; i < zconsMols.size(); i++){
    if (states[i] == zcsMoving){
      zconsMols[i]->getCOMvel(COMvel);
      MVzOfMovingMols_local += massOfZConsMols[i] * COMvel[whichDirection];   
      totalMassOfMovingZMols_local += massOfZConsMols[i];
    }
  }

#ifndef IS_MPI
  MVzOfMovingMols = MVzOfMovingMols_local;
  totalMassOfMovingZMols = totalMassOfMovingZMols_local;
#else
  MPI_Allreduce(&MVzOfMovingMols_local, &MVzOfMovingMols, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&totalMassOfMovingZMols_local, &totalMassOfMovingZMols, 1,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  
#endif

  double vzOfMovingMols;
  vzOfMovingMols = MVzOfMovingMols /
                   (totalMassOfUncons + totalMassOfMovingZMols);

  return vzOfMovingMols;
}

template<typename T> double ZConstraint<T>::calcSysCOMVel(){
  double COMvel[3];
  double tempMVz_local;
  double tempMVz;
  double massOfZCons_local;
  double massOfZCons;


  tempMVz_local = 0;

  for (int i = 0 ; i < nMols; i++){
    molecules[i].getCOMvel(COMvel);
    tempMVz_local += molecules[i].getTotalMass() * COMvel[whichDirection];
  }

  massOfZCons_local = 0;

  for (int i = 0; i < (int) (massOfZConsMols.size()); i++){
    massOfZCons_local += massOfZConsMols[i];
  }
#ifndef IS_MPI
  massOfZCons = massOfZCons_local;
  tempMVz = tempMVz_local;
#else
  MPI_Allreduce(&massOfZCons_local, &massOfZCons, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&tempMVz_local, &tempMVz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  return tempMVz / (totalMassOfUncons + massOfZCons);
}

template<typename T> double ZConstraint<T>::calcTotalForce(){
  double force[3];  
  double totalForce_local;
  double totalForce;

  totalForce_local = 0;

  for (int i = 0; i < nAtoms; i++){
    atoms[i]->getFrc(force);
    totalForce_local += force[whichDirection];
  }

#ifndef IS_MPI
  totalForce = totalForce_local;
#else
  MPI_Allreduce(&totalForce_local, &totalForce, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
#endif

  return totalForce;
}

template<typename T> void ZConstraint<T>::PolicyByNumber::update(){
  //calculate the number of atoms of moving z-constrained molecules
  int nMovingZAtoms_local;
  int nMovingZAtoms;

  nMovingZAtoms_local = 0;
  for (int i = 0; i < (int) ((zconsIntegrator->zconsMols).size()); i++)
    if ((zconsIntegrator->states)[i] == (zconsIntegrator->zcsMoving)){
      nMovingZAtoms_local += (zconsIntegrator->zconsMols)[i]->getNAtoms();
    }

#ifdef IS_MPI
  MPI_Allreduce(&nMovingZAtoms_local, &nMovingZAtoms, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
#else
  nMovingZAtoms = nMovingZAtoms_local;
#endif
  totNumOfMovingAtoms = nMovingZAtoms + zconsIntegrator->totNumOfUnconsAtoms;
}

template<typename T> double ZConstraint<T>::PolicyByNumber::getZFOfFixedZMols(Molecule* mol,
                                                                              Atom* atom,
                                                                              double totalForce){
  return totalForce / mol->getNAtoms();
}

template<typename T> double ZConstraint<T>::PolicyByNumber::getZFOfMovingMols(Atom* atom,
                                                                              double totalForce){
  return totalForce / totNumOfMovingAtoms;
}

template<typename T> double ZConstraint<T>::PolicyByNumber::getHFOfFixedZMols(Molecule* mol,
                                                                              Atom* atom,
                                                                              double totalForce){
  return totalForce / mol->getNAtoms();
}

template<typename T> double ZConstraint<T>::PolicyByNumber::getHFOfUnconsMols(Atom* atom,
                                                                              double totalForce){
  return totalForce / zconsIntegrator->totNumOfUnconsAtoms;
}


template<typename T> void ZConstraint<T>::PolicyByMass::update(){
  //calculate the number of atoms of moving z-constrained molecules
  double massOfMovingZAtoms_local;
  double massOfMovingZAtoms;

  massOfMovingZAtoms_local = 0;
  for (int i = 0; i < (int) ((zconsIntegrator->zconsMols).size()); i++)
    if ((zconsIntegrator->states)[i] == (zconsIntegrator->zcsMoving)){
      massOfMovingZAtoms_local += (zconsIntegrator->zconsMols)[i]->getTotalMass();
    }

#ifdef IS_MPI
  MPI_Allreduce(&massOfMovingZAtoms_local, &massOfMovingZAtoms, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
#else
  massOfMovingZAtoms = massOfMovingZAtoms_local;
#endif
  totMassOfMovingAtoms = massOfMovingZAtoms +
                         zconsIntegrator->totalMassOfUncons;
}

template<typename T> double ZConstraint<T>::PolicyByMass::getZFOfFixedZMols(Molecule* mol,
                                                                            Atom* atom,
                                                                            double totalForce){
  return totalForce * atom->getMass() / mol->getTotalMass();
}

template<typename T> double ZConstraint<T>::PolicyByMass::getZFOfMovingMols(Atom* atom,
                                                                            double totalForce){
  return totalForce * atom->getMass() / totMassOfMovingAtoms;
}

template<typename T> double ZConstraint<T>::PolicyByMass::getHFOfFixedZMols(Molecule* mol,
                                                                            Atom* atom,
                                                                            double totalForce){
  return totalForce * atom->getMass() / mol->getTotalMass();
}

template<typename T> double ZConstraint<T>::PolicyByMass::getHFOfUnconsMols(Atom* atom,
                                                                            double totalForce){
  return totalForce * atom->getMass() / zconsIntegrator->totalMassOfUncons;
}

template<typename T> void ZConstraint<T>::updateZPos(){
  double curTime;
  double COM[3];
  
  curTime = info->getTime();

  for (size_t i = 0; i < zconsMols.size(); i++){

    if (states[i] == zcsFixed && curTime >= endFixTime[i]){
      zPos[i] += zconsGap;

      if (usingSMD){
        zconsMols[i]->getCOM(COM);
        cantPos[i] = COM[whichDirection];
      }
      
    }
    
  }
  
}

template<typename T> void ZConstraint<T>::updateCantPos(){
  double curTime;
  double dt;

  curTime = info->getTime();
  dt = info->dt;

  for (size_t i = 0; i < zconsMols.size(); i++){
    if (states[i] == zcsMoving){
      cantPos[i] += cantVel[i] * dt;
    }
  }

}
