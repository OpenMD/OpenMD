#include <iostream>
#include <stdlib.h>
#include <math.h>
#ifdef IS_MPI
#include "mpiSimulation.hpp"
#include <unistd.h>
#endif //is_mpi

#ifdef PROFILE
#include "mdProfile.hpp"
#endif // profile

#include "Integrator.hpp"
#include "simError.h"


template<typename T> Integrator<T>::Integrator(SimInfo* theInfo,
                                               ForceFields* the_ff){
  info = theInfo;
  myFF = the_ff;
  isFirst = 1;

  molecules = info->molecules;
  nMols = info->n_mol;

  // give a little love back to the SimInfo object

  if (info->the_integrator != NULL){
    delete info->the_integrator;
  }

  nAtoms = info->n_atoms;
  integrableObjects = info->integrableObjects;


  // check for constraints

  constrainedA = NULL;
  constrainedB = NULL;
  constrainedDsqr = NULL;
  moving = NULL;
  moved = NULL;
  oldPos = NULL;

  nConstrained = 0;

  checkConstraints();

}

template<typename T> Integrator<T>::~Integrator(){

  if (nConstrained){
    delete[] constrainedA;
    delete[] constrainedB;
    delete[] constrainedDsqr;
    delete[] moving;
    delete[] moved;
    delete[] oldPos;
  }

}


template<typename T> void Integrator<T>::checkConstraints(void){
  isConstrained = 0;

  Constraint* temp_con;
  Constraint* dummy_plug;
  temp_con = new Constraint[info->n_SRI];
  nConstrained = 0;
  int constrained = 0;

  SRI** theArray;
  for (int i = 0; i < nMols; i++){

	  theArray = (SRI * *) molecules[i].getMyBonds();
    for (int j = 0; j < molecules[i].getNBonds(); j++){
      constrained = theArray[j]->is_constrained();

      if (constrained){
        dummy_plug = theArray[j]->get_constraint();
        temp_con[nConstrained].set_a(dummy_plug->get_a());
        temp_con[nConstrained].set_b(dummy_plug->get_b());
        temp_con[nConstrained].set_dsqr(dummy_plug->get_dsqr());

        nConstrained++;
        constrained = 0;
      }
    }

    theArray = (SRI * *) molecules[i].getMyBends();
    for (int j = 0; j < molecules[i].getNBends(); j++){
      constrained = theArray[j]->is_constrained();

      if (constrained){
        dummy_plug = theArray[j]->get_constraint();
        temp_con[nConstrained].set_a(dummy_plug->get_a());
        temp_con[nConstrained].set_b(dummy_plug->get_b());
        temp_con[nConstrained].set_dsqr(dummy_plug->get_dsqr());

        nConstrained++;
        constrained = 0;
      }
    }

    theArray = (SRI * *) molecules[i].getMyTorsions();
    for (int j = 0; j < molecules[i].getNTorsions(); j++){
      constrained = theArray[j]->is_constrained();

      if (constrained){
        dummy_plug = theArray[j]->get_constraint();
        temp_con[nConstrained].set_a(dummy_plug->get_a());
        temp_con[nConstrained].set_b(dummy_plug->get_b());
        temp_con[nConstrained].set_dsqr(dummy_plug->get_dsqr());

        nConstrained++;
        constrained = 0;
      }
    }
  }


  if (nConstrained > 0){
    isConstrained = 1;

    if (constrainedA != NULL)
      delete[] constrainedA;
    if (constrainedB != NULL)
      delete[] constrainedB;
    if (constrainedDsqr != NULL)
      delete[] constrainedDsqr;

    constrainedA = new int[nConstrained];
    constrainedB = new int[nConstrained];
    constrainedDsqr = new double[nConstrained];

    for (int i = 0; i < nConstrained; i++){
      constrainedA[i] = temp_con[i].get_a();
      constrainedB[i] = temp_con[i].get_b();
      constrainedDsqr[i] = temp_con[i].get_dsqr();
    }


    // save oldAtoms to check for lode balancing later on.

    oldAtoms = nAtoms;

    moving = new int[nAtoms];
    moved = new int[nAtoms];

    oldPos = new double[nAtoms * 3];
  }

  delete[] temp_con;
}


template<typename T> void Integrator<T>::integrate(void){

  double runTime = info->run_time;
  double sampleTime = info->sampleTime;
  double statusTime = info->statusTime;
  double thermalTime = info->thermalTime;
  double resetTime = info->resetTime;

  double difference;
  double currSample;
  double currThermal;
  double currStatus;
  double currReset;

  int calcPot, calcStress;

  tStats = new Thermo(info);
  statOut = new StatWriter(info);
  dumpOut = new DumpWriter(info);

  atoms = info->atoms;

  dt = info->dt;
  dt2 = 0.5 * dt;

  readyCheck();

  // remove center of mass drift velocity (in case we passed in a configuration
  // that was drifting
  tStats->removeCOMdrift();

  // initialize the retraints if necessary
  if (info->useSolidThermInt && !info->useLiquidThermInt) {
    myFF->initRestraints();
  }

  // initialize the forces before the first step

  calcForce(1, 1);

  //execute constraint algorithm to make sure at the very beginning the system is constrained  
  if(nConstrained){
    preMove();
    constrainA();
    calcForce(1, 1);
    constrainB();
  }

  if (info->setTemp){
    thermalize();
  }

  calcPot     = 0;
  calcStress  = 0;
  currSample  = sampleTime + info->getTime();
  currThermal = thermalTime+ info->getTime();
  currStatus  = statusTime + info->getTime();
  currReset   = resetTime  + info->getTime();

  dumpOut->writeDump(info->getTime());
  statOut->writeStat(info->getTime());


#ifdef IS_MPI
  strcpy(checkPointMsg, "The integrator is ready to go.");
  MPIcheckPoint();
#endif // is_mpi

  while (info->getTime() < runTime && !stopIntegrator()){
    difference = info->getTime() + dt - currStatus;
    if (difference > 0 || fabs(difference) < 1e-4 ){
      calcPot = 1;
      calcStress = 1;
    }

#ifdef PROFILE
    startProfile( pro1 );
#endif
    
    integrateStep(calcPot, calcStress);

#ifdef PROFILE
    endProfile( pro1 );

    startProfile( pro2 );
#endif // profile

    info->incrTime(dt);

    if (info->setTemp){
      if (info->getTime() >= currThermal){
        thermalize();
        currThermal += thermalTime;
      }
    }

    if (info->getTime() >= currSample){
      dumpOut->writeDump(info->getTime());
      currSample += sampleTime;
    }

    if (info->getTime() >= currStatus){
      statOut->writeStat(info->getTime());
      calcPot = 0;
      calcStress = 0;
      currStatus += statusTime;
    }

    if (info->resetIntegrator){
      if (info->getTime() >= currReset){
        this->resetIntegrator();
        currReset += resetTime;
      }
    }
    
#ifdef PROFILE
    endProfile( pro2 );
#endif //profile

#ifdef IS_MPI
    strcpy(checkPointMsg, "successfully took a time step.");
    MPIcheckPoint();
#endif // is_mpi
  }

  dumpOut->writeFinal(info->getTime());

  // dump out a file containing the omega values for the final configuration
  if (info->useSolidThermInt && !info->useLiquidThermInt)
    myFF->dumpzAngle();
  

  delete dumpOut;
  delete statOut;
}

template<typename T> void Integrator<T>::integrateStep(int calcPot,
                                                       int calcStress){
  // Position full step, and velocity half step

#ifdef PROFILE
  startProfile(pro3);
#endif //profile

  //save old state (position, velocity etc)
  preMove();
#ifdef PROFILE
  endProfile(pro3);

  startProfile(pro4);
#endif // profile

  moveA();

#ifdef PROFILE
  endProfile(pro4);
  
  startProfile(pro5);
#endif//profile


#ifdef IS_MPI
  strcpy(checkPointMsg, "Succesful moveA\n");
  MPIcheckPoint();
#endif // is_mpi

  // calc forces
  calcForce(calcPot, calcStress);

#ifdef IS_MPI
  strcpy(checkPointMsg, "Succesful doForces\n");
  MPIcheckPoint();
#endif // is_mpi

#ifdef PROFILE
  endProfile( pro5 );

  startProfile( pro6 );
#endif //profile

  // finish the velocity  half step

  moveB();

#ifdef PROFILE
  endProfile(pro6);
#endif // profile

#ifdef IS_MPI
  strcpy(checkPointMsg, "Succesful moveB\n");
  MPIcheckPoint();
#endif // is_mpi
}


template<typename T> void Integrator<T>::moveA(void){
  size_t i, j;
  DirectionalAtom* dAtom;
  double Tb[3], ji[3];
  double vel[3], pos[3], frc[3];
  double mass;
  double omega;
 
  for (i = 0; i < integrableObjects.size() ; i++){
    integrableObjects[i]->getVel(vel);
    integrableObjects[i]->getPos(pos);
    integrableObjects[i]->getFrc(frc);
    
    mass = integrableObjects[i]->getMass();

    for (j = 0; j < 3; j++){
      // velocity half step
      vel[j] += (dt2 * frc[j] / mass) * eConvert;
      // position whole step
      pos[j] += dt * vel[j];
    }

    integrableObjects[i]->setVel(vel);
    integrableObjects[i]->setPos(pos);

    if (integrableObjects[i]->isDirectional()){

      // get and convert the torque to body frame

      integrableObjects[i]->getTrq(Tb);
      integrableObjects[i]->lab2Body(Tb);

      // get the angular momentum, and propagate a half step

      integrableObjects[i]->getJ(ji);

      for (j = 0; j < 3; j++)
        ji[j] += (dt2 * Tb[j]) * eConvert;

      this->rotationPropagation( integrableObjects[i], ji );

      integrableObjects[i]->setJ(ji);
    }
  }

  if(nConstrained)
    constrainA();
}


template<typename T> void Integrator<T>::moveB(void){
  int i, j;
  double Tb[3], ji[3];
  double vel[3], frc[3];
  double mass;

  for (i = 0; i < integrableObjects.size(); i++){
    integrableObjects[i]->getVel(vel);
    integrableObjects[i]->getFrc(frc);

    mass = integrableObjects[i]->getMass();

    // velocity half step
    for (j = 0; j < 3; j++)
      vel[j] += (dt2 * frc[j] / mass) * eConvert;

    integrableObjects[i]->setVel(vel);

    if (integrableObjects[i]->isDirectional()){

      // get and convert the torque to body frame

      integrableObjects[i]->getTrq(Tb);
      integrableObjects[i]->lab2Body(Tb);

      // get the angular momentum, and propagate a half step

      integrableObjects[i]->getJ(ji);

      for (j = 0; j < 3; j++)
        ji[j] += (dt2 * Tb[j]) * eConvert;


      integrableObjects[i]->setJ(ji);
    }
  }

  if(nConstrained)
    constrainB();
}


template<typename T> void Integrator<T>::preMove(void){
  int i, j;
  double pos[3];

  if (nConstrained){
    for (i = 0; i < nAtoms; i++){
      atoms[i]->getPos(pos);

      for (j = 0; j < 3; j++){
        oldPos[3 * i + j] = pos[j];
      }
    }
  }
}

template<typename T> void Integrator<T>::constrainA(){
  int i, j;
  int done;
  double posA[3], posB[3];
  double velA[3], velB[3];
  double pab[3];
  double rab[3];
  int a, b, ax, ay, az, bx, by, bz;
  double rma, rmb;
  double dx, dy, dz;
  double rpab;
  double rabsq, pabsq, rpabsq;
  double diffsq;
  double gab;
  int iteration;

  for (i = 0; i < nAtoms; i++){
    moving[i] = 0;
    moved[i] = 1;
  }

  iteration = 0;
  done = 0;
  while (!done && (iteration < maxIteration)){
    done = 1;
    for (i = 0; i < nConstrained; i++){
      a = constrainedA[i];
      b = constrainedB[i];

      ax = (a * 3) + 0;
      ay = (a * 3) + 1;
      az = (a * 3) + 2;

      bx = (b * 3) + 0;
      by = (b * 3) + 1;
      bz = (b * 3) + 2;

      if (moved[a] || moved[b]){
        atoms[a]->getPos(posA);
        atoms[b]->getPos(posB);

        for (j = 0; j < 3; j++)
          pab[j] = posA[j] - posB[j];

        //periodic boundary condition

        info->wrapVector(pab);

        pabsq = pab[0] * pab[0] + pab[1] * pab[1] + pab[2] * pab[2];

        rabsq = constrainedDsqr[i];
        diffsq = rabsq - pabsq;

        // the original rattle code from alan tidesley
        if (fabs(diffsq) > (tol * rabsq * 2)){
          rab[0] = oldPos[ax] - oldPos[bx];
          rab[1] = oldPos[ay] - oldPos[by];
          rab[2] = oldPos[az] - oldPos[bz];

          info->wrapVector(rab);

          rpab = rab[0] * pab[0] + rab[1] * pab[1] + rab[2] * pab[2];

          rpabsq = rpab * rpab;


          if (rpabsq < (rabsq * -diffsq)){
#ifdef IS_MPI
            a = atoms[a]->getGlobalIndex();
            b = atoms[b]->getGlobalIndex();
#endif //is_mpi
            sprintf(painCave.errMsg,
                    "Constraint failure in constrainA at atom %d and %d.\n", a,
                    b);
            painCave.isFatal = 1;
            simError();
          }

          rma = 1.0 / atoms[a]->getMass();
          rmb = 1.0 / atoms[b]->getMass();

          gab = diffsq / (2.0 * (rma + rmb) * rpab);

          dx = rab[0] * gab;
          dy = rab[1] * gab;
          dz = rab[2] * gab;

          posA[0] += rma * dx;
          posA[1] += rma * dy;
          posA[2] += rma * dz;

          atoms[a]->setPos(posA);

          posB[0] -= rmb * dx;
          posB[1] -= rmb * dy;
          posB[2] -= rmb * dz;

          atoms[b]->setPos(posB);

          dx = dx / dt;
          dy = dy / dt;
          dz = dz / dt;

          atoms[a]->getVel(velA);

          velA[0] += rma * dx;
          velA[1] += rma * dy;
          velA[2] += rma * dz;

          atoms[a]->setVel(velA);

          atoms[b]->getVel(velB);

          velB[0] -= rmb * dx;
          velB[1] -= rmb * dy;
          velB[2] -= rmb * dz;

          atoms[b]->setVel(velB);

          moving[a] = 1;
          moving[b] = 1;
          done = 0;
        }
      }
    }

    for (i = 0; i < nAtoms; i++){
      moved[i] = moving[i];
      moving[i] = 0;
    }

    iteration++;
  }

  if (!done){
    sprintf(painCave.errMsg,
            "Constraint failure in constrainA, too many iterations: %d\n",
            iteration);
    painCave.isFatal = 1;
    simError();
  }

}

template<typename T> void Integrator<T>::constrainB(void){
  int i, j;
  int done;
  double posA[3], posB[3];
  double velA[3], velB[3];
  double vxab, vyab, vzab;
  double rab[3];
  int a, b, ax, ay, az, bx, by, bz;
  double rma, rmb;
  double dx, dy, dz;
  double rvab;
  double gab;
  int iteration;

  for (i = 0; i < nAtoms; i++){
    moving[i] = 0;
    moved[i] = 1;
  }

  done = 0;
  iteration = 0;
  while (!done && (iteration < maxIteration)){
    done = 1;

    for (i = 0; i < nConstrained; i++){
      a = constrainedA[i];
      b = constrainedB[i];

      ax = (a * 3) + 0;
      ay = (a * 3) + 1;
      az = (a * 3) + 2;

      bx = (b * 3) + 0;
      by = (b * 3) + 1;
      bz = (b * 3) + 2;

      if (moved[a] || moved[b]){
        atoms[a]->getVel(velA);
        atoms[b]->getVel(velB);

        vxab = velA[0] - velB[0];
        vyab = velA[1] - velB[1];
        vzab = velA[2] - velB[2];

        atoms[a]->getPos(posA);
        atoms[b]->getPos(posB);

        for (j = 0; j < 3; j++)
          rab[j] = posA[j] - posB[j];

        info->wrapVector(rab);

        rma = 1.0 / atoms[a]->getMass();
        rmb = 1.0 / atoms[b]->getMass();

        rvab = rab[0] * vxab + rab[1] * vyab + rab[2] * vzab;

        gab = -rvab / ((rma + rmb) * constrainedDsqr[i]);

        if (fabs(gab) > tol){
          dx = rab[0] * gab;
          dy = rab[1] * gab;
          dz = rab[2] * gab;

          velA[0] += rma * dx;
          velA[1] += rma * dy;
          velA[2] += rma * dz;

          atoms[a]->setVel(velA);

          velB[0] -= rmb * dx;
          velB[1] -= rmb * dy;
          velB[2] -= rmb * dz;

          atoms[b]->setVel(velB);

          moving[a] = 1;
          moving[b] = 1;
          done = 0;
        }
      }
    }

    for (i = 0; i < nAtoms; i++){
      moved[i] = moving[i];
      moving[i] = 0;
    }

    iteration++;
  }

  if (!done){
    sprintf(painCave.errMsg,
            "Constraint failure in constrainB, too many iterations: %d\n",
            iteration);
    painCave.isFatal = 1;
    simError();
  }
}

template<typename T> void Integrator<T>::rotationPropagation
( StuntDouble* sd, double ji[3] ){

  double angle;
  double A[3][3], I[3][3];
  int i, j, k;

  // use the angular velocities to propagate the rotation matrix a
  // full time step

  sd->getA(A);
  sd->getI(I);

  if (sd->isLinear()) {
    i = sd->linearAxis();
    j = (i+1)%3;
    k = (i+2)%3;
    
    angle = dt2 * ji[j] / I[j][j];
    this->rotate( k, i, angle, ji, A );

    angle = dt * ji[k] / I[k][k];
    this->rotate( i, j, angle, ji, A);

    angle = dt2 * ji[j] / I[j][j];
    this->rotate( k, i, angle, ji, A );

  } else {
    // rotate about the x-axis
    angle = dt2 * ji[0] / I[0][0];
    this->rotate( 1, 2, angle, ji, A );
    
    // rotate about the y-axis
    angle = dt2 * ji[1] / I[1][1];
    this->rotate( 2, 0, angle, ji, A );
    
    // rotate about the z-axis
    angle = dt * ji[2] / I[2][2];
    sd->addZangle(angle);
    this->rotate( 0, 1, angle, ji, A);
    
    // rotate about the y-axis
    angle = dt2 * ji[1] / I[1][1];
    this->rotate( 2, 0, angle, ji, A );
    
    // rotate about the x-axis
    angle = dt2 * ji[0] / I[0][0];
    this->rotate( 1, 2, angle, ji, A );
    
  }
  sd->setA( A  );
}

template<typename T> void Integrator<T>::rotate(int axes1, int axes2,
                                                double angle, double ji[3],
                                                double A[3][3]){
  int i, j, k;
  double sinAngle;
  double cosAngle;
  double angleSqr;
  double angleSqrOver4;
  double top, bottom;
  double rot[3][3];
  double tempA[3][3];
  double tempJ[3];

  // initialize the tempA

  for (i = 0; i < 3; i++){
    for (j = 0; j < 3; j++){
      tempA[j][i] = A[i][j];
    }
  }

  // initialize the tempJ

  for (i = 0; i < 3; i++)
    tempJ[i] = ji[i];

  // initalize rot as a unit matrix

  rot[0][0] = 1.0;
  rot[0][1] = 0.0;
  rot[0][2] = 0.0;

  rot[1][0] = 0.0;
  rot[1][1] = 1.0;
  rot[1][2] = 0.0;

  rot[2][0] = 0.0;
  rot[2][1] = 0.0;
  rot[2][2] = 1.0;

  // use a small angle aproximation for sin and cosine

  angleSqr = angle * angle;
  angleSqrOver4 = angleSqr / 4.0;
  top = 1.0 - angleSqrOver4;
  bottom = 1.0 + angleSqrOver4;

  cosAngle = top / bottom;
  sinAngle = angle / bottom;

  rot[axes1][axes1] = cosAngle;
  rot[axes2][axes2] = cosAngle;

  rot[axes1][axes2] = sinAngle;
  rot[axes2][axes1] = -sinAngle;

  // rotate the momentum acoording to: ji[] = rot[][] * ji[]

  for (i = 0; i < 3; i++){
    ji[i] = 0.0;
    for (k = 0; k < 3; k++){
      ji[i] += rot[i][k] * tempJ[k];
    }
  }

  // rotate the Rotation matrix acording to:
  //            A[][] = A[][] * transpose(rot[][])


  // NOte for as yet unknown reason, we are performing the
  // calculation as:
  //                transpose(A[][]) = transpose(A[][]) * transpose(rot[][])

  for (i = 0; i < 3; i++){
    for (j = 0; j < 3; j++){
      A[j][i] = 0.0;
      for (k = 0; k < 3; k++){
        A[j][i] += tempA[i][k] * rot[j][k];
      }
    }
  }
}

template<typename T> void Integrator<T>::calcForce(int calcPot, int calcStress){
  myFF->doForces(calcPot, calcStress);
}

template<typename T> void Integrator<T>::thermalize(){
  tStats->velocitize();
}

template<typename T> double Integrator<T>::getConservedQuantity(void){
  return tStats->getTotalE();
}
template<typename T> string Integrator<T>::getAdditionalParameters(void){
  //By default, return a null string
  //The reason we use string instead of char* is that if we use char*, we will
  //return a pointer point to local variable which might cause problem
  return string();
}
