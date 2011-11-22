/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include <cmath>


#include "io/StatWriter.hpp"
#include "minimizers/Minimizer.hpp"
#include "primitives/Molecule.hpp"
#ifdef IS_MPI
#include <mpi.h>
#endif
namespace OpenMD {
  RealType dotProduct(const std::vector<RealType>& v1, const std::vector<RealType>& v2) {
    if (v1.size() != v2.size()) {

    }


    RealType result = 0.0;    
    for (unsigned int i = 0; i < v1.size(); ++i) {
      result += v1[i] * v2[i];
    }

    return result;
  }

  Minimizer::Minimizer(SimInfo* rhs) :
    info(rhs), usingShake(false) {

      forceMan = new ForceManager(info);
      paramSet= new MinimizerParameterSet(info), calcDim();
      curX = getCoor();
      curG.resize(ndim);

    }

  Minimizer::~Minimizer() {
    delete forceMan;
    delete paramSet;
  }

  void Minimizer::calcEnergyGradient(std::vector<RealType> &x,
				     std::vector<RealType> &grad, RealType&energy, int&status) {

    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* integrableObject;    
    std::vector<RealType> myGrad;    
    int shakeStatus;

    status = 1;

    setCoor(x);

    if (usingShake) {
      shakeStatus = shakeR();
    }

    energy = calcPotential();

    if (usingShake) {
      shakeStatus = shakeF();
    }

    x = getCoor();

    int index = 0;

    for (mol = info->beginMolecule(i); mol != NULL; mol = info->nextMolecule(i)) {
      for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL;
	   integrableObject = mol->nextIntegrableObject(j)) {

	myGrad = integrableObject->getGrad();
	for (unsigned int k = 0; k < myGrad.size(); ++k) {

	  grad[index++] = myGrad[k];
	}
      }            
    }

  }

  void Minimizer::setCoor(std::vector<RealType> &x) {
    Vector3d position;
    Vector3d eulerAngle;
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* integrableObject;    
    int index = 0;

    for (mol = info->beginMolecule(i); mol != NULL; mol = info->nextMolecule(i)) {
      for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL;
	   integrableObject = mol->nextIntegrableObject(j)) {

	position[0] = x[index++];
	position[1] = x[index++];
	position[2] = x[index++];

	integrableObject->setPos(position);

	if (integrableObject->isDirectional()) {
	  eulerAngle[0] = x[index++];
	  eulerAngle[1] = x[index++];
	  eulerAngle[2] = x[index++];

	  integrableObject->setEuler(eulerAngle);
	}
      }
    }
    
  }

  std::vector<RealType> Minimizer::getCoor() {
    Vector3d position;
    Vector3d eulerAngle;
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* integrableObject;    
    int index = 0;
    std::vector<RealType> x(getDim());

    for (mol = info->beginMolecule(i); mol != NULL; mol = info->nextMolecule(i)) {
      for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL;
	   integrableObject = mol->nextIntegrableObject(j)) {
                
	position = integrableObject->getPos();
	x[index++] = position[0];
	x[index++] = position[1];
	x[index++] = position[2];

	if (integrableObject->isDirectional()) {
	  eulerAngle = integrableObject->getEuler();
	  x[index++] = eulerAngle[0];
	  x[index++] = eulerAngle[1];
	  x[index++] = eulerAngle[2];
	}
      }
    }
    return x;
  }


  /*
    int Minimizer::shakeR() {
    int    i,       j;

    int    done;

    RealType posA[3], posB[3];

    RealType velA[3], velB[3];

    RealType pab[3];

    RealType rab[3];

    int    a,       b,
    ax,      ay,
    az,      bx,
    by,      bz;

    RealType rma,     rmb;

    RealType dx,      dy,
    dz;

    RealType rpab;

    RealType rabsq,   pabsq,
    rpabsq;

    RealType diffsq;

    RealType gab;

    int    iteration;

    for(i = 0; i < nAtoms; i++) {
    moving[i] = 0;

    moved[i] = 1;
    }

    iteration = 0;

    done = 0;

    while (!done && (iteration < maxIteration)) {
    done = 1;

    for(i = 0; i < nConstrained; i++) {
    a = constrainedA[i];

    b = constrainedB[i];

    ax = (a * 3) + 0;

    ay = (a * 3) + 1;

    az = (a * 3) + 2;

    bx = (b * 3) + 0;

    by = (b * 3) + 1;

    bz = (b * 3) + 2;

    if (moved[a] || moved[b]) {
    posA = atoms[a]->getPos();

    posB = atoms[b]->getPos();

    for(j = 0; j < 3; j++)
    pab[j] = posA[j] - posB[j];

    //periodic boundary condition

    info->wrapVector(pab);

    pabsq = pab[0] * pab[0] + pab[1] * pab[1] + pab[2] * pab[2];

    rabsq = constrainedDsqr[i];

    diffsq = rabsq - pabsq;

    // the original rattle code from alan tidesley

    if (fabs(diffsq) > (tol * rabsq * 2)) {
    rab[0] = oldPos[ax] - oldPos[bx];

    rab[1] = oldPos[ay] - oldPos[by];

    rab[2] = oldPos[az] - oldPos[bz];

    info->wrapVector(rab);

    rpab = rab[0] * pab[0] + rab[1] * pab[1] + rab[2] * pab[2];

    rpabsq = rpab * rpab;

    if (rpabsq < (rabsq * -diffsq)) {

    #ifdef IS_MPI

    a = atoms[a]->getGlobalIndex();

    b = atoms[b]->getGlobalIndex();

    #endif //is_mpi

    //std::cerr << "Waring: constraint failure" << std::endl;

    gab = sqrt(rabsq / pabsq);

    rab[0] = (posA[0] - posB[0])
    * gab;

    rab[1] = (posA[1] - posB[1])
    * gab;

    rab[2] = (posA[2] - posB[2])
    * gab;

    info->wrapVector(rab);

    rpab =
    rab[0] * pab[0] + rab[1] * pab[1] + rab[2] * pab[2];
    }

    //rma = 1.0 / atoms[a]->getMass();

    //rmb = 1.0 / atoms[b]->getMass();

    rma = 1.0;

    rmb = 1.0;

    gab = diffsq / (2.0 * (rma + rmb) * rpab);

    dx = rab[0]*
    gab;

    dy = rab[1]*
    gab;

    dz = rab[2]*
    gab;

    posA[0] += rma *dx;

    posA[1] += rma *dy;

    posA[2] += rma *dz;

    atoms[a]->setPos(posA);

    posB[0] -= rmb *dx;

    posB[1] -= rmb *dy;

    posB[2] -= rmb *dz;

    atoms[b]->setPos(posB);

    moving[a] = 1;

    moving[b] = 1;

    done = 0;
    }
    }
    }

    for(i = 0; i < nAtoms; i++) {
    moved[i] = moving[i];

    moving[i] = 0;
    }

    iteration++;
    }

    if (!done) {
    std::cerr << "Waring: can not constraint within maxIteration"
    << std::endl;

    return -1;
    } else
    return 1;
    }

    //remove constraint force along the bond direction


    int Minimizer::shakeF() {
    int    i,       j;

    int    done;

    RealType posA[3], posB[3];

    RealType frcA[3], frcB[3];

    RealType rab[3],  fpab[3];

    int    a,       b,
    ax,      ay,
    az,      bx,
    by,      bz;

    RealType rma,     rmb;

    RealType rvab;

    RealType gab;

    RealType rabsq;

    RealType rfab;

    int    iteration;

    for(i = 0; i < nAtoms; i++) {
    moving[i] = 0;

    moved[i] = 1;
    }

    done = 0;

    iteration = 0;

    while (!done && (iteration < maxIteration)) {
    done = 1;

    for(i = 0; i < nConstrained; i++) {
    a = constrainedA[i];

    b = constrainedB[i];

    ax = (a * 3) + 0;

    ay = (a * 3) + 1;

    az = (a * 3) + 2;

    bx = (b * 3) + 0;

    by = (b * 3) + 1;

    bz = (b * 3) + 2;

    if (moved[a] || moved[b]) {
    posA = atoms[a]->getPos();

    posB = atoms[b]->getPos();

    for(j = 0; j < 3; j++)
    rab[j] = posA[j] - posB[j];

    info->wrapVector(rab);

    atoms[a]->getFrc(frcA);

    atoms[b]->getFrc(frcB);

    //rma = 1.0 / atoms[a]->getMass();

    //rmb = 1.0 / atoms[b]->getMass();

    rma = 1.0;

    rmb = 1.0;

    fpab[0] = frcA[0] * rma - frcB[0] * rmb;

    fpab[1] = frcA[1] * rma - frcB[1] * rmb;

    fpab[2] = frcA[2] * rma - frcB[2] * rmb;

    gab = fpab[0] * fpab[0] + fpab[1] * fpab[1] + fpab[2] * fpab[2];

    if (gab < 1.0)
    gab = 1.0;

    rabsq = rab[0] * rab[0] + rab[1] * rab[1] + rab[2] * rab[2];

    rfab = rab[0] * fpab[0] + rab[1] * fpab[1] + rab[2] * fpab[2];

    if (fabs(rfab) > sqrt(rabsq*gab) * 0.00001) {
    gab = -rfab / (rabsq * (rma + rmb));

    frcA[0] = rab[0]*
    gab;

    frcA[1] = rab[1]*
    gab;

    frcA[2] = rab[2]*
    gab;

    atoms[a]->addFrc(frcA);

    frcB[0] = -rab[0]*gab;

    frcB[1] = -rab[1]*gab;

    frcB[2] = -rab[2]*gab;

    atoms[b]->addFrc(frcB);

    moving[a] = 1;

    moving[b] = 1;

    done = 0;
    }
    }
    }

    for(i = 0; i < nAtoms; i++) {
    moved[i] = moving[i];

    moving[i] = 0;
    }

    iteration++;
    }

    if (!done) {
    std::cerr << "Waring: can not constraint within maxIteration"
    << std::endl;

    return -1;
    } else
    return 1;
    }

  */
    
  //calculate the value of object function

  void Minimizer::calcF() {
    calcEnergyGradient(curX, curG, curF, egEvalStatus);
  }

  void Minimizer::calcF(std::vector < RealType > &x, RealType&f, int&status) {
    std::vector < RealType > tempG;

    tempG.resize(x.size());

    calcEnergyGradient(x, tempG, f, status);
  }

  //calculate the gradient

  void Minimizer::calcG() {
    calcEnergyGradient(curX, curG, curF, egEvalStatus);
  }

  void Minimizer::calcG(std::vector<RealType>& x, std::vector<RealType>& g, RealType&f, int&status) {
    calcEnergyGradient(x, g, f, status);
  }

  void Minimizer::calcDim() {

    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* integrableObject;    
    ndim = 0;

    for (mol = info->beginMolecule(i); mol != NULL; mol = info->nextMolecule(i)) {
      for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL;
	   integrableObject = mol->nextIntegrableObject(j)) {

	ndim += 3;

	if (integrableObject->isDirectional()) {
	  ndim += 3;
	}
      }

    }
  }

  void Minimizer::setX(std::vector < RealType > &x) {
    if (x.size() != ndim) {
      sprintf(painCave.errMsg, "Minimizer Error: dimesion of x and curX does not match\n");
      painCave.isFatal = 1;
      simError();
    }

    curX = x;
  }

  void Minimizer::setG(std::vector < RealType > &g) {
    if (g.size() != ndim) {
      sprintf(painCave.errMsg, "Minimizer Error: dimesion of g and curG does not match\n");
      painCave.isFatal = 1;
      simError();
    }

    curG = g;
  }


  /**

  * In thoery, we need to find the minimum along the search direction
  * However, function evaluation is too expensive. 
  * At the very begining of the problem, we check the search direction and make sure
  * it is a descent direction
  * we will compare the energy of two end points,
  * if the right end point has lower energy, we just take it
  * @todo optimize this line search algorithm
  */

  int Minimizer::doLineSearch(std::vector<RealType> &direction,
			      RealType stepSize) {

    std::vector<RealType> xa;
    std::vector<RealType> xb;
    std::vector<RealType> xc;
    std::vector<RealType> ga;
    std::vector<RealType> gb;
    std::vector<RealType> gc;
    RealType fa;
    RealType fb;
    RealType fc;
    RealType a;
    RealType b;
    RealType c;
    int    status;
    RealType initSlope;
    RealType slopeA;
    RealType slopeB;
    RealType slopeC;
    bool   foundLower;
    int    iter;
    int    maxLSIter;
    RealType mu;
    RealType eta;
    RealType ftol;
    RealType lsTol;

    xa.resize(ndim);
    xb.resize(ndim);
    xc.resize(ndim);
    ga.resize(ndim);
    gb.resize(ndim);
    gc.resize(ndim);

    a = 0.0;

    fa = curF;

    xa = curX;

    ga = curG;

    c = a + stepSize;

    ftol = paramSet->getFTol();

    lsTol = paramSet->getLineSearchTol();

    //calculate the derivative at a = 0 

    slopeA = 0;

    for(size_t i = 0; i < ndim; i++) {
      slopeA += curG[i] * direction[i];
    }
    
    initSlope = slopeA;

    // if  going uphill, use negative gradient as searching direction 

    if (slopeA > 0) {

      for(size_t i = 0; i < ndim; i++) {
	direction[i] = -curG[i];
      }
        
      for(size_t i = 0; i < ndim; i++) {
	slopeA += curG[i] * direction[i];
      }
        
      initSlope = slopeA;
    }

    // Take a trial step

    for(size_t i = 0; i < ndim; i++) {
      xc[i] = curX[i] + direction[i]* c;
    }
    
    calcG(xc, gc, fc, status);

    if (status < 0) {
      if (bVerbose)
	std::cerr << "Function Evaluation Error" << std::endl;
    }

    //calculate the derivative at c

    slopeC = 0;

    for(size_t i = 0; i < ndim; i++) {
      slopeC += gc[i] * direction[i];
    }
    // found a lower point

    if (fc < fa) {
      curX = xc;

      curG = gc;

      curF = fc;

      return LS_SUCCEED;
    } else {
      if (slopeC > 0)
	stepSize *= 0.618034;
    }

    maxLSIter = paramSet->getLineSearchMaxIteration();

    iter = 0;

    do {

      // Select a new trial point.

      // If the derivatives at points a & c have different sign we use cubic interpolate    

      //if (slopeC > 0){     

      eta = 3 * (fa - fc) / (c - a) + slopeA + slopeC;

      mu = sqrt(eta * eta - slopeA * slopeC);

      b = a + (c - a)
	* (1 - (slopeC + mu - eta) / (slopeC - slopeA + 2 * mu));

      if (b < lsTol) {
	break;
      }

      //}

      // Take a trial step to this new point - new coords in xb 

      for(size_t i = 0; i < ndim; i++) {
	xb[i] = curX[i] + direction[i]* b;
      }
        
      //function evaluation

      calcG(xb, gb, fb, status);

      if (status < 0) {
	if (bVerbose)
	  std::cerr << "Function Evaluation Error" << std::endl;
      }

      //calculate the derivative at c

      slopeB = 0;

      for(size_t i = 0; i < ndim; i++) {
	slopeB += gb[i] * direction[i];
      }
        
      //Amijo Rule to stop the line search 

      if (fb <= curF +  initSlope * ftol * b) {
	curF = fb;

	curX = xb;

	curG = gb;

	return LS_SUCCEED;
      }

      if (slopeB < 0 && fb < fa) {

	//replace a by b

	fa = fb;

	a = b;

	slopeA = slopeB;

	// swap coord  a/b 

	std::swap(xa, xb);

	std::swap(ga, gb);
      } else {

	//replace c by b

	fc = fb;

	c = b;

	slopeC = slopeB;

	// swap coord  b/c 

	std::swap(gb, gc);

	std::swap(xb, xc);
      }

      iter++;
    } while ((fb > fa || fb > fc) && (iter < maxLSIter));

    if (fb < curF || iter >= maxLSIter) {

      //could not find a lower value, we might just go uphill.      

      return LS_ERROR;
    }

    //select the end point

    if (fa <= fc) {
      curX = xa;

      curG = ga;

      curF = fa;
    } else {
      curX = xc;

      curG = gc;

      curF = fc;
    }

    return LS_SUCCEED;
  }

  void Minimizer::minimize() {
    int convgStatus;
    int stepStatus;
    int maxIter;
    int writeFreq;
    int nextWriteIter;
    Snapshot* curSnapshot =info->getSnapshotManager()->getCurrentSnapshot();
    DumpWriter dumpWriter(info);     
    StatsBitSet mask;
    mask.set(Stats::TIME);
    mask.set(Stats::POTENTIAL_ENERGY);
    StatWriter statWriter(info->getStatFileName(), mask);

    init();

    writeFreq = paramSet->getWriteFreq();

    nextWriteIter = writeFreq;

    maxIter = paramSet->getMaxIteration();

    for(curIter = 1; curIter <= maxIter; curIter++) {
      stepStatus = step();

      //if (usingShake)
      //    preMove();

      if (stepStatus < 0) {
	saveResult();

	minStatus = MIN_LSERROR;

	std::cerr
	  << "Minimizer Error: line search error, please try a small stepsize"
	  << std::endl;

	return;
      }

      //save snapshot
      info->getSnapshotManager()->advance();
      //increase time
      curSnapshot->increaseTime(1);    
        
      if (curIter == nextWriteIter) {
	nextWriteIter += writeFreq;
	calcF();
	dumpWriter.writeDumpAndEor();
	statWriter.writeStat(curSnapshot->statData);
      }

      convgStatus = checkConvg();

      if (convgStatus > 0) {
	saveResult();

	minStatus = MIN_CONVERGE;

	return;
      }

      prepareStep();
    }

    if (bVerbose) {
      std::cout << "Minimizer Warning: " << minimizerName
		<< " algorithm did not converge within " << maxIter << " iteration"
		<< std::endl;
    }

    minStatus = MIN_MAXITER;

    saveResult();
  }


  RealType Minimizer::calcPotential() {
    forceMan->calcForces();

    Snapshot* curSnapshot = info->getSnapshotManager()->getCurrentSnapshot();
    RealType potential_local = curSnapshot->statData[Stats::LONG_RANGE_POTENTIAL] + 
      curSnapshot->statData[Stats::SHORT_RANGE_POTENTIAL] ;    
    RealType potential;

#ifdef IS_MPI
    MPI_Allreduce(&potential_local, &potential, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
#else
    potential = potential_local;
#endif

    //save total potential
    curSnapshot->statData[Stats::POTENTIAL_ENERGY] = potential;
    return potential;
  }

}
