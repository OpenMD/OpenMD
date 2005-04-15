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
 
#ifndef MINIMIZERS_OOPSEMINIMIZER_HPP
#define MINIMIZERS_OOPSEMINIMIZER_HPP

#include <iostream>

#include "io/DumpWriter.hpp"
#include "io/StatWriter.hpp"
#include "minimizers/MinimizerParameterSet.hpp"
#include "brains/ForceManager.hpp"

// base class of minimizer

namespace oopse {

  /** @todo need refactorying */
  const int MIN_LSERROR = -1;
  const int MIN_MAXITER = 0;
  const int MIN_CONVERGE = 1;

  const int CONVG_UNCONVG = 0;
  const int CONVG_FTOL = 1;
  const int CONVG_GTOL = 2;
  const int CONVG_ABSGTOL = 3;
  const int CONVG_STEPTOL = 4;

  const int LS_SUCCEED =1;
  const int LS_ERROR  = -1;

  /** @todo move to math module */
  double dotProduct(const std::vector<double>& v1, const std::vector<double>& v2);

  /**
   * @class Minimizer
   * base minimizer class
   */
  class Minimizer {
  public:

    Minimizer(SimInfo *rhs);

    virtual ~Minimizer();

    //
    virtual void init() {}

    //driver function of minimization method
    virtual void minimize();

    //
    virtual int step() = 0;

    //
    virtual void prepareStep() {};

    //line search algorithm, for the time being, we use back track algorithm
    virtual int doLineSearch(std::vector<double>& direction, double stepSize);

    virtual int checkConvg() = 0;

    //save the result when minimization method is done
    virtual void saveResult(){}

    //get the status of minimization
    int getMinStatus() {return minStatus;}

    // get the dimension of the model
    int getDim() {  return ndim;  }

    //get the name of minimizer method
    std::string getMinimizerName() {  return minimizerName;  }

    //return number of  current Iteration  
    int getCurIter() {  return curIter;  }

    // set the verbose mode of minimizer
    void setVerbose(bool verbose) {  bVerbose = verbose;}

    //get and set the coordinate
    std::vector<double> getX() {  return curX;  }
    void setX(std::vector<double>& x);

    //get and set the value of object function
    double getF() {  return curF;  }
    void setF(double f)  { curF = f;  }

    std::vector<double> getG() {  return curG;  }
    void setG(std::vector<double>& g);

    //get and set the gradient
    std::vector<double> getGrad() {  return curG;  }
    void setGrad(std::vector<double>& g) {  curG = g;  }

    //interal function to evaluate the energy and gradient in OOPSE
    void calcEnergyGradient(std::vector<double>& x,  std::vector<double>& grad, double&
			    energy, int& status);

    //calculate the value of object function
    virtual void calcF();
    virtual void calcF(std::vector<double>& x, double&f, int& status);

    //calculate the gradient
    virtual void calcG();
    virtual void calcG(std::vector<double>& x,  std::vector<double>& g, double& f, int& status);

    //calculate the hessian
    //virtual void calcH(int& status);
    //virtual void calcH(vector<double>& x,  std::vector<dobule>& g, SymMatrix& h, int& status);

    friend std::ostream& operator<<(std::ostream& os, const Minimizer& minimizer);

  protected:

    // transfrom cartesian and rotational coordinates into minimization coordinates 
    std::vector<double> getCoor();

    // transfrom minimization coordinates into cartesian and rotational coordinates  
    void setCoor(std::vector<double>& x);



    //constraint the bonds;
    int shakeR() { return 0;}

    //remove the force component along the bond direction
    int shakeF() { return 0;}

    double calcPotential();
        
    SimInfo* info;

    ForceManager* forceMan;
        
    //parameter set of minimization method
    MinimizerParameterSet* paramSet;

    //flag of turning on shake algorithm 
    bool usingShake;
        
    // dimension of the model
    int ndim;

    //name of the minimizer
    std::string minimizerName;

    // current iteration number
    int curIter;
    //status of minimization
    int minStatus;

    //flag of verbose mode
    bool bVerbose;

    //status of energy and gradient evaluation
    int egEvalStatus;

    //initial coordinates
    //vector<double> initX;

    //current value  of the function
    double curF;
        
    // current coordinates
    std::vector<double> curX;

    //gradient at curent coordinates
    std::vector<double> curG;

    //hessian at current coordinates
    //SymMatrix curH;

  private:

    //calculate the dimension od the model for minimization
    void calcDim();

  };

}
#endif

