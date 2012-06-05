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
 
#ifndef MINIMIZERS_MINIMIZER_HPP
#define MINIMIZERS_MINIMIZER_HPP

#include <iostream>

#include "io/DumpWriter.hpp"
#include "io/StatWriter.hpp"
#include "minimizers/MinimizerParameterSet.hpp"
#include "brains/ForceManager.hpp"

// base class of minimizer

namespace OpenMD {
  
  /** minimizer stop codes */
  enum{MIN_MAXITER,
       MIN_MAXEVAL,
       MIN_ETOL,
       MIN_FTOL,
       MIN_DOWNHILL,
       MIN_ZEROALPHA,
       MIN_ZEROFORCE,
       MIN_ZEROQUAD};
  
  typedef int (*FnPtr)(std::vector<RealType> &, std::vector<RealType> &, RealType &);    

  /**
   * @class Minimizer
   * base minimizer class
   */
  class Minimizer {
  public:
    Minimizer(SimInfo *rhs);
    virtual ~Minimizer();
    virtual void init() {}
    //driver function of minimization method
    virtual void minimize();
    virtual int step() = 0;
    virtual void prepareStep() {};

    //line search algorithm, for the time being, we use a back track
    //algorithm
    virtual int doLineSearch(std::vector<RealType>& direction, RealType stepSize);
    virtual int checkConvergence() = 0;

    //save the result when minimization method is done
    virtual void saveResult(){}
    
    //get the status of minimization
    int getMinimizerStatus() {return minStatus;}

    // get the dimension of the model
    int getDim() {  return ndim;  }

    //get the name of minimizer method
    std::string getMinimizerName() {  return minimizerName;  }

    //return number of the current Iteration  
    int getCurrentIteration() {  return curIter;  }

    // set the verbose mode of minimizer
    void setVerbose(bool verbose) {  bVerbose = verbose;}

    //get and set the coordinates
    std::vector<RealType> getX() {  return curX;  }
    void setX(std::vector<RealType>& x);

    //get and set the value of object function
    RealType getF() {  return curF;  }
    void setF(RealType f)  { curF = f;  }

    std::vector<RealType> getG() {  return curG;  }
    void setG(std::vector<RealType>& g);

    //get and set the gradient
    std::vector<RealType> getGrad() {  return curG;  }
    void setGrad(std::vector<RealType>& g) {  curG = g;  }

    void setGradientFunction(FnPtr efunc) { calcEnergyGradient = efunc; }

    //calculate the value of object function
    virtual void calcF();
    virtual void calcF(std::vector<RealType>& x, RealType&f, int& status);

    //calculate the gradient
    virtual void calcG();
    virtual void calcG(std::vector<RealType>& x,  std::vector<RealType>& g, RealType& f, int& status);

    //calculate the hessian
    //virtual void calcH(int& status);
    //virtual void calcH(vector<RealType>& x,  std::vector<dobule>& g, SymMatrix& h, int& status);

    friend std::ostream& operator<<(std::ostream& os, const Minimizer& minimizer);

  protected:

    typedef int (*FnPtr)(std::vector<RealType> &, std::vector<RealType> &, RealType &);    
    FnPtr calcEnergyGradient;

    // transfrom cartesian and rotational coordinates into minimization coordinates 
    std::vector<RealType> getCoor();

    // transfrom minimization coordinates into cartesian and rotational coordinates  
    void setCoor(std::vector<RealType>& x);

    //constraint the bonds;
    int shakeR() { return 0;}

    //remove the force component along the bond direction
    int shakeF() { return 0;}
    RealType calcPotential();
        
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
    //vector<RealType> initX;

    //current value  of the function
    RealType curF;
        
    // current coordinates
    std::vector<RealType> curX;

    //gradient at curent coordinates
    std::vector<RealType> curG;

    //hessian at current coordinates
    //SymMatrix curH;

  private:
    //calculate the dimension of the model for minimization
    void calcDim();
  };

}
#endif

