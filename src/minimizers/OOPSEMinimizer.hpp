#ifndef _OOPSEMINIMIZER_H_
#define _OOPSEMINIMIZER_H_

#include <iostream>

#include "integrators/Integrator.hpp"
#include "minimizers/MinimizerParameterSet.hpp"


using namespace std;

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

// base class of minimizer

class OOPSEMinimizer : public RealIntegrator{
  public:

    OOPSEMinimizer(SimInfo *theInfo, ForceFields* the_ff, 
                                     MinimizerParameterSet* param);

    virtual ~OOPSEMinimizer();

    //
    virtual void init() {}
    
    //driver function of minimization method
    virtual void minimize();

    //
    virtual int step() = 0;

    //
    virtual void prepareStep() {};

    //line search algorithm, for the time being, we use back track algorithm
    virtual int doLineSearch(vector<double>& direction, double stepSize);

    virtual int checkConvg() = 0;

    //print detail information of the minimization method
    virtual void printMinimizerInfo();
    
    //save the result when minimization method is done
    virtual void saveResult(){}

    //write out the trajectory
    virtual void writeOut(vector<double>&x, double curIter);


    //get the status of minimization
    int getMinStatus() {return minStatus;}

    // get the dimension of the model
    int getDim() {  return ndim;  }

    //get the name of minimizer method
    string getMinimizerName() {  return minimizerName;  }

    //return number of  current Iteration  
    int getCurIter() {  return curIter;  }

    // set the verbose mode of minimizer
    void setVerbose(bool verbose) {  bVerbose = verbose;}

    //get and set the coordinate
    vector<double> getX() {  return curX;  }
    void setX(vector<double>& x);

     //get and set the value of object function
    double getF() {  return curF;  }
    void setF(double f)  { curF = f;  }

    vector<double> getG() {  return curG;  }
    void setG(vector<double>& g);
    
    //get and set the gradient
    vector<double> getGrad() {  return curG;  }
    void setGrad(vector<double>& g) {  curG = g;  }

    //interal function to evaluate the energy and gradient in OOPSE
    void calcEnergyGradient(vector<double>& x, vector<double>& grad, double&
                                                 energy, int& status);

    //calculate the value of object function
    virtual void calcF();
    virtual void calcF(vector<double>& x, double&f, int& status);

    //calculate the gradient
    virtual void calcG();
    virtual void calcG(vector<double>& x, vector<double>& g, double& f, int& status);

    //calculate the hessian
    //virtual void calcH(int& status);
    //virtual void calcH(vector<double>& x, vector<dobule>& g, SymMatrix& h, int& status);


  protected:

    // transfrom cartesian and rotational coordinates into minimization coordinates 
    vector<double> getCoor();
    
    // transfrom minimization coordinates into cartesian and rotational coordinates  
    void setCoor(vector<double>& x);

    //flag of turning on shake algorithm 
    bool bShake;

    //constraint the bonds;
    int shakeR();
    
    //remove the force component along the bond direction
    int shakeF();
   
    // dimension of the model
    int ndim;

    //name of the minimizer
    string minimizerName;

    // current iteration number
    int curIter;
    //status of minimization
    int minStatus;

    //flag of verbose mode
    bool bVerbose;


    //parameter set of minimization method
    MinimizerParameterSet* paramSet;

    //status of energy and gradient evaluation
    int egEvalStatus;

    //initial coordinates
    //vector<double> initX;
    
    //current value  of the function
    double curF;
    // current coordinates
    vector<double> curX;
    //gradient at curent coordinates
    vector<double> curG;

    //hessian at current coordinates
    //SymMatrix curH;

  private:

    //calculate the dimension od the model for minimization
    void calcDim();

};

// steepest descent minimizer
class SDMinimizer : public OOPSEMinimizer{
  public:
    SDMinimizer(SimInfo *theInfo, ForceFields* the_ff, 
                             MinimizerParameterSet* param); 

    virtual void init();    
    virtual int step();
    virtual void prepareStep();
    virtual int checkConvg();
  protected:

    vector<double> direction;
    double prevF;
    double stepSize;

};

// conjugate gradient famlily minimzier
class CGFamilyMinmizer : public OOPSEMinimizer{
  public:
    CGFamilyMinmizer(SimInfo *theInfo, ForceFields* the_ff, 
                                          MinimizerParameterSet* param);

    //check the convergence
    virtual int checkConvg();

    
  protected:

    vector<double> direction;
    vector<double> prevX;
    vector<double> prevG;
    double prevF;
    double stepSize;
};

//Polak-Ribiere  Conjugate Gradient Method 
class PRCGMinimizer : public CGFamilyMinmizer{
  public:
      PRCGMinimizer(SimInfo *theInfo, ForceFields* the_ff,  
                                    MinimizerParameterSet* param)
                                   :CGFamilyMinmizer(theInfo, the_ff, param) {}
    virtual int step();
    virtual void init();
    virtual void prepareStep();
};

#endif
