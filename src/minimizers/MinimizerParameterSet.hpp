#ifndef _MININIZERPARAMETERSET_H_
#define _MININIZERPARAMETERSET_H_

const double DEFAULTTOLERANCE = 1.0e-8;

// base class of minimizer's parameter set
class MinimizerParameterSet{
  public:

    MinimizerParameterSet() {setDefaultParameter();}
    MinimizerParameterSet(MinimizerParameterSet& param) {*this = param;}
    
    void operator= (MinimizerParameterSet& param) {

      maxIteration = param.maxIteration;
      stepSize = param.stepSize;
      stepTol = param.stepTol;
      fTol = param.fTol;
      gTol = param.gTol;

      writeFrq = param.writeFrq;

      lsMaxIteration = param.lsMaxIteration;
      lsTol = param.lsTol;

}
    
    virtual void setDefaultParameter(){
      maxIteration = 200;
      stepSize = 0.01;
      stepTol = DEFAULTTOLERANCE;
      fTol = DEFAULTTOLERANCE;
      gTol = DEFAULTTOLERANCE;

      writeFrq = maxIteration;
      
      lsMaxIteration = 50;
      lsTol = DEFAULTTOLERANCE;
    } 
      
    void setStepTol(double tol) { stepTol = tol;}
    double getStepTol() { return stepTol;}

    void setStepSize(double size) {  stepSize = size;}
    double getStepSize() { return stepSize;}

    void setMaxIteration(int iter) { maxIteration = iter;}
    int getMaxIteration() {return maxIteration;}

    void setFTol(double tol) {fTol = tol;} 
    double getFTol() {return fTol;}

    void setGTol(double tol) {gTol = tol;}
    double getGTol() {return gTol;}

    void setLineSearchTol(double tol) {lsTol = tol;}
    double getLineSearchTol() {return lsTol;}

    void setLineSearchMaxIteration(int iter) {lsMaxIteration = iter;}
    int getLineSearchMaxIteration() {return lsMaxIteration;}

    void setWriteFrq(int frq) {writeFrq = frq;}
    int getWriteFrq() {return writeFrq;}

  protected:    

    int maxIteration;
    double stepTol;
    double fTol;
    double gTol;
    double stepSize;

    int lsMaxIteration;
    double lsTol;

    int writeFrq;
    //int resetFrq;
/*
    // Absolute tolerance
    vector<double> absTol;

    // Relative tolerance
    vector<double> relTol;

    // Total specified tolerance at convergence 
    vector<double> totTol;

    // Tolerance achieved at convergence.
    vector<double> achTol;
*/
};


#endif
