#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include <string>
#include <vector>
#include "primitives/Atom.hpp"
#include "primitives/StuntDouble.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/SRI.hpp"
#include "primitives/AbstractClasses.hpp"
#include "brains/SimInfo.hpp"
#include "UseTheForce/ForceFields.hpp"
#include "brains/Thermo.hpp"
#include "io/ReadWrite.hpp"
#include "io/ZConsWriter.hpp"
#include "restraints/Restraints.hpp"

using namespace std;
const double kB = 8.31451e-7;// boltzmann constant amu*Ang^2*fs^-2/K
const double eConvert = 4.184e-4; // converts kcal/mol -> amu*A^2/fs^2
const double p_convert = 1.63882576e8; //converts amu*fs^-2*Ang^-1 -> atm
const int maxIteration = 300;
const double tol = 1.0e-6;

class VelVerletConsFramework;
template<typename T = BaseIntegrator> class Integrator : public T {

public:
  Integrator( SimInfo *theInfo, ForceFields* the_ff );
  virtual ~Integrator();
  void integrate( void );
  virtual double  getConservedQuantity(void);
  virtual string getAdditionalParameters(void);

protected:

  virtual void integrateStep( int calcPot, int calcStress );
  virtual void preMove( void );
  virtual void moveA( void );
  virtual void moveB( void );
  virtual void constrainA( void );
  virtual void constrainB( void );
  virtual int  readyCheck( void ) { return 1; }

  virtual void resetIntegrator( void ) { }

  virtual void calcForce( int calcPot, int calcStress );
  virtual void thermalize();

  virtual bool stopIntegrator() {return false;}

  virtual void rotationPropagation( StuntDouble* sd, double ji[3] );

  void checkConstraints( void );
  void rotate( int axes1, int axes2, double angle, double j[3],
         double A[3][3] );

  ForceFields* myFF;

  SimInfo *info; // all the info we'll ever need
  vector<StuntDouble*> integrableObjects;
  int nAtoms;  /* the number of atoms */
  int oldAtoms;
  Atom **atoms; /* array of atom pointers */
  Molecule* molecules;
  int nMols;


  int isConstrained; // boolean to know whether the systems contains constraints.
  int nConstrained;  // counter for number of constraints
  int *constrainedA; // the i of a constraint pair
  int *constrainedB; // the j of a constraint pair
  double *constrainedDsqr; // the square of the constraint distance

  int* moving; // tells whether we are moving atom i
  int* moved;  // tells whether we have moved atom i
  double* oldPos; // pre constrained positions

  short isFirst; /*boolean for the first time integrate is called */

  double dt;
  double dt2;

  Thermo *tStats;
  StatWriter*  statOut;
  DumpWriter*  dumpOut;

};

typedef Integrator<BaseIntegrator> RealIntegrator;

// ansi instantiation
// template class Integrator<BaseIntegrator>;


template<typename T> class NVE : public T {

public:
  NVE ( SimInfo *theInfo, ForceFields* the_ff ):
    T( theInfo, the_ff ){}
  virtual ~NVE(){}
};


template<typename T> class NVT : public T {

public:

  NVT ( SimInfo *theInfo, ForceFields* the_ff);
  virtual ~NVT();

  void setTauThermostat(double tt) {tauThermostat = tt; have_tau_thermostat=1;}
  void setTargetTemp(double tt) {targetTemp = tt; have_target_temp = 1;}
  void setChiTolerance(double tol) {chiTolerance = tol;}
  virtual double  getConservedQuantity(void);
  virtual string getAdditionalParameters(void);

protected:

  virtual void moveA( void );
  virtual void moveB( void );

  virtual int readyCheck();

  virtual void resetIntegrator( void );

  // chi is a propagated degree of freedom.

  double chi;

  //integral of chi(t)dt
  double integralOfChidt;

  // targetTemp must be set.  tauThermostat must also be set;

  double targetTemp;
  double tauThermostat;

  short int have_tau_thermostat, have_target_temp;

  double *oldVel;
  double *oldJi;

  double chiTolerance;
  short int have_chi_tolerance;

};



template<typename T> class NPT : public T{

public:

  NPT ( SimInfo *theInfo, ForceFields* the_ff);
  virtual ~NPT();

  virtual void integrateStep( int calcPot, int calcStress ){
    calcStress = 1;
    T::integrateStep( calcPot, calcStress );
  }

  virtual double getConservedQuantity(void) = 0;
  virtual string getAdditionalParameters(void) = 0;
  
  double myTauThermo( void ) { return tauThermostat; }
  double myTauBaro( void ) { return tauBarostat; }

  void setTauThermostat(double tt) {tauThermostat = tt; have_tau_thermostat=1;}
  void setTauBarostat(double tb) {tauBarostat = tb; have_tau_barostat=1;}
  void setTargetTemp(double tt) {targetTemp = tt; have_target_temp = 1;}
  void setTargetPressure(double tp) {targetPressure = tp; have_target_pressure = 1;}
  void setChiTolerance(double tol) {chiTolerance = tol; have_chi_tolerance = 1;}
  void setPosIterTolerance(double tol) {posIterTolerance = tol; have_pos_iter_tolerance = 1;}
  void setEtaTolerance(double tol) {etaTolerance = tol; have_eta_tolerance = 1;}

protected:

  virtual void  moveA( void );
  virtual void moveB( void );

  virtual int readyCheck();

  virtual void resetIntegrator( void );

  virtual void getVelScaleA( double sc[3], double vel[3] ) = 0;
  virtual void getVelScaleB( double sc[3], int index ) = 0;
  virtual void getPosScale(double pos[3], double COM[3],
			   int index, double sc[3]) = 0;

  virtual void calcVelScale( void ) = 0;

  virtual bool chiConverged( void );
  virtual bool etaConverged( void ) = 0;

  virtual void evolveChiA( void );
  virtual void evolveEtaA( void ) = 0;
  virtual void evolveChiB( void );
  virtual void evolveEtaB( void ) = 0;

  virtual void scaleSimBox( void ) = 0;

  void accIntegralOfChidt(void) { integralOfChidt += dt * chi;}

  // chi and eta are the propagated degrees of freedom

  double oldChi;
  double prevChi;
  double chi;
  double NkBT;
  double fkBT;

  double tt2, tb2;
  double instaTemp, instaPress, instaVol;
  double press[3][3];

  int Nparticles;

  double integralOfChidt;

  // targetTemp, targetPressure, and tauBarostat must be set.
  // One of qmass or tauThermostat must be set;

  double targetTemp;
  double targetPressure;
  double tauThermostat;
  double tauBarostat;

  short int have_tau_thermostat, have_tau_barostat, have_target_temp;
  short int have_target_pressure;

  double *oldPos;
  double *oldVel;
  double *oldJi;

  double chiTolerance;
  short int have_chi_tolerance;
  double posIterTolerance;
  short int have_pos_iter_tolerance;
  double etaTolerance;
  short int have_eta_tolerance;

};

template<typename T> class NPTi : public T{

public:
  NPTi( SimInfo *theInfo, ForceFields* the_ff);
  ~NPTi();

  virtual double getConservedQuantity(void);
  virtual void resetIntegrator(void);
  virtual string getAdditionalParameters(void);
protected:



  virtual void evolveEtaA(void);
  virtual void evolveEtaB(void);

  virtual bool etaConverged( void );

  virtual void scaleSimBox( void );

  virtual void getVelScaleA( double sc[3], double vel[3] );
  virtual void getVelScaleB( double sc[3], int index );
  virtual void getPosScale(double pos[3], double COM[3],
			   int index, double sc[3]);

  virtual void calcVelScale( void );

  double eta, oldEta, prevEta;
  double vScale;
};

template<typename T> class NPTf : public T{

public:

  NPTf ( SimInfo *theInfo, ForceFields* the_ff);
  virtual ~NPTf();

  virtual double getConservedQuantity(void);
  virtual string getAdditionalParameters(void);
  virtual void resetIntegrator(void);

protected:

  virtual void evolveEtaA(void);
  virtual void evolveEtaB(void);

  virtual bool etaConverged( void );

  virtual void scaleSimBox( void );

  virtual void getVelScaleA( double sc[3], double vel[3] );
  virtual void getVelScaleB( double sc[3], int index );
  virtual void getPosScale(double pos[3], double COM[3],
			   int index, double sc[3]);

  virtual void calcVelScale( void );

  double eta[3][3];
  double oldEta[3][3];
  double prevEta[3][3];
  double vScale[3][3];
};

template<typename T> class NPTxyz : public T{

public:

  NPTxyz ( SimInfo *theInfo, ForceFields* the_ff);
  virtual ~NPTxyz();

  virtual double getConservedQuantity(void);
  virtual string getAdditionalParameters(void);
  virtual void resetIntegrator(void);

protected:

  virtual void evolveEtaA(void);
  virtual void evolveEtaB(void);

  virtual bool etaConverged( void );

  virtual void scaleSimBox( void );

  virtual void getVelScaleA( double sc[3], double vel[3] );
  virtual void getVelScaleB( double sc[3], int index );
  virtual void getPosScale(double pos[3], double COM[3],
			   int index, double sc[3]);

  virtual void calcVelScale( void );

  double eta[3][3];
  double oldEta[3][3];
  double prevEta[3][3];
  double vScale[3][3];
};


template<typename T> class ZConstraint : public T {

  public:
  class ForceSubtractionPolicy{
    public:
      ForceSubtractionPolicy(ZConstraint<T>* integrator) {zconsIntegrator = integrator;}

      virtual void update() = 0;
      virtual double getZFOfFixedZMols(Molecule* mol, Atom* atom, double totalForce) = 0;
      virtual double getZFOfMovingMols(Atom* atom, double totalForce) = 0;
      virtual double getHFOfFixedZMols(Molecule* mol, Atom* atom, double totalForce) = 0;
      virtual double getHFOfUnconsMols(Atom* atom, double totalForce) = 0;

   protected:
     ZConstraint<T>* zconsIntegrator;
  };

  class PolicyByNumber : public ForceSubtractionPolicy{

    public:
      PolicyByNumber(ZConstraint<T>* integrator) :ForceSubtractionPolicy(integrator) {}
      virtual void update();
      virtual double getZFOfFixedZMols(Molecule* mol, Atom* atom, double totalForce) ;
      virtual double getZFOfMovingMols(Atom* atom, double totalForce) ;
      virtual double getHFOfFixedZMols(Molecule* mol, Atom* atom, double totalForce);
      virtual double getHFOfUnconsMols(Atom* atom, double totalForce);

    private:
      int totNumOfMovingAtoms;
  };

  class PolicyByMass : public ForceSubtractionPolicy{

    public:
      PolicyByMass(ZConstraint<T>* integrator) :ForceSubtractionPolicy(integrator) {}

      virtual void update();
      virtual double getZFOfFixedZMols(Molecule* mol, Atom* atom, double totalForce) ;
      virtual double getZFOfMovingMols(Atom* atom, double totalForce) ;
      virtual double getHFOfFixedZMols(Molecule* mol, Atom* atom, double totalForce);
      virtual double getHFOfUnconsMols(Atom* atom, double totalForce);

   private:
     double totMassOfMovingAtoms;
  };

public:

  ZConstraint( SimInfo *theInfo, ForceFields* the_ff);
  ~ZConstraint();

  void setZConsTime(double time)                  {this->zconsTime = time;}
  void getZConsTime()                             {return zconsTime;}

  void setIndexOfAllZConsMols(vector<int> index) {indexOfAllZConsMols = index;}
  void getIndexOfAllZConsMols()                  {return indexOfAllZConsMols;}

  void setZConsOutput(const char * fileName)          {zconsOutput = fileName;}
  string getZConsOutput()                         {return zconsOutput;}

  virtual void integrate();


#ifdef IS_MPI
  virtual void update();                      //which is called to indicate the molecules' migration
#endif

  enum ZConsState {zcsMoving, zcsFixed};

  vector<Molecule*> zconsMols;              //z-constraint molecules array
  vector<ZConsState> states;                 //state of z-constraint molecules



  int totNumOfUnconsAtoms;              //total number of uncontraint atoms
  double totalMassOfUncons;                //total mas of unconstraint molecules


protected:



  virtual void calcForce( int calcPot, int calcStress );
  virtual void thermalize(void);

  void zeroOutVel();
  void doZconstraintForce();
  void doHarmonic(vector<double>& resPos);
  bool checkZConsState();

  bool haveFixedZMols();
  bool haveMovingZMols();

  double calcZSys();

  int isZConstraintMol(Molecule* mol);


  double zconsTime;                              //sample time
  double zconsTol;                                 //tolerance of z-contratint
  double zForceConst;                           //base force constant term
                                                          //which is estimate by OOPSE


  vector<double> massOfZConsMols;       //mass of z-constraint molecule
  vector<double> kz;                              //force constant array

  vector<double> zPos;                          //


  vector<Molecule*> unconsMols;           //unconstraint molecules array
  vector<double> massOfUnconsMols;    //mass array of unconstraint molecules


  vector<ZConsParaItem>* parameters; //

  vector<int> indexOfAllZConsMols;     //index of All Z-Constraint Molecuels

  vector<int> indexOfZConsMols;                   //index of local Z-Constraint Molecules
  vector<double> fz;
  vector<double> curZPos;

  bool usingSMD;
  vector<double> prevCantPos;
  vector<double> cantPos;
  vector<double> cantVel;

  double zconsFixTime;  
  double zconsGap;
  bool hasZConsGap;
  vector<double> endFixTime;
  
  int whichDirection;                           //constraint direction

private:

  string zconsOutput;                         //filename of zconstraint output
  ZConsWriter* fzOut;                         //z-constraint writer

  double curZconsTime;

  double calcMovingMolsCOMVel();
  double calcSysCOMVel();
  double calcTotalForce();
  void updateZPos();
  void updateCantPos();
  
  ForceSubtractionPolicy* forcePolicy; //force subtraction policy
  friend class ForceSubtractionPolicy;

};


//Sympletic quaternion Scheme Integrator
//Reference:
// T.F. Miller, M. Eleftheriou, P. Pattnaik, A. Ndirango, D. Newns and G.J. Martyna
//Symplectic quaternion Scheme for biophysical molecular dynamics
//116(20), 8649, J. Chem. Phys. (2002)
template<typename T> class SQSIntegrator : public T{
  public:
    virtual void moveA();
    virtual void moveB();
  protected:
    void freeRotor();
    void rotate(int k, double dt);
    
};
#endif
