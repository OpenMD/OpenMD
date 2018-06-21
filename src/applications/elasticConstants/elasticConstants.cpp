/*
 * Copyright (c) 2006 The University of Notre Dame. All Rights Reserved.
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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <fstream>

#include "elasticConstantsCmd.hpp"
#include "brains/Register.hpp"
#include "brains/SimInfo.hpp"
#include "brains/SimCreator.hpp"
#include "brains/Thermo.hpp"
#include "brains/ForceManager.hpp"
#include "brains/Velocitizer.hpp"
#include "constraints/Shake.hpp"
#include "flucq/FluctuatingChargeConstraints.hpp"
#include "flucq/FluctuatingChargeDamped.hpp"
#include "io/DumpReader.hpp"
#include "utils/StringUtils.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/LU.hpp"
#include "math/DynamicRectMatrix.hpp"
// #include "optimization/OptimizationFactory.hpp"
// #include "optimization/Method.hpp"
// #include "optimization/Constraint.hpp"
// #include "optimization/Problem.hpp"
// #include "optimization/BoxObjectiveFunction.hpp"

using namespace OpenMD;
#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>
#include <iomanip>


template<typename Real>
class Vector6 : public Vector<Real, 6>{
public:
  typedef Real ElemType;
  typedef Real* ElemPoinerType;
  Vector6() : Vector<Real, 6>(){}
  
  /** Constructs and initializes a Vector6d from individual coordinates */
  inline Vector6( RealType v0, RealType v1, RealType v2, RealType v3,
                  RealType v4, RealType v5) {
    this->data_[0] = v0;
    this->data_[1] = v1;
    this->data_[2] = v2;
    this->data_[3] = v3;
    this->data_[4] = v4;
    this->data_[5] = v5;
  }
  /** Constructs and initializes from an array*/
  inline Vector6(Real* array) : Vector<Real, 6>(array) {}
  
  inline Vector6(const Vector<Real, 6>& v) : Vector<Real, 6>(v) {}
  
  inline Vector6<Real>& operator = (const Vector<Real, 6>& v) {
    if (this == &v) { return *this; }
    Vector<Real, 6>::operator=(v);
    return *this;
  }

};
  
  
typedef Vector6<RealType> Vector6d;

RealType slope(const std::vector<RealType>& x, const std::vector<RealType>& y) {  
  const size_t n    = x.size();  
  const RealType s_x  = std::accumulate(x.begin(), x.end(), 0.0);
  const RealType s_y  = std::accumulate(y.begin(), y.end(), 0.0);
  const RealType s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
  const RealType s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
  const RealType a    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  return a;
}

void quadraticFit(const std::vector<RealType>& x,
                  const std::vector<RealType>& y,
                  RealType& a, RealType& b, RealType& c) {
  
  RealType s00 = RealType(x.size());
  RealType s10(0.0), s20(0.0), s30(0.0), s40(0.0);
  RealType s01(0.0), s11(0.0), s21(0.0);
  
  for (size_t i = 0; i < x.size(); i++) {
    // std::cerr << x[i] << "\t" << y[i] << "\n";
    s10 += x[i];
    s20 += pow(x[i], 2);
    s30 += pow(x[i], 3);
    s40 += pow(x[i], 4);
    s01 += y[i];
    s11 += x[i]*y[i];
    s21 += pow(x[i],2) * y[i];
  }
  // std::cerr << "&\n";

  RealType D = (s40 * (s20 * s00 - s10 * s10) - 
                s30 * (s30 * s00 - s10 * s20) + 
                s20 * (s30 * s10 - s20 * s20));
  
  a = (s21*(s20 * s00 - s10 * s10) - 
       s11*(s30 * s00 - s10 * s20) + 
       s01*(s30 * s10 - s20 * s20)) / D;
  
  b = (s40*(s11 * s00 - s01 * s10) - 
       s30*(s21 * s00 - s01 * s20) + 
       s20*(s21 * s10 - s11 * s20)) / D;
  
  c = (s40*(s20 * s01 - s10 * s11) - 
       s30*(s30 * s01 - s10 * s21) + 
       s20*(s30 * s11 - s20 * s21)) / D;

  return;
}

void writeMatrix(DynamicRectMatrix<RealType> M, std::string title,
                 std::string units) {

  std::cout << left << title << " (" << units << "):" << std::endl;
  std::cout << std::endl;
    
  std::cout << " ";
  std::cout << right << setw(12) << M(0,0) << " ";
  std::cout << right << setw(12) << M(0,1) << " ";
  std::cout << right << setw(12) << M(0,2) << " ";
  std::cout << right << setw(12) << M(0,3) << " ";
  std::cout << right << setw(12) << M(0,4) << " ";                
  std::cout << right << setw(12) << M(0,5);
  std::cout << std::endl;

  std::cout << " ";
  std::cout << right << setw(12) << M(1,0) << " ";
  std::cout << right << setw(12) << M(1,1) << " ";
  std::cout << right << setw(12) << M(1,2) << " ";
  std::cout << right << setw(12) << M(1,3) << " ";
  std::cout << right << setw(12) << M(1,4) << " ";                
  std::cout << right << setw(12) << M(1,5);
  std::cout << std::endl;

  std::cout << " ";
  std::cout << right << setw(12) << M(2,0) << " ";
  std::cout << right << setw(12) << M(2,1) << " ";
  std::cout << right << setw(12) << M(2,2) << " ";
  std::cout << right << setw(12) << M(2,3) << " ";
  std::cout << right << setw(12) << M(2,4) << " ";
  std::cout << right << setw(12) << M(2,5);
  std::cout << std::endl;

  std::cout << " ";
  std::cout << right << setw(12) << M(3,0) << " ";
  std::cout << right << setw(12) << M(3,1) << " ";
  std::cout << right << setw(12) << M(3,2) << " ";
  std::cout << right << setw(12) << M(3,3) << " ";
  std::cout << right << setw(12) << M(3,4) << " ";                
  std::cout << right << setw(12) << M(3,5);
  std::cout << std::endl;

  std::cout << " ";
  std::cout << right << setw(12) << M(4,0) << " ";
  std::cout << right << setw(12) << M(4,1) << " ";
  std::cout << right << setw(12) << M(4,2) << " ";
  std::cout << right << setw(12) << M(4,3) << " ";
  std::cout << right << setw(12) << M(4,4) << " ";                
  std::cout << right << setw(12) << M(4,5);
  std::cout << std::endl;
    
  std::cout << " ";
  std::cout << right << setw(12) << M(5,0) << " ";
  std::cout << right << setw(12) << M(5,1) << " ";
  std::cout << right << setw(12) << M(5,2) << " ";
  std::cout << right << setw(12) << M(5,3) << " ";
  std::cout << right << setw(12) << M(5,4) << " ";
  std::cout << right << setw(12) << M(5,5);
  std::cout << std::endl;
  std::cout << std::endl;    
}

void writeMaterialProperties(DynamicRectMatrix<RealType> C,
                             DynamicRectMatrix<RealType> S) {
  
  RealType C11 = C(0,0);
  RealType C22 = C(1,1);
  RealType C33 = C(2,2);
  RealType C12 = C(0,1);
  RealType C23 = C(1,2);
  RealType C31 = C(2,0);
  RealType C44 = C(3,3);
  RealType C55 = C(4,4);
  RealType C66 = C(5,5);
  
  RealType S11 = S(0,0);
  RealType S22 = S(1,1);
  RealType S33 = S(2,2);
  RealType S12 = S(0,1);
  RealType S23 = S(1,2);
  RealType S31 = S(2,0);
  RealType S44 = S(3,3);
  RealType S55 = S(4,4);
  RealType S66 = S(5,5);


  // Material Properties defined in:
  //
  // "Charting the complete elastic properties of inorganic
  // crystalline compounds," Maarten de Jong, Wei Chen, Thomas
  // Angsten, Anubhav Jain, Randy Notestine, Anthony Gamst, Marcel
  // Sluiter, Chaitanya Krishna Ande, Sybrand van der Zwaag, Jose J
  // Plata, Cormac Toher, Stefano Curtarolo, Gerbrand Ceder, Kristin
  // A. Persson & Mark Asta,
  //
  // Scientific Data volume 2, Article number: 150009 (2015)
  // doi:10.1038/sdata.2015.9
  //
  // And in
  //
  // "ELATE: an open-source online application for analysis and
  // visualization of elastic tensors," Romain Gaillac, Pluton
  // Pullumbi and François-Xavier Coudert,
  //
  // J. Phys.: Condens. Matter 28 275201 (2016)
  // doi:10.1088/0953-8984/28/27/275201
  
  // Bulk modulus (Voigt average):
  RealType Kv = ((C11+C22+C33) + 2.0*(C12+C23+C31)) / 9.0;
  // Bulk modulus (Reuss average):
  RealType Kr = 1.0 / ((S11+S22+S33) + 2.0*(S12+S23+S31));
  // Shear modulus (Voigt average):
  RealType Gv = ((C11+C22+C33) - (C12+C23+C31) + 3.0*(C44+C55+C66)) / 15.0;
  // Shear modulus (Reuss average):
  RealType Gr = 15.0 / (4.0*(S11+S22+S33) - 4.0*(S12+S23+S31) + 3.0*(S44+S55+S66));
  // Bulk modulus (Hill average):
  RealType Kh = (Kv + Kr) / 2.0;
  // Shear modulus (Hill average):
  RealType Gh = (Gv + Gr) / 2.0;
  // Universal elastic anisotropy
  RealType Au = 5.0 * (Gv/Gr) + (Kv/Kr) - 6.0;
  // Isotropic Poisson ratio
  RealType muv = (1.0 - 3.0*Gv/(3.0*Kv + Gv)) / 2.0;
  RealType mur = (1.0 - 3.0*Gr/(3.0*Kr + Gr)) / 2.0;
  RealType muh = (1.0 - 3.0*Gh/(3.0*Kh + Gh)) / 2.0;
  // Isotropic Young's modulus
  RealType Ev = 1.0/(1.0/(3.0*Gv) + 1.0/(9.0*Kv));
  RealType Er = 1.0/(1.0/(3.0*Gr) + 1.0/(9.0*Kr));
  RealType Eh = 1.0/(1.0/(3.0*Gh) + 1.0/(9.0*Kh));

  std::cout << "                              "
            << setw(12) << "Voigt" << " "
            << setw(12) << "Reuss" << " "
            << setw(12) << "Hill\n";
  std::cout << "Bulk modulus:                 "
            << setw(12) << Kv << " "
            << setw(12) << Kr << " "
            << setw(12) << Kh << " (GPa)\n";

  std::cout << "Shear modulus:                "
            << setw(12) << Gv << " "
            << setw(12) << Gr << " "
            << setw(12) << Gh << " (GPa)\n";

  std::cout << "Young\'s modulus (isotropic):  "
            << setw(12) << Ev << " "
            << setw(12) << Er << " "
            << setw(12) << Eh << " (GPa)\n";

  std::cout << "Poisson\'s Ratio:              "
            << setw(12) << muv << " "
            << setw(12) << mur << " "
            << setw(12) << muh << "\n";
  
  std::cout << "Universal elastic Anisotropy: " << setw(12) << Au << "\n";

    
  // Assume a cubic crystal, and use symmetries:

  C11 = (C(0,0) + C(1,1) + C(2,2)) / 3.0;
  C12 = (C(0,1) + C(0,2) + C(1,2)) / 3.0;
  C44 = (C(3,3) + C(4,4) + C(5,5)) / 3.0;
  S11 = (S(0,0) + S(1,1) + S(2,2)) / 3.0;
  S12 = (S(0,1) + S(0,2) + S(1,2)) / 3.0;
  S44 = (S(3,3) + S(4,4) + S(5,5)) / 3.0;
      
  // Anisotropy factor
  RealType A1 = 2.0*C44 / (C11 - C12);
  RealType A2 = 2.0*(S11 - S12) / S44;
  
  // std::cout << "Anisotropy factor = " << A1 << " " << A2 << "\n";
    
  // Effective Elastic constants for propagation in Cubic Crystals
  RealType kL_100 = C11;
  RealType kT_100 = C44;
  RealType kL_110 = 0.5 * (C11 + C12 + 2.0*C44);
  RealType kT1_110 = C44;
  RealType kT2_110 = 0.5*(C11 - C12);
  RealType kL_111 = (C11 + 2*C12 + 4*C44) / 3.0;
  RealType kT_111 = (C11 - C12 + C44) / 3.0;
   
}
int main(int argc, char *argv []) {
  std::string method;
  std::string inputFileName;
  std::string outputFileName;

  gengetopt_args_info args_info;

  // parse command line arguments
  if (cmdline_parser(argc, argv, &args_info) != 0) {
    cmdline_parser_print_help();
    exit(1);
  }

  // get input file name
  if (args_info.input_given) {
    inputFileName = args_info.input_arg;
  } else {
    if (args_info.inputs_num){
      inputFileName = args_info.inputs[0];
    }
    else{
      sprintf(painCave.errMsg,
              "No input file name was specified on the command line");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  }

  method = "energy";
  if (args_info.method_given) {    
    method = args_info.method_arg;
    toLower(method);
  }  

  int nMax = args_info.npoints_arg;

  RealType dmax = args_info.delta_arg;
  
  std::vector<Vector6d > Ls;
  Ls.push_back(Vector6d( 1., 1., 1., 0., 0., 0.));
  Ls.push_back(Vector6d( 1., 0., 0., 0., 0., 0.));
  Ls.push_back(Vector6d( 0., 1., 0., 0., 0., 0.));
  Ls.push_back(Vector6d( 0., 0., 1., 0., 0., 0.));
  Ls.push_back(Vector6d( 0., 0., 0., 2., 0., 0.));
  Ls.push_back(Vector6d( 0., 0., 0., 0., 2., 0.));
  Ls.push_back(Vector6d( 0., 0., 0., 0., 0., 2.));
  Ls.push_back(Vector6d( 1., 1., 0., 0., 0., 0.));
  Ls.push_back(Vector6d( 1., 0., 1., 0., 0., 0.));
  Ls.push_back(Vector6d( 1., 0., 0., 2., 0., 0.));
  Ls.push_back(Vector6d( 1., 0., 0., 0., 2., 0.));
  Ls.push_back(Vector6d( 1., 0., 0., 0., 0., 2.));
  Ls.push_back(Vector6d( 0., 1., 1., 0., 0., 0.));
  Ls.push_back(Vector6d( 0., 1., 0., 2., 0., 0.));
  Ls.push_back(Vector6d( 0., 1., 0., 0., 2., 0.));
  Ls.push_back(Vector6d( 0., 1., 0., 0., 0., 2.));
  Ls.push_back(Vector6d( 0., 0., 1., 2., 0., 0.));
  Ls.push_back(Vector6d( 0., 0., 1., 0., 2., 0.));
  Ls.push_back(Vector6d( 0., 0., 1., 0., 0., 2.));
  Ls.push_back(Vector6d( 0., 0., 0., 2., 2., 0.));
  Ls.push_back(Vector6d( 0., 0., 0., 2., 0., 2.));
  Ls.push_back(Vector6d( 0., 0., 0., 0., 2., 2.));
  Ls.push_back(Vector6d( 0., 0., 0., 2., 2., 2.));
  Ls.push_back(Vector6d(-1., .5, .5, 0., 0., 0.));
  Ls.push_back(Vector6d( .5,-1., .5, 0., 0., 0.));
  Ls.push_back(Vector6d( .5, .5,-1., 0., 0., 0.));
  Ls.push_back(Vector6d( 1.,-1., 0., 0., 0., 0.));
  Ls.push_back(Vector6d( 1.,-1., 0., 0., 0., 2.));
  Ls.push_back(Vector6d( 0., 1.,-1., 0., 0., 2.));
  Ls.push_back(Vector6d( .5, .5,-1., 0., 0., 2.));
  Ls.push_back(Vector6d( 1., 0., 0., 2., 2., 0.));
  Ls.push_back(Vector6d( 1., 1.,-1., 0., 0., 0.));
  Ls.push_back(Vector6d( 1., 1., 1.,-2.,-2.,-2.));
  Ls.push_back(Vector6d( .5, .5,-1., 2., 2., 2.));
  Ls.push_back(Vector6d( 0., 0., 0., 2., 2., 4.));
  Ls.push_back(Vector6d( 1., 2., 3., 4., 5., 6.));
  Ls.push_back(Vector6d(-2., 1., 4.,-3., 6.,-5.));
  Ls.push_back(Vector6d( 3.,-5.,-1., 6., 2.,-4.));
  Ls.push_back(Vector6d(-4.,-6., 5., 1.,-3., 2.));
  Ls.push_back(Vector6d( 5., 4., 6.,-2.,-1.,-3.));
  Ls.push_back(Vector6d(-6., 3.,-2., 5.,-4., 1.));
  
  std::vector<int> Lag_strain_list;
  
  if (!method.compare("energy")) {
    for (int i = 1; i < 22; i++) {
      Lag_strain_list.push_back(i);   
    }
  } else {
    for (int i = 35; i < 42; i++) {
      Lag_strain_list.push_back(i);
    }
  }


  //register forcefields, integrators and minimizers
  registerAll();

  // Parse the input file, set up the system, and read the last frame:
  SimCreator creator;
  SimInfo* info = creator.createSim(inputFileName, true);
  Globals* simParams = info->getSimParams();
  ForceManager* forceMan = new ForceManager(info);
  Velocitizer* veloSet = new Velocitizer(info);

  forceMan->initialize();
  info->update();
  
  Shake* shake = new Shake(info);
  bool hasFlucQ = false;
  FluctuatingChargePropagator* flucQ;
    
  if (info->usesFluctuatingCharges()) {
    if (info->getNFluctuatingCharges() > 0) {
      hasFlucQ = true;
      flucQ = new FluctuatingChargeDamped(info);
      flucQ->setForceManager(forceMan);
      flucQ->initialize();
    }
  }    
  
  // Important utility classes for computing system properties:
  Thermo thermo(info);

  // Just in case we were passed a system that is on the move:
  veloSet->removeComDrift();
 
  // std::cerr << "Volume before = " << thermo.getVolume();

  // MinimizerParameters* miniPars = simParams->getMinimizerParameters();
  // OptimizationMethod* minim =OptimizationFactory::getInstance()->createOptimization(toUpperCopy(miniPars->getMethod()), info);
  
  // if (minim == NULL) {
  //   sprintf(painCave.errMsg,
  //           "Optimization Factory can not create %s OptimizationMethod\n",
  //           miniPars->getMethod().c_str());
  //   painCave.isFatal = 1;
  //   simError();
  // }
  
  // BoxObjectiveFunction boxObjf(info, forceMan); 
  // DumpStatusFunction dsf(info);
  // DynamicVector<RealType> initCoords = boxObjf.setInitialCoords();
  // Problem problem(boxObjf, *(new NoConstraint()), dsf, initCoords);
    
  // int maxIter = miniPars->getMaxIterations();
  // int mssIter = miniPars->getMaxStationaryStateIterations();
  // RealType rEps = miniPars->getRootEpsilon();
  // RealType fEps = miniPars->getFunctionEpsilon();
  // RealType gnEps = miniPars->getGradientNormEpsilon();
  
  // EndCriteria endCriteria(maxIter, mssIter, rEps, fEps, gnEps); 

  // minim->minimize(problem, endCriteria);
  // delete minim;

  // std::cerr << "Volume after = " << thermo.getVolume();
  Snapshot* snap = info->getSnapshotManager()->getCurrentSnapshot();
  Mat3x3d refHmat = snap->getHmat();

  forceMan->calcForces();
  Mat3x3d ptRef = thermo.getPressureTensor();
  RealType V0 = thermo.getVolume();
  ptRef.negate();
  ptRef *= Constants::elasticConvert;
  Vector6d stressRef = ptRef.toVoigtTensor();
  
  Vector6d stress(0.0);
  Vector6d strain(0.0);
  Mat3x3d epsilon(0.0);

  std::vector<std::vector<RealType> > stressStrain;
  stressStrain.resize(6);
  std::vector<RealType> strainValues;
  std::vector<RealType> energyValues;
  DynamicRectMatrix<RealType> C(6, 6, 0.0);
  DynamicRectMatrix<RealType> S(6, 6, 0.0);
  
  Vector3d pos;
  SimInfo::MoleculeIterator miter;
  Molecule * mol;
  RealType de;
  Vector3d delta;

  Vector6d Ls_list;
  Vector6d L(0.0);
  Mat3x3d eta(0.0);
  Mat3x3d eps(0.0);
  Mat3x3d x(0.0);
  Mat3x3d test(0.0);
  Mat3x3d deformation(0.0);
  std::vector<RealType> A2;
  RealType norm;
  RealType a, b, c;
  RealType energy;
  Mat3x3d pressureTensor;
  

  for(std::vector<int>::iterator it = Lag_strain_list.begin();
      it !=  Lag_strain_list.end(); ++it) {
    
    Ls_list = Ls[*it];

    // std::cerr << "doing deformation " << Ls_list << "\n";
    strainValues.clear();
    energyValues.clear();
    
    for (int n = 0; n < nMax; n++) {

      de = -0.5*dmax + dmax * RealType(n) / RealType(nMax-1);

      L = Ls_list;
      L *= de;
      eta.setupVoigtTensor(L[0], L[1], L[2], L[3]/2., L[4]/2., L[5]/2.);
      norm = 1.0;
      eps = eta;
      if (eta.frobeniusNorm() > 0.7) {
        std::cerr << "Too large deformation!\n";
      }

      while (norm > 1.0e-10) {
        x = eta - eps*eps / 2.0;
        test = x - eps;
        norm = test.frobeniusNorm();       
        eps = x;
      }
      deformation = SquareMatrix3<RealType>::identity() + eps;
     
      info->getSnapshotManager()->advance();      
      
      for (mol = info->beginMolecule(miter); mol != NULL; 
           mol = info->nextMolecule(miter)) {
        pos = mol->getCom();
        delta = deformation * pos;
        mol->moveCom(delta - pos);
      }

      Mat3x3d Hmat = deformation * refHmat;
      snap->setHmat(Hmat);
      
      shake->constraintR();
      forceMan->calcForces();
      if (hasFlucQ) flucQ->applyConstraints();
      shake->constraintF();

      if (!method.compare("energy")) {
        energy = thermo.getPotential();
        energyValues.push_back(energy);

      } else {
        pressureTensor = thermo.getPressureTensor();
        pressureTensor.negate();
        pressureTensor *= Constants::elasticConvert;
        stress = pressureTensor.toVoigtTensor();

        for (int j = 0; j < 6; j++) {
          stressStrain[j].push_back(stress[j]);
        }

      }
      
      strainValues.push_back(de);
      info->getSnapshotManager()->resetToPrevious();

    }

    if (!method.compare("energy")) {
      quadraticFit(strainValues, energyValues, a, b, c);
      A2.push_back( a * Constants::energyElasticConvert / V0 );
    } else {
    }
  }


  if (!method.compare("energy")) {
    C(0,0) = 2.*A2[0];
    C(0,1) = 1.*(-A2[0]-A2[1]+A2[6]);
    C(0,2) = 1.*(-A2[0]-A2[2]+A2[7]);
    C(0,3) = .5*(-A2[0]-A2[3]+A2[8]) ;
    C(0,4) = .5*(-A2[0]+A2[9]-A2[4]);
    C(0,5) = .5*(-A2[0]+A2[10]-A2[5]);
    C(1,1) = 2.*A2[1];
    C(1,2) = 1.*(A2[11]-A2[1]-A2[2]);
    C(1,3) = .5*(A2[12]-A2[1]-A2[3]);
    C(1,4) = .5*(A2[13]-A2[1]-A2[4]);
    C(1,5) = .5*(A2[14]-A2[1]-A2[5]);
    C(2,2) = 2.*A2[2] ;
    C(2,3) = .5*(A2[15]-A2[2]-A2[3]);
    C(2,4) = .5*(A2[16]-A2[2]-A2[4]);
    C(2,5) = .5*(A2[17]-A2[2]-A2[5]);
    C(3,3) = .5*A2[3];
    C(3,4) = .25*(A2[18]-A2[3]-A2[4]);
    C(3,5) = .25*(A2[19]-A2[3]-A2[5]);
    C(4,4) = .5*A2[4];
    C(4,5) = .25*(A2[20]-A2[4]-A2[5]);
    C(5,5) = .5*A2[5];
  }
    
    
  //     Mat3x3d pt = thermo.getPressureTensor();
  //     pt.negate();
  //     pt *= Constants::elasticConvert;
  //     stress = pt.toVoigtTensor();

  //     for (int j = 0; j < 6; j++) {
  //       stressStrain[j].push_back(stress[j]);
  //     }
      
  //     info->getSnapshotManager()->resetToPrevious();
  //   }
  //   for (int j = 0; j < 6; j++) {
  //     C(i,j) = slope(strainValues, stressStrain[j]);
  //   }
  // }

  // Symmetrize C:

  for (int i = 0; i < 5; i++) {
    for (int j = i+1; j < 6; j++) {
      RealType ctmp = 0.5 * (C(i,j) + C(j,i));
      C(i,j) = ctmp;
      C(j,i) = ctmp;      
    }
  }

  // matrix is destroyed during inversion:
  DynamicRectMatrix<RealType> tmpMat(C);
  invertMatrix(tmpMat, S);

  writeMatrix(C, "Elastic Tensor", "GPa");
  writeMatrix(S, "Compliance Tensor", "GPa^-1");
  
  writeMaterialProperties(C, S);  

  return 0;
}
