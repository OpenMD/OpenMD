/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include "hydrodynamics/HydroIO.hpp"

#include <fstream>
#include <iomanip>

using namespace std;

namespace OpenMD {

  HydroIO::~HydroIO() {}

  void HydroIO::openWriter(std::ostream& os) {
    std::string h = "OpenMD-Hydro";

    j_[h] = ordered_json::array();
  }

  void HydroIO::writeHydroProp(HydroProp* hp, RealType viscosity,
                               RealType temperature, std::ostream& os) {
    std::string h = "OpenMD-Hydro";
    ordered_json o;

    o["name"]      = hp->getName();
    o["viscosity"] = viscosity;

    Vector3d cor            = hp->getCenterOfResistance();
    o["centerOfResistance"] = {cor[0], cor[1], cor[2]};

    Mat6x6d Xi            = hp->getResistanceTensor();
    o["resistanceTensor"] = json::array();
    for (unsigned int i = 0; i < 6; i++) {
      o["resistanceTensor"][i] = {Xi(i, 0), Xi(i, 1), Xi(i, 2),
                                  Xi(i, 3), Xi(i, 4), Xi(i, 5)};
    }

    Vector3d cod = hp->getCenterOfDiffusion(temperature);
    Mat6x6d Xid  = hp->getDiffusionTensorAtPos(cod, temperature);

    o["temperature"] = temperature;

    o["centerOfDiffusion"] = {cod[0], cod[1], cod[2]};

    o["diffusionTensor"] = json::array();
    for (unsigned int i = 0; i < 6; i++) {
      o["diffusionTensor"][i] = {Xid(i, 0), Xid(i, 1), Xid(i, 2),
                                 Xid(i, 3), Xid(i, 4), Xid(i, 5)};
    }

    Vector3d cop = hp->getCenterOfPitch();
    Mat3x3d pitchAxes;
    Vector3d pitches;
    RealType pitchScalar;

    hp->pitchAxes(pitchAxes, pitches, pitchScalar);

    o["pitch"]         = pitchScalar;
    o["centerOfPitch"] = {cop[0], cop[1], cop[2]};
    o["pitches"]       = {pitches[0], pitches[1], pitches[2]};

    o["pitchAxes"] = json::array();
    for (unsigned int i = 0; i < 3; i++) {
      o["pitchAxes"][i] = {pitchAxes(i, 0), pitchAxes(i, 1), pitchAxes(i, 2)};
    }

    j_[h].push_back(o);
  }

  void HydroIO::closeWriter(std::ostream& os) { os << j_.dump(2) << std::endl; }

  map<string, HydroProp*> HydroIO::parseHydroFile(const string& f) {
    map<string, HydroProp*> props;

    ifstream ifs(f);
    json ij = json::parse(ifs);

    auto& entries = ij["OpenMD-Hydro"];

    for (auto& entry : entries) {
      HydroProp* hp = new HydroProp();
      std::string name;
      Vector3d cor;
      Mat6x6d Xi;

      name = entry["name"].get<std::string>();

      for (unsigned int i = 0; i < 3; i++) {
        cor[i] = entry["centerOfResistance"].get<vector<RealType>>()[i];
      }

      for (unsigned int i = 0; i < 6; i++) {
        for (unsigned int j = 0; j < 6; j++) {
          Xi(i, j) =
              entry["resistanceTensor"].get<vector<vector<RealType>>>()[i][j];
        }
      }

      hp->setName(name);
      hp->setCenterOfResistance(cor);
      hp->setResistanceTensor(Xi);
      props.insert(map<string, HydroProp*>::value_type(name, hp));
    }
    return props;
  }

  void HydroIO::interpretHydroProp(HydroProp* hp, RealType viscosity,
                                   RealType temperature) {
    Vector3d ror = hp->getCenterOfResistance();

    Mat6x6d Xi;
    Xi = hp->getResistanceTensor();
    Mat3x3d Xirtt;
    Mat3x3d Xirrt;
    Mat3x3d Xirtr;
    Mat3x3d Xirrr;

    Xi.getSubMatrix(0, 0, Xirtt);
    Xi.getSubMatrix(0, 3, Xirrt);
    Xi.getSubMatrix(3, 0, Xirtr);
    Xi.getSubMatrix(3, 3, Xirrr);

    Mat6x6d Dr;
    Dr = hp->getDiffusionTensor(temperature);

    Mat3x3d Drtt;
    Mat3x3d Drrt;
    Mat3x3d Drtr;
    Mat3x3d Drrr;

    Dr.getSubMatrix(0, 0, Drtt);
    Dr.getSubMatrix(0, 3, Drrt);
    Dr.getSubMatrix(3, 0, Drtr);
    Dr.getSubMatrix(3, 3, Drrr);

    std::cout << "\n";
    std::cout << "-----------------------------------------\n";
    std::cout << "viscosity = " << viscosity << " Poise" << std::endl;
    std::cout << "temperature = " << temperature << " K" << std::endl;
    std::cout << "-----------------------------------------\n";
    std::cout << "The centers are based on the beads generated by Hydro (.xyz "
                 "file), i.e.,"
              << std::endl;
    std::cout << "not based on the geometry in the .omd file.\n" << std::endl;
    std::cout << "-----------------------------------------\n\n";
    std::cout << "Center of resistance :" << std::endl;
    std::cout << ror << "\n" << std::endl;
    std::cout << "-----------------------------------------\n\n";
    std::cout << "Resistance tensor at center of resistance\n" << std::endl;
    std::cout << "translation [kcal.fs/(mol.A^2)]:" << std::endl;
    std::cout << Xirtt << std::endl;
    std::cout << "rotation-translation [kcal.fs/(mol.A.radian)]:" << std::endl;
    std::cout << Xirtr.transpose() << std::endl;
    std::cout << "translation-rotation [kcal.fs/(mol.A.radian)]:" << std::endl;
    std::cout << Xirtr << std::endl;
    std::cout << "rotation [kcal.fs/(mol.radian^2)]:" << std::endl;
    std::cout << Xirrr << std::endl;
    std::cout << "-----------------------------------------\n\n";
    std::cout << "Diffusion tensor at center of resistance\n" << std::endl;
    std::cout << "translation (A^2 / fs):" << std::endl;
    std::cout << Drtt << std::endl;
    std::cout << "rotation-translation (A.radian / fs):" << std::endl;
    std::cout << Drrt << std::endl;
    std::cout << "translation-rotation (A.radian / fs):" << std::endl;
    std::cout << Drtr << std::endl;
    std::cout << "rotation (radian^2 / fs):" << std::endl;
    std::cout << Drrr << std::endl;
    std::cout << "-----------------------------------------\n\n";

    // calculate center of diffusion using the same arbitrary origin as above
    // (from the generated geometry file .xyz)

    Vector3d cod = hp->getCenterOfDiffusion(temperature);
    Mat6x6d Xid  = hp->getResistanceTensorAtPos(cod);

    Mat3x3d Xidtt;
    Mat3x3d Xidrt;
    Mat3x3d Xidtr;
    Mat3x3d Xidrr;

    Xid.getSubMatrix(0, 0, Xidtt);
    Xid.getSubMatrix(0, 3, Xidrt);
    Xid.getSubMatrix(3, 0, Xidtr);
    Xid.getSubMatrix(3, 3, Xidrr);

    // calculate Diffusion Tensor at center of diffusion
    Mat6x6d Dd = hp->getDiffusionTensorAtPos(cod, temperature);

    Mat3x3d Ddtt;
    Mat3x3d Ddtr;
    Mat3x3d Ddrt;
    Mat3x3d Ddrr;

    Dd.getSubMatrix(0, 0, Ddtt);
    Dd.getSubMatrix(0, 3, Ddrt);
    Dd.getSubMatrix(3, 0, Ddtr);
    Dd.getSubMatrix(3, 3, Ddrr);

    std::cout << "Center of diffusion: " << std::endl;
    std::cout << cod << "\n" << std::endl;
    std::cout << "-----------------------------------------\n\n";
    std::cout << "Diffusion tensor at center of diffusion \n " << std::endl;
    std::cout << "translation (A^2 / fs) :" << std::endl;
    std::cout << Ddtt << std::endl;
    std::cout << "rotation-translation (A.radian / fs):" << std::endl;
    std::cout << Ddtr.transpose() << std::endl;
    std::cout << "translation-rotation (A.radian / fs):" << std::endl;
    std::cout << Ddtr << std::endl;
    std::cout << "rotation (radian^2 / fs):" << std::endl;
    std::cout << Ddrr << std::endl;
    std::cout << "-----------------------------------------\n\n";
    std::cout << "Resistance tensor at center of diffusion \n " << std::endl;
    std::cout << "translation [kcal.fs/(mol.A^2)]:" << std::endl;
    std::cout << Xidtt << std::endl;
    std::cout << "rotation-translation [kcal.fs/(mol.A.radian)]:" << std::endl;
    std::cout << Xidrt << std::endl;
    std::cout << "translation-rotation [kcal.fs/(mol.A.radian)]:" << std::endl;
    std::cout << Xidtr << std::endl;
    std::cout << "rotation [kcal.fs/(mol.radian^2)]:" << std::endl;
    std::cout << Xidrr << std::endl;
    std::cout << "-----------------------------------------\n\n";
  }
}  // namespace OpenMD
