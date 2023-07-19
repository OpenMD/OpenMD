/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
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
#include "utils/simError.h"

#include <fstream>
#include <iomanip>

using namespace std;

namespace OpenMD {

  HydroIO::~HydroIO() {}

  void HydroIO::openWriter(std::ostream& os) {
    std::string h = "OpenMD-Hydro";
#if defined(NLOHMANN_JSON)    
    j_[h] = ordered_json::array();
    writerOpen_ = true;
#elif defined(RAPID_JSON)
    osw_ = new OStreamWrapper(os);
    w_.Reset(*osw_);
    writerOpen_ = true; 

    w_.SetMaxDecimalPlaces(7);
    w_.SetIndent(' ', 2);

    w_.StartObject();
    w_.Key(h.c_str());
    w_.StartArray();
#endif
  }

  void HydroIO::writeHydroProp(HydroProp* hp, RealType viscosity,
                               RealType temperature, std::ostream& os) {

    if (!writerOpen_) openWriter(os);

    std::string h = "OpenMD-Hydro";
    std::string name = hp->getName();
    Vector3d cor = hp->getCenterOfResistance();
    Mat6x6d Xi   = hp->getResistanceTensor();
    Vector3d cod = hp->getCenterOfDiffusion(temperature);
    Mat6x6d Xid  = hp->getDiffusionTensorAtPos(cod, temperature);
    Vector3d cop = hp->getCenterOfPitch();
    Mat3x3d pitchAxes;
    Vector3d pitches;
    RealType pitchScalar;

    hp->pitchAxes(pitchAxes, pitches, pitchScalar);

#if defined(NLOHMANN_JSON)    
    
    ordered_json o;
    o["name"] = name;
    o["viscosity"] = viscosity;
    o["centerOfResistance"] = {cor[0], cor[1], cor[2]};
    o["resistanceTensor"] = json::array();
    
    for (unsigned int i = 0; i < 6; i++) {
      o["resistanceTensor"][i] = {Xi(i, 0), Xi(i, 1), Xi(i, 2),
                                  Xi(i, 3), Xi(i, 4), Xi(i, 5)};
    }

    o["temperature"] = temperature;
    o["centerOfDiffusion"] = {cod[0], cod[1], cod[2]};
    o["diffusionTensor"] = json::array();
    
    for (unsigned int i = 0; i < 6; i++) {
      o["diffusionTensor"][i] = {Xid(i, 0), Xid(i, 1), Xid(i, 2),
                                 Xid(i, 3), Xid(i, 4), Xid(i, 5)};
    }
    
    o["pitch"] = pitchScalar;
    o["centerOfPitch"] = {cop[0], cop[1], cop[2]};
    o["pitches"] = {pitches[0], pitches[1], pitches[2]};

    o["pitchAxes"] = json::array();
    for (unsigned int i = 0; i < 3; i++) {
      o["pitchAxes"][i] = {pitchAxes(i, 0), pitchAxes(i, 1), pitchAxes(i, 2)};
    }

    j_[h].push_back(o);
    
#elif defined(RAPID_JSON)

    w_.StartObject();
    w_.Key("name");
    w_.String(name.c_str());

    w_.Key("viscosity");
    w_.Double(viscosity);
    w_.Key("centerOfResistance");
    w_.StartArray();
    w_.SetFormatOptions(kFormatSingleLineArray);

    for (unsigned i = 0; i < 3; i++)
      w_.Double( cor[i] );
    w_.EndArray();
    w_.SetFormatOptions(kFormatDefault);

    w_.Key("resistanceTensor");
    w_.StartArray();
    for (unsigned i = 0; i < 6; i++) {
      w_.StartArray();
      w_.SetFormatOptions(kFormatSingleLineArray);
      
      for (unsigned j = 0; j < 6; j++) {
        w_.Double( Xi(i, j) );
      }
      w_.EndArray();
      w_.SetFormatOptions(kFormatDefault);

    }
    w_.EndArray();

    w_.Key("temperature");
    w_.Double(temperature);
    w_.Key("centerOfDiffusion");
    w_.StartArray();
    w_.SetFormatOptions(kFormatSingleLineArray);

    for (unsigned i = 0; i < 3; i++)
      w_.Double( cod[i] );
    w_.EndArray();
    w_.SetFormatOptions(kFormatDefault);

    w_.Key("diffusionTensor");
    w_.StartArray();
    for (unsigned i = 0; i < 6; i++) {
      w_.StartArray();
      w_.SetFormatOptions(kFormatSingleLineArray);
      
      for (unsigned j = 0; j < 6; j++) {
        w_.Double( Xid(i, j) );
      }
      w_.EndArray();
      w_.SetFormatOptions(kFormatDefault);

    }
    w_.EndArray();

    w_.Key("pitch");
    w_.Double(pitchScalar);
    w_.Key("centerOfPitch");
    w_.StartArray();
    w_.SetFormatOptions(kFormatSingleLineArray);

    for (unsigned i = 0; i < 3; i++)
      w_.Double( cop[i] );
    w_.EndArray();
    w_.SetFormatOptions(kFormatDefault);
    
    w_.Key("pitchAxes");
    w_.StartArray();
    for (unsigned i = 0; i < 3; i++) {
      w_.StartArray();
      w_.SetFormatOptions(kFormatSingleLineArray);
      for (unsigned j = 0; j < 3; j++) {
        w_.Double( pitchAxes(i, j) );
      }
      w_.EndArray();
      w_.SetFormatOptions(kFormatDefault);
    }
    w_.EndArray();
    
    w_.EndObject();
#endif
  }

  void HydroIO::closeWriter(std::ostream& os) {
#if defined(NLOHMANN_JSON)    
    os << j_.dump(2) << std::endl;
#elif defined(RAPID_JSON)
    w_.EndArray();
    w_.EndObject();
    delete osw_;
#endif
    writerOpen_ = false;
  }

  map<string, HydroProp*> HydroIO::parseHydroFile(const string& f) {
    map<string, HydroProp*> props;

    ifstream ifs(f);
    
    if (!ifs.good()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
	       "HydroIO: Cannot open file: %s\n", f.c_str());
      painCave.isFatal = 1;
      simError();
    }
          
#if defined(NLOHMANN_JSON)    
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
#elif defined(RAPID_JSON)
    // Parse entire file into memory once, then reuse d_ for subsequent
    // hydroProps.
    
    if (ifs.peek() != EOF) {
      rapidjson::IStreamWrapper isw(ifs);
      d_.ParseStream(isw);
      if (d_.HasParseError()) {
	snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
		 "HydroIO: JSON parse error in file %s\n"
		 "\tError: %zu : %s\n", f.c_str(), d_.GetErrorOffset(),
		 rapidjson::GetParseError_En(d_.GetParseError()));
	painCave.isFatal = 1;
	simError();	
      }
      if (!d_.IsObject()) {
	snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
		 "HydroIO: OpenMD-Hydro should be a single object.\n");
	painCave.isFatal = 1;
	simError();
      }
      // OpenMD-Hydro has a single object, but check that it's really
      // OpenMD-Hydro
      if (!d_.HasMember("OpenMD-Hydro")) {
	snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
		 "HydroIO: File %s does not have a OpenMD-Hydro object.\n",
		 f.c_str());
	painCave.isFatal = 1;
	simError();
      }
    }
    const Value& entries = d_["OpenMD-Hydro"];
    for (auto& entry : entries.GetArray()) {
      HydroProp* hp = new HydroProp();
      std::string name;
      Vector3d cor;
      Mat6x6d Xi;
      
      name = entry["name"].GetString();
      
      for (unsigned int i = 0; i < 3; i++) {
	cor[i] = entry["centerOfResistance"][i].GetDouble();
      }
      
      for (unsigned int i = 0; i < 6; i++) {
	for (unsigned int j = 0; j < 6; j++) {
	  Xi(i, j) =
	    entry["resistanceTensor"][i][j].GetDouble();
	}
      }
      
      hp->setName(name);
      hp->setCenterOfResistance(cor);
      hp->setResistanceTensor(Xi);
      props.insert(map<string, HydroProp*>::value_type(name, hp));
    }               
#endif
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
    std::cout << "The centers are based on the elements generated by Hydro " << std::endl;
    std::cout << "which have been placed in an .xyz or .stl file." << std::endl;
    std::cout << "They are not based on the geometry in the .omd file.\n" << std::endl;
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
