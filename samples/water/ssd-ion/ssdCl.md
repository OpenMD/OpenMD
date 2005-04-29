#include "water.md"

molecule{
  name = "Cl-";
  nAtoms = 1;
  atom[0]{
    type = "Cl-";
    position(0.0, 0.0, 0.0);
  }
}

molecule{
  name = "Na+";
  nAtoms = 1;
  atom[0]{
    type = "Na+";
    position(0.0, 0.0, 0.0);
  }
}

nComponents = 2;
component{
  type = "SSD_E";
  nMol = 500;
}
component{
  type = "Cl-";
  nMol = 1;
}

initialConfig = "./ssdCl.in";

ensemble = NVT;
forceField = "DUFF";
useReactionField = "false";
dielectric = 80.0;
cutoffRadius = 10.5;
switchingRadius = 8.925;

targetTemp = 298;
targetPressure = 1.0;

tauThermostat = 1e3;
tauBarostat = 1e4;

dt = 2.0;
runTime = 5e3;

//tempSet = "true";
//thermalTime = 200;
sampleTime = 100;
statusTime = 10;
