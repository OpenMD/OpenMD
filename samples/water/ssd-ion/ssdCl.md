#include "water.md"

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
