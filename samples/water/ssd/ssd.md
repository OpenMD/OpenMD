#include "water.md"


component{
  type = "SSD";
  nMol = 256;
}

initialConfig = "./ssd.in";

ensemble = NVE;
forceField = "DUFF";
electrostaticSummationMethod = "none";
cutoffRadius = 9.0;
switchingRadius = 7.65;

targetTemp = 300;
targetPressure = 1.0;

tauThermostat = 1e3;
tauBarostat = 1e4;

dt = 2.0;
runTime = 1e4;

//tempSet = "true";
//thermalTime = 200;
sampleTime = 100;
statusTime = 10;
