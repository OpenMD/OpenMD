#include "water.md"


component{
  type = "SSD_E";
  nMol = 256;
}

initialConfig = "./ssde.in";

ensemble = NVT;
forceField = "DUFF";
electrostaticSummationMethod = "none";
cutoffRadius = 9.0;
switchingRadius = 7.7;

targetTemp =2.0;
targetPressure = 1.0;

tauThermostat = 1e3;
tauBarostat = 1e4;

dt = 2.0;
runTime = 5e4;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

//tempSet = "true";
//thermalTime = 10;
sampleTime = 1000;
statusTime = 10;
