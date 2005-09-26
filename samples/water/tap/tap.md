#include "water.md"

nComponents = 1;
component{
  type = "TAP";
  nMol = 256;
}

initialConfig = "./tap.in";

ensemble = NVT;
forceField = "DUFF";
electrostaticSummationMethod = "none"
dielectric = 80.0;
cutoffRadius = 9.0;
switchingRadius = 7.7;

targetTemp = 300;
targetPressure = 1.0;

tauThermostat = 1e3;
tauBarostat = 1e4;

dt = 2.0;
runTime = 2.5e4;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

//tempSet = "true";
//thermalTime = 10;
sampleTime = 100;
statusTime = 10;
