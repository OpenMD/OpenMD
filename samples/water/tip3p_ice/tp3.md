#include "water.md"

nComponents = 1;
component{
  type = "TIP3P";
  nMol = 1000;
}

initialConfig = "./tp3.in";

ensemble = NVT;
forceField = "DUFF";
electrostaticSummationMethod = "none";
dielectric = 80.0;
cutoffRadius = 9.0;
switchingRadius = 7.7;

density = 0.0334;

targetTemp = 200;
targetPressure = 1.0;

tauThermostat = 1e3;
tauBarostat = 5e3;

dt = 2.0;
runTime = 1e4;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

//tempSet = "true";
//thermalTime = 10;
sampleTime = 100;
statusTime = 2;
