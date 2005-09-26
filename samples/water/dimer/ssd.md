#include "water.md"

nComponents = 1;
component{
  type = "SSD_E";
  nMol = 2;
}

initialConfig = "./ssd.in";


ensemble = NVT;
forceField = "DUFF";
electrostaticSummationMethod = "none"
dielectric = 80.0;
cutoffRadius = 9.0;
switchingRadius = 7.8;

density = 0.0334;

targetTemp = 0.1;
targetPressure = 1.0;

tauThermostat = 1e4;
tauBarostat = 1e4;

dt = 1.0;
runTime = 1e4;

sampleTime = 1e2;
statusTime = 10;
useInitialTime = "false";
useInitialExtendedSystemState = "false";
