#include "water.md"

nComponents = 1;
component{
  type = "SSD";
  nMol = 256;
}

initialConfig = "./ssd.in";

ensemble = NPTi;
forceField = "DUFF";
useReactionField = "false";
dielectric = 80.0;
cutoffRadius = 9.0;
switchingRadius = 7.8;

density = 0.0334;

targetTemp = 300;
targetPressure = 1.0;

tauThermostat = 1e3;
tauBarostat = 1e4;

dt = 2.0;
runTime = 1e3;

sampleTime = 200;
statusTime = 2;
