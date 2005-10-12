#include "lj.md"

nComponents = 2;
component{
  type = "Kr";
  nMol = 9;
}
component{
  type = "He";
  nMol = 855;
}

initialConfig = "./cutoff_test.in";

forceField = "LJ";
targetTemp = 50;
targetPressure = 1;
tauThermostat = 1e3;
tauBarostat = 1e3;
cutoffRadius = 10.0;
switchingRadius = 10.0;

ensemble = "NPTi";
dt = 1;
runTime = 1e5;

sampleTime = 1000;
statusTime = 10;

thermalTime = 1000.0;
tempSet = "false";
