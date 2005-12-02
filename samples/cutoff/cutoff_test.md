#include "lj.md"

component{
  type = "Kr";
  nMol = 9;
}
component{
  type = "He";
  nMol = 855;
}

initialConfig = "./cutoff_test.in";
cutoffPolicy = "mix";

forceField = "LJ";
targetTemp = 5;
targetPressure = 1;
tauThermostat = 1e3;
tauBarostat = 1e3;

ensemble = "NVE";
dt = 2;
runTime = 1e5;

sampleTime = 100;
statusTime = 1;

thermalTime = 1000.0;
tempSet = "false";
