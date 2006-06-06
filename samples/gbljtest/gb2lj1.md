#include "gb.md"
#include "lj.md"


component{
  type = "GBlinear";
  nMol = 2;
}
component{
  type = "C";
  nMol =1;
}


initialConfig = "./gb2lj1.in";

ensemble = NVT;
targetTemp = 1;
tauThermostat = 1000;
forceField = "DUFF";

cutoffRadius = 20.0;
switchingRadius = 18.0;
dt = 1.0;
runTime = 1e5;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

sampleTime = 100;
statusTime = 100;
