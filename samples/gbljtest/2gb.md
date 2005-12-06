#include "gb.md"
#include "lj.md"


component{
  type = "linear";
  nMol = 2;
}


initialConfig = "./2gb.in";

ensemble = NVT;
targetTemp = 1;
tauThermostat = 1000;
forceField = "DUFF";

cutoffRadius = 20.0;
switchingRadius = 18.0;
dt = 1.0;
runTime = 6e4;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

thermalTime = 10;
sampleTime = 100;
statusTime = 1;
