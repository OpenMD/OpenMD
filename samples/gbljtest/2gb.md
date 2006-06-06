#include "gb.md"
#include "lj.md"


component{
  type = "GBlinear";
  nMol = 2;
}


initialConfig = "./2gb.in";

ensemble = NVT;
targetTemp = 10;
tauThermostat = 30;
forceField = "DUFF";

cutoffRadius = 20.0;
switchingRadius = 18.0;
dt = 1.0;
runTime = 1e5;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

thermalTime = 10;
sampleTime = 100;
statusTime = 1;
