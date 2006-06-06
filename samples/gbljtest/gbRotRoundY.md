#include "gb.md"
#include "lj.md"


component{
  type = "C";
  nMol =1;
}
component{
  type = "GBlinear";
  nMol = 1;
}


initialConfig = "./gbRotRoundY.in";

ensemble = NVE;
targetTemp = 1;
tauThermostat = 1000;
forceField = "DUFF";

cutoffRadius = 20.0;
switchingRadius = 18.0;
dt = 1;
runTime = 100000;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

thermalTime = 10;
sampleTime = 100;
statusTime = 1;
