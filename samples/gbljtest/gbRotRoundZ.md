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


initialConfig = "./gbRotRoundZ.in";

ensemble = NVT;
targetTemp = 0.001;
tauThermostat = 1000;
forceField = "DUFF";

cutoffRadius = 20.0;
switchingRadius = 18.0;
dt = 1.0;
runTime = 1e5;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

//thermalTime = 10;
sampleTime = 100;
statusTime = 100;
