#include "gb.md"
#include "lj.md"

nComponents = 2;
component{
  type = "C";
  nMol =1;
}
component{
  type = "linear";
  nMol = 1;
}


initialConfig = "./gbRotRoundY.in";

ensemble = NVE;
forceField = "DUFF";

cutoffRadius = 20.0;
switchingRadius = 18.0;
dt = 1.0;
runTime = 1e5;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

thermalTime = 10;
sampleTime = 100;
statusTime = 100;
