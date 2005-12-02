#include "gb.md"
#include "lj.md"

component{
  type = "linear";
  nMol = 1;
}
component{
  type = "Ar";
  nMol =863;
}


initialConfig = "./gb-ar.in";

ensemble = NVE;
forceField = "DUFF";
cutoffPolicy = "traditional";

dt = 1;
runTime = 1e4;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

sampleTime = 100;
statusTime = 1;
