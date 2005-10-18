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
cutoffPolicy = "traditional";

forceField = "LJ";

ensemble = "NVE";

dt = 1;
runTime = 1e3;
sampleTime = 10;
statusTime = 1;
