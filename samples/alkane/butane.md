#include "alkanes.md"

nComponents = 1;
component{
  type = "butane";
  nMol = 32;
}

initialConfig = "./butane.in";
dt = 1.0;
forceField = "DUFF";
ensemble = NVE;
runTime = 10000.0;
sampleTime = 10.0;
