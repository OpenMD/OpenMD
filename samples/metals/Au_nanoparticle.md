#include "metals.md"

nComponents = 1;
component{
  type = "Au";
  nMol = 1926;
}

initialConfig = "./Au_nanoparticle.in";

forceField = "EAM";
forceFieldVariant = "VC";
targetTemp = 1000.0;
density = 0.02143659;

ensemble = "NVE";
dt = 1.0;
runTime = 1000.0;

usePeriodicBoundaryConditions = "false";

sampleTime = 50.0;

thermalTime = 100.0;
tempSet = "false";
