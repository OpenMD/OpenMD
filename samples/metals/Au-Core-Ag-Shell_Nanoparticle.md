#include "metals.md"

nComponents = 2;
component{
  type = "Ag";
  nMol = 1436;
}
component{
  type = "Au";
  nMol = 490;
}

initialConfig = "./Au-Core-Ag-Shell_Nanoparticle.in";

forceField = "EAM";
forceFieldVariant = "u3";
targetTemp = 1000.0;
density = 0.02143659;

ensemble = "NVE";
dt = 1.0;
runTime = 1000.0;

usePeriodicBoundaryConditions = "false";

sampleTime = 5.0;

thermalTime = 100.0;
tempSet = "false";
