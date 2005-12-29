#include "metals.md"


component{
  type = "Au";
	nMol = 9179;
}

initialConfig = "/Au_nanoparticle.in";

forceField = "EAM";
forceFieldVariant = "VC";
targetTemp = 1000.0;


ensemble = "NVE";
dt = 1.0;
runTime = 1000.0;

usePeriodicBoundaryConditions = "false";

sampleTime = 50.0;

thermalTime = 100.0;
tempSet = "false";
