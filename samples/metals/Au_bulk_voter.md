#include "metals.md"

nComponents = 1;
component{
  type = "Au";
  nMol = 1372;
}

initialConfig = "Au_bulk.in";
targetTemp = 500.0;
targetPressure = 1.0;

forceField = "EAM";
forceFieldVariant="VC";

ensemble = "NVE";
dt = 4.0;
runTime = 1e5;


sampleTime = 200.0;
seed = 985456376;

useInitialExtendedSystemState="false";
useInitialTime="false";
tauThermostat = 1E3;
tauBarostat = 5E3;
