#include "metals.md"


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
runTime = 6e3;


sampleTime = 200.0;
statusTime = 40.0;
seed = 985456376;

useInitialExtendedSystemState="false";
useInitialTime="false";
tauThermostat = 1E3;
tauBarostat = 5E3;
