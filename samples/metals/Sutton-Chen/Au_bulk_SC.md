#include "metals.md"


component{
  type = "Au";
  nMol = 1372;
}

initialConfig = "./Au_bulk.in";
targetTemp = 400.0;
targetPressure = 1.0;

forceField = "SC";
forceFieldFileName = "SuttonChen.frc";
//forceFieldVariant="SC";

ensemble = "NVE";
dt = 4.0;
runTime = 3e3;


sampleTime = 200.0;
statusTime = 4.0;
seed = 985456376;

tempSet="true";
thermalTime = 500;
useInitialExtendedSystemState="false";
useInitialTime="false";
tauThermostat = 1E3;
tauBarostat = 5E3;
