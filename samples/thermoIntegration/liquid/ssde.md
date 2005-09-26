#include "water.md"

nComponents = 1;
component{
  type = "SSD_E";
  nMol = 256;
}

initialConfig = "./ssde.in";

ensemble = NVT;
forceField = "DUFF";
electrostaticSummationMethod = "none"
dielectric = 80.0;
cutoffRadius = 9.0;
switchingRadius = 7.7;

density = 0.0334;

targetTemp = 300;
targetPressure = 1.0;

tauThermostat = 1e3;
tauBarostat = 5e3;

dt = 2.0;
runTime = 1e3;
useInitialTime = "false";
useInitialExtendedSystemState = "false";
useLiquidThermInt = "true";
thermodynamicIntegrationLambda = 1.0;
thermodynamicIntegrationK = 1.0;

//tempSet = "true";
//thermalTime = 10;
sampleTime = 100;
statusTime = 10;
