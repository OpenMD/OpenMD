#include "water.md"


component{
  type = "SPCE";
  nMol = 256;
}

initialConfig = "./spce.in";

ensemble = NVE;
forceField = "DUFF";
electrostaticSummationMethod = "none";
//electrostaticScreeningMethod = "damped";
dielectric = 80.0;
dampingAlpha = 0.15;
cutoffRadius = 9.0;
//switchingRadius = 9.0;
switchingFunctionType = "fifth_order_polynomial";

targetTemp = 300;
targetPressure = 1.0;

tauThermostat = 1e3;
tauBarostat = 1e4;

dt = 2.0;
runTime = 1e3;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

//tempSet = "true";
//thermalTime = 10;
sampleTime = 100;
statusTime = 2;
