#include "water.md"

nComponents = 1;
component{
  type = "SSD";
  nMol = 256;
}

initialConfig = "./zconsFixed.in";

ensemble = NVE;
forceField = "DUFF";

zconsTime = 1;
zconsTol = 0.01; 

nZconstraints = 1;
zConstraint[0]{
  molIndex =0;
  kRatio = 0.5;
}
	
cutoffRadius = 9.0;
switchingRadius = 7.8;

electrostaticSummationMethod = "none";

dielectric = 80.0;

density = 0.0334;

targetTemp = 300;
targetPressure = 1.0;

dt = 1;
tauBarostat = 1E5;
tauThermostat = 1E3;

runTime = 1e4;
tempSet = false;
sampleTime = 100;
