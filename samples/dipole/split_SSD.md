#include "splitDipole.md"
#include "water.md"

nComponents = 2;
component{
  type = "HEAD";
  nMol = 1;
}

component{
  type = "SSD";
  nMol = 1;
}

initialConfig = "./split_SSD.in";


ensemble = NVT;
forceField = "DUFF";
forceFieldFileName = "Lipid.frc";
useReactionField = "false";
dielectric = 80.0;
cutoffRadius = 9.0;
switchingRadius = 7.8;

density = 0.0334;

targetTemp = 100.0;
targetPressure = 1.0;

tauThermostat = 1e3;
tauBarostat = 1e4;

dt = 1.0;
runTime = 1e5;

sampleTime = 1e2;
statusTime = 10;
useInitialTime = "false";
useInitialExtendedSystemState = "false";
