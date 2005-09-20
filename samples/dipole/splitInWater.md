#include "splitDipole.md"
#include "water.md"

nComponents = 2;
component{
  type = "HEAD";
  nMol = 1;
}

component{
  type = "SSD";
  nMol = 248;
}

initialConfig = "./splitInWater.in";

ensemble = NVT;
forceField = "DUFF";
forceFieldFileName = "Lipid.frc";
electrostaticSummationMethod = "none";
dielectric = 80.0;
cutoffRadius = 9.0;
switchingRadius = 7.8;

density = 0.0334;

targetTemp = 300;
targetPressure = 1.0;

tauThermostat = 1e3;
tauBarostat = 1e4;

dt = 2.0;
runTime = 1e4;

sampleTime = 100;
statusTime = 2;
