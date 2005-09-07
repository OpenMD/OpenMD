#include "lj.md"

nComponents = 1;
component{
  type = "Ar";
  nMol = 108;
}

initialConfig = "./argonNVT.in";

forceField = "LJ";
targetTemp = 119.8;
density = 0.02143659;

ensemble = "NVT";
tauThermostat = 1e3;
dt = 1.0;
runTime = 1e5;

sampleTime = 100;
statusTime = 10;

thermalTime = 100.0;
tempSet = "false";
