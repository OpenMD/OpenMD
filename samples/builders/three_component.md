<OOPSE version=4>
  <MetaData>
#include "lj.md"

component{
  type = "Ar";
  nMol = 1;
}
component{
  type = "Kr";
  nMol = 1;
}
component{
  type = "He";
  nMol = 1;
}


forceField = "LJ";
targetTemp = 135.1344;
tauThermostat = 1e3;
cutoffRadius = 15.0;
switchingRadius = 15.0;

ensemble = "NVT";
dt = 1.0;
runTime = 1e5;

sampleTime = 100;
statusTime = 1;

thermalTime = 1000.0;
tempSet = "true";
  </MetaData>
</OOPSE>
