#include "water.md"


component{
  type = "SSD";
  nMol = 256;
}

initialConfig = "./ssdEM.in";

minimizer = CG;
minimizerMaxIter = 500;
minimizerWriteFrq = 1;
minimizerStepSize = 0.01;
forceField = "DUFF";
electrostaticSummationMethod = "none";
dielectric = 80.0;
cutoffRadius = 9.0;
switchingRadius = 7.0;
