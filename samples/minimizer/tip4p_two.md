#include "water.md"

nComponents = 1;
component{
  type = "TIP4P";
  nMol = 2;
}

initialConfig = "./tip4p_two.in";

minimizer = CG;
minimizerMaxIter = 5000;
minimizerWriteFrq = 1;
minimizerStepSize = 0.0001;

forceField = "WATER";
useReactionField = "false";
dielectric = 80.0;
cutoffRadius = 12.0;
switchingRadius = 9.0;
