molecule{
  name = "linear";
  nAtoms = 1;
  atom[0]{
     type="linear";
     position( 0.0, 0.0, 0.0 );
     orientation( 0.0, 0.0, 1.0);
  }
}

nComponents = 1;
component{
  type = "linear";
  nMol = 2;
}

initialConfig = "./linear.in";

forceField = "SHAPES";
targetTemp = 119.8;

ensemble = "NVE";
dt = 1.0;
runTime = 1e2;

sampleTime = 1;
statusTime = 1;
