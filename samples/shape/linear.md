molecule{
  name = "linear";
  
  atom[0]{
     type="linear";
     position( 0.0, 0.0, 0.0 );
     orientation( 0.0, 0.0, 1.0);
  }
}
molecule{
  name = "Ar";
  
  atom[0]{
     type="Ar";
     position( 0.0, 0.0, 0.0 );
     orientation( 0.0, 0.0, 1.0);
  }
}


component{
  type = "linear";
  nMol = 1;
}

initialConfig = "./linear.in";

forceField = "SHAPES";
cutoffRadius = 12.0;
switchingRadius = 10.2;

targetTemp = 119.8;

ensemble = "NVE";
dt = 1.0;
runTime = 1e3;

sampleTime = 1;
statusTime = 1;
