<OpenMD version=2>
  <MetaData>


molecule{
  name = "P+";
  
  atom[0]{
    type = "P+";
    position(0.0, 0.0, 0.0);
  }
}

molecule{
  name = "N-";
  
  atom[0]{
    type = "N-";
    position(0.0, 0.0, 0.0);
  }
}

component{
  type = "N-";
  nMol = 1;
}
component{
  type = "P+";
  nMol = 1;
}

ensemble = NVE;
forceField = "Multipole";

cutoffMethod = "shifted_force";
electrostaticScreeningMethod = "damped";
dielectric = 80.0;
cutoffRadius = 9.0;

dt = 1.0;
runTime = 1.0;

sampleTime = 1.0;
statusTime = 1.0;

outputForceVector = true;
outputParticlePotential = true;
outputElectricField = true;
outputFluctuatingCharges = true;

  </MetaData>
  <Snapshot>
    <FrameData>
        Time: 0
        Hmat: {{ 100.0, 0, 0 }, { 0, 100.0, 0 }, { 0, 0, 100.0 }}
  Thermostat: 0 , 0
    Barostat: {{ 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }}
    </FrameData>
    <StuntDoubles>
         0      pv     0 0 0 0 0 0
         1      pv     0 0 4 0 0 0
    </StuntDoubles>
  </Snapshot>
</OpenMD>
