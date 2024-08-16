# Metallic nanoparticles


## Background Information

In these examples, we’ll build gold nanoparticles and equilibrate them to a temperature of 300K. 

## Instructions
### Building a spherical metal nanoparticle

1. We need to start with an input `.omd` file for the metal of your choice. This directory has a `gold.omd` for this example. It looks like:
```
<OpenMD>
  <MetaData>
molecule{
  name = "Au";

  atom[0]{
    type = "Au";
    position(0.0, 0.0, 0.0);
  }
}

component{
  type = "Au";
  nMol = 1;
}

forceField = "SC";
forceFieldFileName = "SuttonChen.QSC.frc";
  </MetaData>
</OpenMD>
```

2. Now to build the nanoparticle. We’ll choose a radius of 15 Å and use the lattice constant of gold, 4.08 Å. We’ll call our initial structure `NP15.omd`. 
```
nanoparticleBuilder --latticeConstant=4.08 --radius=15 gold.omd -o NP15.omd
```

3. To take a look at the structure we’ve just created, you can use the following command:
```
Dump2XYZ -i NP15.omd
```
to create a file called `NP15.xyz`, which can be viewed in VMD, Jmol, or any other chemical structure viewer.

4. Add the following lines to the new NP15.omd file below the `forceFieldFileName` line. This sets some details for the simulation. 
```
ensemble = "LHull";
targetTemp = 5;
targetPressure = 1;
viscosity = 0.1;
dt = 4.0;
runTime = 2E5;

sampleTime = 2000.0;
statusTime = 4;
seed = 985456376;

usePeriodicBoundaryConditions = "false";
tauThermostat = 1E3;
tauBarostat = 5E3;
```

5. NanoparticleBuilder carves a nanoparticle of our chosen radius out of a perfect gold crystal. We need to give the atoms some initial velocities before we start equilibrating. We’ll start it out at 5K:
```
thermalizer -t 5 NP15.omd -o NP15_5K.omd
```

6. For the first step in the equilibration we need to let the gold lattice structurally relax. `NP15_5K.omd` can now be run:
```
openmd NP15_5K.omd
```

7. Running the simulation will create several new files. `NP15_5K.dump` contains the trajectory of the simulation. Statistics such as temperature, pressure, and energy will be recorded in the `NP15_5K.stat` file and can be viewed using:
```
xmgrace -nxy NP15_5K.stat
```

8. The end-of-run file `NP15_5K.eor` stores the last configuration of the simulation. We’ll copy it to a new .omd file.
```
cp NP15_5K.eor NP15_100K.omd
```

9. To continue with the equilibration we need to change the `targetTemp` of `NP15_100K.omd`. We’ll increase it to 100 and run the `NP15_100K.omd` file.
    
10. We’ll continue the procedure of copying the .eor file to a new .omd file and increasing the temperature until we’ve reached 300K. Temperature increases of 50 – 100K and simulation times of 100 – 200 ps are reasonable.

### Building an icosahedral metal nanoparticle

1. To start the process, a metal file (`gold.omd`) is needed to describe the material composition of the particle.  We'll use the same file from the spherical example above.

2. To create the particle coordinates, use the command:
```
icosahedralBuilder -o file.omd --shell=8 --latticeConstant=4.08 --ico gold.omd
```
where `file.omd` is the name of your icosahedra, 8 is the number of shells (in this example), 4.08 is the lattice constant in Angstroms for the metal (4.08 is the correct value for gold), and the skeletal OpenMD file above is called `gold.omd`

3. To make sure this icosahedra is the right size for your purposes; use Dump2XYZ to create an xyz file:
```
Dump2XYZ -i file.omd
```
and view the file in VMD or Jmol (where you might want to measure the diameter of the particle).

4. To run the icosahedra in a Langevin Hull and heat the system, we’ll start by inserting the following excerpt after the `forceFieldFileName` in `file.omd`.
```
ensemble = "LHull";
targetTemp = 5;
targetPressure = 1;
viscosity = 0.1;
dt = 4.0;
runTime = 2E5;
sampleTime = 2000.0;
statusTime = 4;
seed = 985456376;
usePeriodicBoundaryConditions = "false";
tauThermostat = 1E3;
tauBarostat = 5E3;
```

5. Before we can equilibrate the system, the atoms need initial velocities, which we can set using the thermalizer command.
```
thermalizer -t 5 -o file-5K.omd file.omd
```
where -t is followed by the desired temperature in Kelvin, and -o is followed by the output file name.

6. To relax the particle at 5K, we’ll run the simulation described in file-5K.omd using the openmd command:
```
openmd file-5K.omd
```

7.  To ensure that the temperature of the system has reached 5K, we can use xmgrace to view the stat file:
```
xmgrace -nxy file-5K.stat
```
In a typical stat file, set 3 (the fourth column) contains the information about the temperature of the system. You can view only the temperature for the trajectory by double clicking any line and bringing up the Set Appearance box. All the sets are initially displayed; turn off the sets you don’t want to see by highlighting and right clicking the set to click the Hide option.

8. The "end of run" file that has been relaxed can then be copied to a new simulation that will be heated to 100K.
```
cp file-5K.eor file-100K.omd
```
Then change the `targetTemp` from 5 to 100.

9. We’ll continue the procedure of copying the .eor file to a new .omd file and increasing the temperature until we’ve reached 250K. Temperature increases of 50 – 100K and simulation times of 100 – 200 ps are reasonable. Heating the system in intervals is relatively effective in keeping the geometry of the particle stable.

## Expected Output
