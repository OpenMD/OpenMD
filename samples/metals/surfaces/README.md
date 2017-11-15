# Sample metal surfaces

Contained here are examples of (111) cut surfaces of the coinage
metals, as well as a few of the catalytically active metals.

## slabBuilder

OpenMD also has a utility script which makes creation of these types
of systems trivial. **slabBuilder** is a python script which generates
*sc*, *bcc*, and *fcc* lattices and cleaves the crystals along a
desired *(hkl)* plane. The systems are then reoriented such that the
cleaved facet is presented to the z-dimension of the simulation
box. **slabBuilder** comes with a help message when passed *-h* or
*--help*.

