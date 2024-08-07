# Compiling OpenMD

OpenMD is written in C++. Compiling is the process of turning this
C++ into instructions that the computer’s processor can understand.

## Requirements

To build OpenMD, you need the following:

* The source code for the latest release of OpenMD
* A C++ compiler
* CMake 3.15 or newer

OpenMD uses CMake as its build system. CMake is an open source
cross-platform build system from KitWare.

You need to install CMake 3.15 or newer. This is available as a
binary package from the KitWare website; alternatively, it may be
available through your package manager (on Linux). If necessary, you
can also compile it yourself from the source code.

The following are optional when compiling OpenMD, but if they are not
available some features will be missing:

* OpenMPI – A very good implementation of the MPI-2 specification
for parallel computing.  A version of the MPI library is required
if you want to run the multi-processor version of OpenMD

* python - An interpreted scripting language that some of the OpenMD 
utilities use to parse and process data files. Some python
scripts also depend on NumPy and SciPy.

* qhull – A computational geometry toolbox for computing convex
hulls and Delaunay triangulations.  qhull is required for the
LangevinHull integrator and for any of the tools that compute the
Hull atoms or hull volumes of nanoparticles and clusters.

* openbabel – a chemical toolbox for converting between different
data formats.  This is required for building the atom2md program
which helps prepare initial "metadata" or md files for
simulations.

* fftw - a library for computing discrete Fourier transforms.  This
is required for surface undulation spectra (Hxy in
staticProps). Get version 3.

* zlib - required to support reading gzipped trajectory files

You’ll also likely want to download and compile the following useful
tools for interacting with the data:

  * Jmol
  * xmgr
  * grace
  * NumPy
  * SciPy
  * vmd

If you are going to be extending or developing OpenMD, you’ll need
the following tools:

* antlr – our tool for parsing meta-data files.  You’ll want
  version 2, not 3.  

* gengetopt - a tool to generate C code to parse the command line
  arguments argc and argv that are part of every C or C++ program

## Basic build procedure

The recommended way to build OpenMD is to use a separate source and
build directory; for example, openmd-3.0 and build. The first step
is to create these directories:

```bash
$ tar zxf openmd-3.0.tar.gz   # (this creates openmd-3.0)
$ mkdir build
```

Now you need to run cmake to configure the build. The following will
configure the build to use all of the default options:

```bash
$ cd build
$ cmake ../openmd-3.0
```

If you need to specify a particular compiler, you can do that with
environment variables before the cmake line

```bash
$ export CXX=/opt/local/lib/openmpi/bin/mpic++
$ cmake ../openmd-3.0
```

If you need to specify an option, use the -D switch to cmake. For
example, the following line sets the value of `CMAKE_INSTALL_PREFIX`
and `CMAKE_BUILD_TYPE`:

```bash
$ cmake ../openmd-3.0 -DCMAKE_INSTALL_PREFIX=~/Tools -DCMAKE_BUILD_TYPE=DEBUG
```

We will discuss various possible options later.

At this point, it would be a good idea to compile OpenMD:

```bash
$ make
```

Have a coffee while the magic happens. If you have a multi-processor
machine and would prefer an espresso, try a parallel build instead:

```bash
$ make -j 4  
```

And finally, as root (or using sudo) you should install it:

```bash
$ umask 0022 && make install
```

Or:
  
```bash
$ umask 0022
$ sudo make install
```

### Local build

With the right sort of environment variable magic (see below), you
can actually use OpenMD straight from the build folder. But life is
a bit easier if you install it somewhere, either system-wide or
locally.

By default, OpenMD is installed in /usr/local on a Unix-like
system. This requires root access (or sudo). Even if you do have
root access, you may not want to overwrite an existing installation
or you may want to avoid conflicts with a version of OpenMD
installed by your package manager.

The solution to all of these problems is to do a local install into
a directory somewhere in your home folder. An additional advantage
of a local install is that if you ever want to uninstall it, all you
need to do is delete the installation directory; removing the files
from a global install is more work.

To configure cmake to install into `~/Tools/openmd-install`, for
example, you would do the following:

```bash
$ cmake ../openmd-3.0 -DCMAKE_INSTALL_PREFIX=~/Tools/openmd-install
```

Then you can run make and make install without needing root access:

```bash
$ make && make install
```

### Troubleshooting build problems

* CMake caches some variables from run-to-run. How can I wipe the
  cache to start from scratch?

    > Delete CMakeCache.txt in the build directory. This is also a very
  useful file to look into if you have any problems.

* What environment variables affect how OpenMD finds force field and
  data files?

  > `FORCE_PARAM_PATH` - Used to find the location of the data files
                        used for force fields and atom sizes, etc.
  >
  > If you get errors about not being able to find some .txt files,
  then you should set this to the name of the folder containing
  files such as Amber.frc and element.txt. These are typically
  installed to `/usr/local/openmd/forceFields`.

* CMake honors user umask for creating directories as of now:
  [http://public.kitware.com/Bug/view.php?id=9620](http://public.kitware.com/Bug/view.php?id=9620).
  To get predictable results please set umask explicitly before
  running the make install command.

### Advanced build options 

* How do I do a debug build?

  >`-DCMAKE_BUILD_TYPE=Debug` does a debug build (gcc -g). 
  To revert to a regular build use<br>
  `-DCMAKE_BUILD_TYPE=Release`.

* How do I see what commands cmake is using to build?

  > Run Make as follows:
  > 
  > ```bash
  > $ VERBOSE=1 make
  > ```

* How do I build the Doxygen documentation?

  > If CMake found the "doxygen" program in your PATH, an optional
  build target called "docs" is created.  If the Doxygen executable
  was not on the PATH, you will need to specify its location with<br>
  `-DDOXYGEN_EXECUTABLE=wherever`.  To build the documentation, type:
  > 
  > ```bash
  > $ make docs
  > ```
