# Compiling OpenMD

OpenMD is written in C++. Compiling is the process of turning this
C++ into instructions that the computer’s processor can understand.

## Requirements

To build OpenMD, you need the following:

* The source code for the latest release of OpenMD
* A C++ compiler
* CMake 3.20 or newer

OpenMD uses CMake as its build system. CMake is an open source
cross-platform build system from KitWare.

You need to install CMake 3.20 or newer. This is available as a
binary package from the KitWare website; alternatively, it may be
available through your package manager (on Linux). If necessary, you
can also compile it yourself from the source code.

The following are optional when compiling OpenMD, but if they are not
available some features will be missing:

* OpenMPI or MPICH – Very good implementations of the MPI-2 specification
for parallel computing.  A version of the MPI library is required
if you want to run the multi-processor version of OpenMD. 

* python - An interpreted scripting language that some of the OpenMD 
utilities use to parse and process data files. Some python
scripts also depend on NumPy and SciPy. You'll want version 3 or higher,
and we test on 3.12

* qhull – A computational geometry toolbox for computing convex
hulls and Delaunay triangulations.  qhull is required for the
LangevinHull integrator and for any of the tools that compute the
Hull atoms or hull volumes of nanoparticles and clusters.  You'll
want version qhull 2020.2 or newer.

* openbabel – a chemical toolbox for converting between different
data formats.  This is required for building the atom2md program
which helps prepare initial "metadata" or md files for
simulations.  You'll want openbabel version 3.1.1 or newer, but note
that when openbabel is found, there are some compiler warnings that
will appear.

* fftw - a library for computing discrete Fourier transforms.  This
is required for surface undulation spectra (Hxy in
staticProps). Get version 3.3.10 or newer.

* zlib - required to support reading gzipped trajectory files (this 
is usually included on most unix / macOS machines)

You’ll also likely want to download and compile the following useful
tools for interacting with the data:

  * [Jmol](https://jmol.sourceforge.net/)
  * [xmgrace/ grace](https://plasma-gate.weizmann.ac.il/Grace/)
  * [NumPy](https://numpy.org/)
  * [SciPy](https://scipy.org/)
  * [vmd](https://www.ks.uiuc.edu/Research/vmd/)

If you are going to be extending or developing OpenMD, you’ll need
the following tools:

* [antlr](https://www.antlr2.org/ – our tool for parsing meta-data files.
You’ll want version 2, not 3.  

* [gengetopt](https://www.gnu.org/software/gengetopt/gengetopt.html) - a
tool to generate C code to parse the command line arguments argc and argv
that are part of every C or C++ program.

## Basic build procedure

The recommended way to build OpenMD is to use a separate source and
build directory; for example, openmd-3.1 and build. The first step
is to create these directories:

```bash
$ tar zxf openmd-3.1.tar.gz   # (this creates openmd-3.1)
$ mkdir build
```

Now you need to run cmake to configure the build. The following will
configure the build to use all of the default options:

```bash
$ cd build
$ cmake ../openmd-3.1
```

If you need to specify a particular compiler, you can do that with
environment variables before the cmake line

```bash
$ export CXX=/opt/local/lib/openmpi/bin/mpic++
$ cmake ../openmd-3.1
```

If you need to specify an option, use the -D switch to cmake. For
example, the following line sets the value of `CMAKE_INSTALL_PREFIX`
and `CMAKE_BUILD_TYPE`:

```bash
$ cmake ../openmd-3.1 -DCMAKE_INSTALL_PREFIX=~/Tools -DCMAKE_BUILD_TYPE=DEBUG
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

By default, OpenMD installs directories into `/usr/local/openmd`,
so in order for the command line options to be in your path, you'll 
need:
```bash
export PATH=${PATH}:/usr/local/openmd/bin
```
or if you use csh or tcsh as your shell:
```
setenv PATH ${PATH}:/usr/local/openmd/bin
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
$ cmake ../openmd-3.1 -DCMAKE_INSTALL_PREFIX=~/Tools/openmd-install
```

Then you can run make and make install without needing root access:

```bash
$ make && make install
```

Once you have installed OpenMD in a specified location, a 
`bin` subdirectory will contain all of the command line tools.
In order for these command line tools to be accessible commands, 
you'll need:
```bash
export PATH=${PATH}:~/Tools/openmd-install/bin
```
or if you use csh or tcsh as your shell:
```
setenv PATH ${PATH}:~/Tools/openmd-install/bin
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
