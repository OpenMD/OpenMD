DEV_ROOT=/home/maul/gezelter/tim/code/OOPSE-2.0
IS_UNIX=1

Packages = \
	utils \
	math \
	types \
	primitives \
	visitors \
	UseTheForce/DarkSide \
	UseTheForce \
	brains \
	io \
	integrators \
	minimizers \
	constraints \
	profiling \
	restraints \
	applications \

IncludeDirs = \
	/usr/include \
	/usr/local/include \
	
LibraryDirs = \
        /usr/lib 

Libraries = \
        mpich 

include $(DEV_ROOT)/make/Makefile
