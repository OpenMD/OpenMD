DEV_ROOT=/home/maul/gezelter/gezelter/OOPSE-1.0
IS_UNIX=1

Packages = \
	applications \
	brains \
	constraints \
	integrators \
	io \
	math \
	minimizers \
	primitives \
	profiling \
	restraints \
	types \
	UseTheForce \
	utils \
	visitors \

IncludeDirs = \
	/usr/include 
	
LibraryDirs = \
        /usr/lib 

Libraries = \
        mpich 

include $(DEV_ROOT)/make/Makefile
