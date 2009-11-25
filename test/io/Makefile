Main = IOTest

TestMain = $(Main)

Source = IOTest.cpp \
	 IfstrstreamTestCase.cpp

DEV_ROOT               = ../..
UseMPI                 = yes

ParallelDeclare        = -DIS_MPI
CppCompiler            = icpc
Linker                 = icpc
Delete                 = /bin/rm -rf
MakeDepend             = makedepend

MPIRun                 = /usr/local/mpich/bin/mpirun
MPINumProc             = -np 2

CppOptions             = 
CppParallelOptions     = $(ParallelDeclare)

OpenMDSrcDir             = $(DEV_ROOT)/src
SourceDir                = $(DEV_ROOT)/test
Package                  = $(subst $(shell cd $(SourceDir); pwd)/,,$(shell pwd))
PackageSourceDir         = $(SourceDir)/$(Package)
PackageTargetDir         = $(DEV_ROOT)/obj
PackageParallelTargetDir = $(DEV_ROOT)/MPIobj

MakeOptions            = -k

IncludeDirs = \
	/usr/local/include \
	/usr/local/mpich/include \
        /home/maul/gezelter/tim/cppunit/include \
        $(OpenMDSrcDir)

IncludePath      = -I$(SourceDir) $(IncludeDirs:%=-I%)

LibraryDirs = \
	/usr/local/lib \
	/usr/local/mpich/lib \
        /home/maul/gezelter/tim/cppunit/lib 

LibDirs          = -L$(LibDir)    $(LibraryDirs:%=-L%)

Libraries = \
	-lsprng \
	-lmpich \
	-lfmpich \
        -lcppunit \
        -ldl

#LinkOptions = \
	 -L/usr/local/intel/compiler81/lib -L/usr/lib -lifport -lifcore -limf -lm -lcxa -lirc -lunwind -lirc_s 

#ParallelLinkOptions = \
	 -L/usr/local/intel/compiler81/lib -L/usr/lib -lifport -lifcore -limf -lm -lcxa -lirc -lunwind -lirc_s 


ObjectFiles = $(Source:%.cpp= $(PackageTargetDir)/%.o)
ParallelObjectFiles = $(Source:%.cpp= $(PackageParallelTargetDir)/%.o)

ifneq  "$(words $(ObjectFiles) $(ParallelObjectFiles))" "0"
#  DependencyFile    = $(PackageSourceDir)/Makedepend
endif

ifneq "$(words $(Main))" "0"
    Executable             = $(PackageSourceDir)/$(Main)
    ParallelExecutable     = $(PackageSourceDir)/$(Main)_MPI
endif
#Rules
default : test

# .cpp -> .o
$(PackageTargetDir)/%.o : %.cpp
	@echo $(shell cd $(SourceDir); pwd)
	@echo $(Package)
	$(CppCompiler) $(CppOptions) -c $(IncludePath) $< -o $@

$(PackageParallelTargetDir)/%.o : %.cpp
	$(CppCompiler) $(CppParallelOptions) -c $(IncludePath) $< -o $@

ifeq "$(UseMPI)" "yes"
%.o : %.cpp
	$(MAKE) $(MakeOptions) $(PackageTargetDir)/$@
	$(MAKE) $(MakeOptions) $(PackageParallelTargetDir)/$@
else
%.o : %.cpp
	$(MAKE) $(MakeOptions) $(PackageTargetDir)/$@
endif


# Executable
$(Executable) : $(ObjectFiles)

	$(Linker) $(ObjectFiles) $(LinkOptions) $(LibDirs) $(Libraries) -o $@

$(ParallelExecutable) : $(ParallelObjectFiles)
	$(Linker) $(ParallelObjectFiles) $(ParallelLinkOptions) $(LibDirs) $(Libraries) -o $@


$(DependencyFile): 
	$(Print) $@ 

	$(DEV_ROOT)/scripts/filepp  -I $(DEV_ROOT)/src -od '$$(TargetDir)/' -MM $(Source)>> Make.cpptemp
	@cat Make.cpptemp  >> $(DependencyFile)
	$(Delete) Make.cpptemp

  ifeq "$(UseMPI)" "yes"
	$(DEV_ROOT)/scripts/filepp  -I $(DEV_ROOT)/src -od '$$(ParallelTargetDir)/' $(ParallelDeclare) -MM $(Source)>> Make.cpptemp
	@cat Make.cpptemp  >> $(DependencyFile)
	@$(Delete) Make.cpptemp
  endif





ifeq "$(UseMPI)" "yes"
test : \
	$(DependencyFile) \
	$(ObjectFiles) \
	$(ParallelObjectFiles) \
	$(Executable) \
	$(ParallelExecutable) 
else
test : \
	$(DependencyFile) \
	$(ObjectFiles) \
	$(Executable) 
endif

ifeq "$(UseMPI)" "yes"
runtest : test
	@echo
	@echo test single version
	$(Executable)
	@echo
	@echo test MPI version
	$(MPIRun) $(MPINumProc) $(ParallelExecutable) 
else
runtest : test
	@echo test single version
	$(Executable)
endif
clean:
	$(Delete) $(ObjectFiles) $(ParallelObjectFiles) $(Executable) $(ParallelExecutable)
 
ifdef $(DependencyFile)
-include $(DependencyFile)
endif
