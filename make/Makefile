#---------------------------------------------------------------------------
#  (C) 1999 - 2002 Jacob Dreyer - Geotechnical Software Services
#  jacob.dreyer@geosoft.no - http://geosoft.no
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, 
#  MA  02111-1307, USA.
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#
#  GnuMake crash course:
#
#  target : depends
#    rule
#
#  target  - the parameter given to make. I.e. what to build
#  depends - file or other targets target depends on
#  rule    - how to create target (note that rule is preceeded by a TAB char)
#  $(VAR)  - environment variable or variable defined above
#  $@      - Current target
#  $*      - Current target without extension
#  $<      - Current dependency
#
#---------------------------------------------------------------------------
IS_UNIX=1
#packages to be made
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
	applications/oopse \
	applications/dump2Xyz \
	applications/simpleBuilder \

#packages contain libraries
PackageLibs = \
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

IncludeDirs = \
	/usr/local/include \
	/usr/local/mpich/include

LibraryDirs = \
	/usr/local/lib \
	/usr/local/mpich/lib 

Libraries = \
	-lsprng \
	-lmpich \
	-lfmpich

OopseHome       = /usr/local/oopse
ForceParamDir   = /usr/local/oopse/share/forceFields/
SampleSimDir    = /usr/local/oopse/share/samples/
FrcDeclare      = -DFRC_PATH="$(ForceParamDir)"
ParallelDeclare = -DIS_MPI
UseMPI          = yes

#---------------------------------------------------------------------------
#
#  Directories
#
#---------------------------------------------------------------------------

SourceDir         = $(DEV_ROOT)/src
TargetDir         = $(DEV_ROOT)/obj
ParallelTargetDir = $(DEV_ROOT)/MPIobj
LibDir            = $(DEV_ROOT)/lib
ParallelLibDir    = $(DEV_ROOT)/MPIlib
MakeDir           = $(DEV_ROOT)/make
BinDir            = $(DEV_ROOT)/bin
DocsDir           = $(DEV_ROOT)/docs
CurrentDir        = $(CURDIR)
RPath = $(shell cd $(DEV_ROOT); pwd)/lib:$(OopseHome)/lib
ParallelRPath = $(shell cd $(DEV_ROOT); pwd)/MPIlib:$(OopseHome)/MPIlib
LinkOptions = \
	 -L/usr/local/intel/compiler81/lib -L/usr/lib -lifport -lifcore -limf -lm -lcxa -lirc -lunwind -lirc_s \
	-Wl,-R$(RPath)

ParallelLinkOptions = \
	 -L/usr/local/intel/compiler81/lib -L/usr/lib -lifport -lifcore -limf -lm -lcxa -lirc -lunwind -lirc_s \
	-Wl,-R$(ParallelRPath)

ifdef Source
#get the relative path of current package to source directory
# /home/maul/gezelter/src/code/src/UseTheForce/Darkside --> UseTheForce/Darkside
#Package          = $(shell echo $(CurrentDir) | sed -e 's/^.*\/src\/\(.*\)/\1/g')
#use shell script to get the absolute path and then rip it off from $(CurrentDir)
Package          = $(subst $(shell cd $(SourceDir); pwd)/,,$(CurrentDir))

PackageList              = $(Package)
PackageSourceDir         = $(SourceDir)/$(Package)
PackageTargetDir         = $(TargetDir)
PackageParallelTargetDir = $(ParallelTargetDir)
JavaMainClass            = $(subst /,.,$(Package)).$(Main)
else
PackageList              = $(Packages) $(JavaPackages)
endif

PackageListLoop          = $(patsubst %,$(SourceDir)/%/.loop,$(PackageList))

JRE                      = $(JAVA_HOME)/jre/lib/rt.jar

ifdef IS_UNIX
X = :
else
X = \;
endif

#---------------------------------------------------------------------------
#
#  Classification of files
#
#---------------------------------------------------------------------------

# Source
JavaFiles              = $(filter %.java,  $(Source))
CppFiles               = $(filter %.cpp,   $(Source))
CFiles                 = $(filter %.c,     $(Source))
FortranFiles           = $(filter %.f,     $(Source))
F90Files               = $(filter %.F90,   $(Source))
CorbaFiles             = $(filter %.idl,   $(Source))
LexFiles               = $(filter %.l,     $(Source))
YaccFiles              = $(filter %.y,     $(Source))
OtherSourceFiles       = $(filter-out $(JavaFiles) $(CppFiles) $(CFiles) \
		                      $(FortranFiles) $(F90Files) $(LexFiles) \
	                              $(YaccFiles) $(CorbaFiles), \
	                              $(Source))
ManifestFile           = $(PackageSourceDir)/Manifest

SourceFiles            = $(JavaFiles)\
                         $(CppFiles)\
                         $(CFiles)\
                         $(FortranFiles)\
                         $(F90Files)\
                         $(LexFiles)\
                         $(YaccFiles)

# Target
JavaClassFiles         = $(JavaFiles:%.java=  $(PackageTargetDir)/%.class)
JavaClassFilesRel      = $(JavaFiles:%.java=  $(Package)/%.class)
RmiStubFiles           = $(RmiSource:%.java=  $(PackageTargetDir)/%_Stub.class)
RmiSkeletonFiles       = $(RmiSource:%.java=  $(PackageTargetDir)/%_Skel.class)
JniClassFiles          = $(JniSource:%.java=  $(PackageTargetDir)/%.class)
JniHeaders             = $(JniSource:%.java=  $(PackageSourceDir)/%.h)
ObjectFiles            = $(CFiles:%.c=        $(PackageTargetDir)/%.o)\
                         $(CppFiles:%.cpp=    $(PackageTargetDir)/%.o)\
		         $(FortranFiles:%.f=  $(PackageTargetDir)/%.o)\
		         $(F90Files:%.F90=    $(PackageTargetDir)/%.o)\
		         $(LexFiles:%.l=      $(PackageTargetDir)/%.o)\
		         $(YaccFiles:%.y=     $(PackageTargetDir)/%.o)
ParallelObjectFiles    = $(CFiles:%.c=        $(PackageParallelTargetDir)/%.o)\
		         $(CppFiles:%.cpp=    $(PackageParallelTargetDir)/%.o)\
		         $(FortranFiles:%.f=  $(PackageParallelTargetDir)/%.o)\
		         $(F90Files:%.F90=    $(PackageParallelTargetDir)/%.o)\
		         $(LexFiles:%.l=      $(PackageParallelTargetDir)/%.o)\
		         $(YaccFiles:%.y=     $(PackageParallelTargetDir)/%.o)

DerivedSource          = $(YaccFiles:%.y=    %.h) \
	                 $(YaccFiles:%.y=    %.c) \
	                 $(LexFiles:%.l=     %.c)

DerivedCFiles          = $(YaccFiles:%.y=    %.c) \
	                 $(LexFiles:%.l=     %.c)

OtherTargetFiles       = $(OtherSourceFiles:%=$(PackageTargetDir)/%)

ThirdPartyJarsTmp = $(patsubst %,$(LibDir)/%,$(JavaLibraries))
ThirdPartyJars    = $(subst $(Space),$(X),$(ThirdPartyJarsTmp))

ifneq "$(words $(JavaFiles))" "0"
  JavaPackageName   = $(subst /,.,$(Package))
  JarFile           = $(LibDir)/$(subst /,,$(Package)).jar
endif

#if Main is defined, do not build library. It may not be true sometimes
ifneq  "$(words $(ObjectFiles) $(ParallelObjectFiles))" "0"
  DependencyFile    = $(PackageSourceDir)/Makedepend
  ifneq "$(Main)" ""
    Executable             = $(BinDir)/$(Main)
    ParallelExecutable     = $(BinDir)/$(Main)_MPI
  else
    SharedLibrary          = $(LibDir)/lib$(subst /,,$(patsubst %,oopse_%,$(Package))).so
    StaticLibrary          = $(LibDir)/lib$(subst /,,$(patsubst %,oopse_%,$(Package))).a
    ParallelSharedLibrary  = $(ParallelLibDir)/lib$(subst /,,$(patsubst %,oopse_%,$(Package))).so
    ParallelStaticLibrary  = $(ParallelLibDir)/lib$(subst /,,$(patsubst %,oopse_%,$(Package))).a
  endif
endif
#
# Misc
#
ClassPath        = $(JRE)$(X)$(TargetDir)$(X)$(ThirdPartyJars)
JavaPackageNames = $(subst /,.,$(JavaPackages))
IncludePath      = -I$(SourceDir) $(IncludeDirs:%=-I%)
LibDirs          = -L$(LibDir)    $(LibraryDirs:%=-L%)
ParallelLibDirs  = -L$(ParallelLibDir)  $(LibraryDirs:%=-L%)
LocalLibs        = $(subst /,,$(patsubst %, oopse_%, $(PackageLibs))) 
ParallelLocalLibs= $(subst /,,$(patsubst %, oopse_%, $(PackageLibs))) 
LibList          = $(LocalLibs:%=-l%) $(Libraries)
ParallelLibList  = $(ParallelLocalLibs:%=-l%) $(Libraries)


#---------------------------------------------------------------------------
#
#  Tools & Options
#
#---------------------------------------------------------------------------
Print                  = @echo
Move                   = mv
Copy                   = cp
CCompiler              = icc
CppCompiler            = icpc
Linker                 = icpc
MakeDepend             = makedepend
Install                = /usr/bin/install -c
InstallProgram         = ${INSTALL}
InstallData            = ${INSTALL} -m 644
MakeDir                = /home/maul/gezelter/tim/code/OOPSE-2.0/ac-tools/shtool mkdir -p -f
Delete                 = rm -fr
StaticArchiver         = ar
DynamicArchiver        = icc
FortranCompiler        = @FC@
F90Compiler            = ifort
JavaCompiler           = $(JAVA_HOME)/bin/javac
JavaArchiver           = $(JAVA_HOME)/bin/jar
JarSigner              = $(JAVA_HOME)/bin/jarsigner
JavadocGenerator       = $(JAVA_HOME)/bin/javadoc
JniCompiler            = $(JAVA_HOME)/bin/javah
RmiCompiler            = $(JAVA_HOME)/bin/rmic
JavaExecute            = $(JAVA_HOME)/bin/java
Purify                 = purify
WordCount              = wc
List                   = cat
Yacc                   = bison -y
Lex                    = flex


MakeOptions            = -k
MakeDependOptions      = 
StaticArchiverOptions  = rc
DynamicArchiverOptions = -shared
JavaArchiverOptions    =
JniOptions             =
RmiOptions             = -d $(TargetDir) -classpath $(ClassPath) \
			 -sourcepath $(SourceDir)
COptions               = $(FrcDeclare)
CParallelOptions       = $(FrcDeclare) $(ParallelDeclare)
CppOptions             = $(FrcDeclare)
CppParallelOptions     = $(FrcDeclare) $(ParallelDeclare)
FortranOptions         =
F90Options             =  -fpp1  -O -I$(SourceDir) -module $(TargetDir) 
F90ParallelOptions     =  -fpp1  -O -I$(SourceDir) -module $(ParallelTargetDir)   $(ParallelDeclare)
JavaCompilerOptions    = -d $(TargetDir) -classpath $(ClassPath) \
			 -sourcepath $(SourceDir) -deprecation
JavaRunOptions         = -classpath $(ClassPath) 
PurifyOptions          = 
JavadocOptions         = -d $(DocsDir) \
			 -sourcepath $(SourceDir) \
			 -classpath $(ClassPath) \
			 -author \
			 -package \
			 -use \
			 -splitIndex \
			 -version \
			 -link file:$(JAVA_HOME)/docs/api \
			 -windowtitle $(JavadocWindowTitle) \
			 -doctitle $(JavadocDocTitle) \
			 -header $(JavadocHeader) \
			 -bottom $(JavadocFooter)
WordCountOptions       = --lines

Empty                  =
Space                  = $(Empty) $(Empty)


#---------------------------------------------------------------------------
#
# Rules
#
#---------------------------------------------------------------------------
default : build

%.loop :
	@$(MAKE) $(MakeOptions) -C $(subst .loop,,$@) _$(MAKECMDGOALS)all

# Create target directory
$(PackageTargetDir) :
	$(MakeDir) $@

# .c -> .o
$(PackageTargetDir)/%.o : %.c
	$(Print) $@
	$(CCompiler) $(COptions) -c $(IncludePath) $< -o $@

$(PackageParallelTargetDir)/%.o : %.c
	$(Print) $@
	$(CCompiler) $(CParallelOptions) -c $(IncludePath) $< -o $@

ifeq "$(UseMPI)" "yes"
%.o : %.c
	$(MAKE) $(MakeOptions) $(PackageTargetDir)/$@
	$(MAKE) $(MakeOptions) $(PackageParallelTargetDir)/$@
else
%.o : %.c
	$(MAKE) $(MakeOptions) $(PackageTargetDir)/$@
endif

# .cpp -> .o
$(PackageTargetDir)/%.o : %.cpp
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

# .f -> .o
$(PackageTargetDir)/%.o : %.f
	$(FortranCompiler) $(FortranOptions) -c $< -o $@

$(PackageParallelTargetDir)/%.o : %.f
	$(FortranCompiler) $(FortranParallelOptions) -c $< -o $@

ifeq "$(UseMPI)" "yes"
%.o : %.f
	$(MAKE) $(MakeOptions) $(PackageTargetDir)/$@
	$(MAKE) $(MakeOptions) $(PackageParallelTargetDir)/$@
else
%.o : %.f
	$(MAKE) $(MakeOptions) $(PackageTargetDir)/$@
endif

# .F90 -> .o
$(PackageTargetDir)/%.o : %.F90
	$(F90Compiler) $(F90Options) $(IncludePath) -c $< -o $@

$(PackageParallelTargetDir)/%.o : %.F90
	$(F90Compiler) $(F90ParallelOptions) $(IncludePath) -c $< -o $@

ifeq "$(UseMPI)" "yes"
%.o : %.F90
	$(MAKE) $(MakeOptions) $(PackageTargetDir)/$@
	$(MAKE) $(MakeOptions) $(PackageParallelTargetDir)/$@
else
%.o : %.F90
	$(MAKE) $(MakeOptions) $(PackageTargetDir)/$@
endif


# .java -> .class
$(PackageTargetDir)/%.class : $(PackageSourceDir)/%.java
	$(JavaCompiler) $(JavaCompilerOptions) $< 

%.class : $(PackageSourceDir)/%.java
	@$(MAKE) $(MakeOptions) $(PackageTargetDir)/$@

# .class -> .h
$(PackageSourceDir)/%.h : $(PackageTargetDir)/%.class
	$(JniCompiler) $(JniOptions) $(JavaPackageName).$*

%.h : %.class
	$(MAKE) $(MakeOptions) $(PackageSourceDir)/$@

#.y -> .h
%.h : %.y
	$(Yacc) -d  $?
	@$(Move) y.tab.h $*.h
	@$(Delete) y.tab.c

#.y -> .c
%.c : %.y
	$(Yacc) -d  $?
	@$(Move) y.tab.c $*.c
	@$(Delete) y.tab.h

# .l -> .c
%.c : %.l
	$(Print) $@
	$(Print) $(Lex) -o$@ $?
	@$(Lex) -o$@ $?

# .o -> .a

$(LibDir)/%.a : $(ObjectFiles)
	$(StaticArchiver) $(StaticArchiverOptions) $@ $(ObjectFiles)

$(ParallelLibDir)/%.a: $(ParallelObjectFiles)
	$(StaticArchiver) $(StaticArchiverOptions) $@ $(ParallelObjectFiles)

ifeq "$(UseMPI)" "yes"
%.a : $(ObjectFiles) $(ParallelObjectFiles)
	$(MAKE) $(MakeOptions) $(LibDir)/$@
	$(MAKE) $(MakeOptions) $(ParallelLibDir)/$@
else
%.a : $(ObjectFiles)
	$(MAKE) $(MakeOptions) $(LibDir)/$@
endif

# .o -> .so
$(LibDir)/%.so : $(ObjectFiles)
	$(DynamicArchiver) $(ObjectFiles) $(DynamicArchiverOptions) -o $@

$(ParallelLibDir)/%.so : $(ParallelObjectFiles)
	$(DynamicArchiver) $(ParallelObjectFiles) $(DynamicArchiverOptions) -o $@

ifeq "$(UseMPI)" "yes"
%.so : $(ObjectFiles) $(ParallelObjectFiles)
	$(MAKE) $(MakeOptions) $(LibDir)/$@
	$(MAKE) $(MakeOptions) $(ParallelLibDir)/$@
else
%.so : $(ObjectFiles)
	$(MAKE) $(MakeOptions) $(LibDir)/$@
endif

# .class -> .jar
$(LibDir)/%.jar : $(JavaClassFiles) $(OtherTargetFiles)
	$(Print) $@
	@cd $(TargetDir); $(JavaArchiver) -cf $@ \
	$(JavaClassFilesRel) $(OtherTargetFiles)

%.jar : $(JavaClassFiles) $(OtherTargetFiles)
	$(MAKE) $(MakeOptions) $(LibDir)/$@

# .class -> JavaDoc
javadoc : 
	$(Print) $(JavaPackageNames) > $(DEV_ROOT)/packages.tmp
	$(JavadocGenerator) $(JavadocOptions) @$(DEV_ROOT)/packages.tmp
	$(Delete) $(DEV_ROOT)/packages.tmp
	$(Print) Done JavaDoc.

# .class -> _Stub.class 
$(PackageTargetDir)/%_Stub.class : $(PackageTargetDir)/%.class
	$(Print) $@
	$(RmiCompiler) $(RmiOptions) $(JavaPackageName).$*

%_Stub.class : %.class
	$(MAKE) $(MakeOptions) $(PackageTargetDir)/$@

# .class -> _Skel.class
$(PackageTargetDir)/%_Skel.class : $(PackageTargetDir)/%.class
	$(Print) $@
	$(RmiCompiler) $(RmiOptions) $(JavaPackageName).$*

%_Skel.class : %.class
	$(MAKE) $(MakeOptions) $(PackageTargetDir)/$@

$(OopseStaticLibrary) : .stamp
	touch .stamp

$(OopseStaticLibrary) : $(ObjectFiles)
	$(StaticArchiver)  $(StaticArchiverOptions) $(OopseStaticLibrary) $(ObjectFiles)

.MPIstamp : $(OopseParallelStaticLibrary)
	touch .MPIstamp

$(OopseParallelStaticLibrary) : $(ParallelObjectFiles)
	$(StaticArchiver)  $(StaticArchiverOptions) $(OopseParallelStaticLibrary) $(ParallelObjectFiles)


# Executable
$(Executable) : $(ObjectFiles)
	$(Linker) $(ObjectFiles) $(LinkOptions) $(LibDirs) $(LibList) -o $@

$(ParallelExecutable) : $(ParallelObjectFiles)
	$(Linker) $(ParallelObjectFiles)  $(ParallelLinkOptions) $(ParallelLibDirs) $(LibList) -o $@

# Anything else is just copied from source to target
$(PackageTargetDir)/% : $(PackageSourceDir)/%
	$(Print) $@
	$(Copy) $< $@

# make (or make build)
build : $(PackageListLoop)
	$(Print) Done build.

_all : _buildall

ifeq "$(UseMPI)" "yes"
_buildall : \
	$(DependencyFile) \
	$(PackageTargetDir) \
	$(ObjectFiles) \
	$(ParallelObjectFiles) \
	$(JavaClassFiles) \
	$(RmiStubFiles) \
	$(RmiSkeletonFiles) \
	$(OtherTargetFiles) \
	$(SharedLibrary) \
	$(StaticLibrary) \
	$(ParallelSharedLibrary) \
	$(ParallelStaticLibrary) \
	$(JarFile) \
	$(Executable) \
	$(ParallelExecutable) 
else
_buildall : \
	$(DependencyFile) \
	$(PackageTargetDir) \
	$(ObjectFiles) \
	$(JavaClassFiles) \
	$(RmiStubFiles) \
	$(RmiSkeletonFiles) \
	$(OtherTargetFiles) \
	$(SharedLibrary) \
	$(StaticLibrary) \
	$(JarFile) \
	$(Executable) 
endif


# make clean
clean : $(PackageListLoop)
	$(Print) Done clean.	

_cleanall :
	$(Delete) $(ObjectFiles) $(ParallelObjectFiles)

# make distclean
distclean : $(PackageListLoop)
	$(Print) Done clean.	

_distcleanall :
	$(Delete) $(ObjectFiles) \
	          $(ParallelObjectFiles) \
		        $(JarFile) \
	          $(SharedLibrary) \
	          $(StaticLibrary) \
	          $(ParallelSharedLibrary) \
	          $(ParallelStaticLibrary) \
	          $(Executable) \
	          $(ParallelExecutable) \
	          $(DependencyFile) 


# make depend
depend : $(PackageListLoop)
	$(Print) Done dependencies.

_dependall : $(DependencyFile)

$(DependencyFile) : $(DerivedSource)
	$(Print) $@ 
	@cd $(PackageSourceDir)

	@touch Make.ctemp
	@touch Make.ctemp

ifneq "$(words $(CppFiles))" "0"
	$(CppCompiler) $(IncludePath) -MM $(CppFiles)  > Make.cpptemp
	cat Make.cpptemp | sed 's/^[a-zA-Z0-9]/$$\(DEV_ROOT\)\/obj\/&/g' >> $(DependencyFile)
	cat Make.cpptemp | sed 's/^[a-zA-Z0-9]/$$\(DEV_ROOT\)\/MPIobj\/&/g' >> $(DependencyFile)
endif

ifneq "$(words $(CFiles))" "0" 
	$(CCompiler)   $(IncludePath) -MM $(CFiles) $(DerivedCFiles)   > Make.ctemp
	cat Make.ctemp | sed 's/^[a-zA-Z0-9]/$$\(DEV_ROOT\)\/obj\/&/g' >> $(DependencyFile)
	cat Make.ctemp | sed 's/^[a-zA-Z0-9]/$$\(DEV_ROOT\)\/MPIobj\/&/g' >> $(DependencyFile)

endif

ifneq "$(words $(F90Files))" "0"
	$(DEV_ROOT)/scripts/filepp  -I $(DEV_ROOT)/src -od '$$(DEV_ROOT)/obj/' -D__FORTRAN90  $(F90Files) > Make.ftemp
	$(DEV_ROOT)/scripts/filepp  -I $(DEV_ROOT)/src -od '$$(DEV_ROOT)/MPIobj/' -DIS_MPI -D__FORTRAN90  $(F90Files) >> Make.ftemp
	cat Make.ftemp >> $(DependencyFile)
endif
	$(Delete) Make.cpptemp Make.ctemp Make.ftemp

# make lib
lib    : $(PackageListLoop)
	$(Print) Libraries built.

_liball : $(JarFile) $(SharedLibrary) $(StaticLibrary)

jar : $(JarFile)

jarsign : $(JarFile)
	$(JarSigner) -keystore GeoSoftKeystore $(JarFile) myself

# make statistics
_statisticsall :
	@$(Print) $(SourceFiles) >> $(DEV_ROOT)/files.tmp

statistics : $(PackageListLoop)
	@$(List) $(DEV_ROOT)/files.tmp | xargs $(WordCount) $(WordCountOptions)
	@$(Delete) $(DEV_ROOT)/files.tmp
	$(Print) Done statistics.

# make pure
$(Executable).pure :
	$(Purify) $(PurifyOptions) $(CppCompiler) $(LinkOptions) $(LibDirs) \
	$(LibList) $(ObjectFiles) -o $@

pure : $(Executable).pure

#make cvslog
cvslog:
	$(DEV_ROOT)/scripts/cvs2cl.pl

# Execute
_runexe :
	$(Executable) $(RunParameters)

_runjava :
	$(JavaExecute) $(JavaRunOptions) $(JavaMainClass) $(RunParameters)

run : _runjava


ifdef $(DependencyFile)
-include $(DependencyFile)
endif
