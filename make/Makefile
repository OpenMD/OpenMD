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
#---------------------------------------------------------------------------
#
#  Directories
#
#---------------------------------------------------------------------------
SourceDir        = $(DEV_ROOT)/src
TargetDir        = $(DEV_ROOT)/obj
LibDir           = $(DEV_ROOT)/lib
MakeDir          = $(DEV_ROOT)/make
BinDir           = $(DEV_ROOT)/bin
DocsDir          = $(DEV_ROOT)/docs
CurrentDir       = $(CURDIR)

ifdef Source
Package          = $(subst $(SourceDir)/,,$(CurrentDir))
PackageList      = $(Package)
PackageSourceDir = $(SourceDir)/$(Package)
PackageTargetDir = $(TargetDir)
JavaMainClass    = $(subst /,.,$(Package)).$(Main)
else
PackageList      = $(Packages) $(JavaPackages)
endif

PackageListLoop  = $(patsubst %,$(SourceDir)/%/.loop,$(PackageList))

JRE              = $(JAVA_HOME)/jre/lib/rt.jar

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
SourceFiles            = $(JavaFiles:%.java=  $(PackageSourceDir)/%.java)\
                         $(CppFiles:%.cpp=    $(PackageSourceDir)/%.cpp)\
                         $(CFiles:%.c=        $(PackageSourceDir)/%.c)\
                         $(FortranFiles:%.f=  $(PackageSourceDir)/%.f)\
                         $(F90Files:%.F90=    $(PackageSourceDir)/%.F90)\
                         $(LexFiles:%.l=      $(PackageSourceDir)/%.l)\
                         $(YaccFiles:%.y=     $(PackageSourceDir)/%.y)\


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
ifneq  "$(words $(ObjectFiles))" "0"
  DependencyFile    = $(PackageSourceDir)/Makedepend
  SharedLibrary     = $(LibDir)/lib$(subst /,,$(Package)).so 
  StaticLibrary     = $(LibDir)/lib$(subst /,,$(Package)).a
  ifneq "$(Main)" ""
    Executable        = $(BinDir)/$(Main)
  endif
endif
#
# Misc
#
ClassPath        = $(JRE)$(X)$(TargetDir)$(X)$(ThirdPartyJars)
JavaPackageNames = $(subst /,.,$(JavaPackages))
IncludePath      = -I$(SourceDir) $(IncludeDirs:%=-I%)
LibDirs          = -L$(LibDir)    $(LibraryDirs:%=-L%)
LocalLibs        = $(subst /,,$(Packages))
LibList          = $(LocalLibs:%=-l%) $(Libraries:%=-l%)


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
MakeDir                = mkdir -p
Delete                 = rm -fr
StaticArchiver         = ar
DynamicArchiver        = icc
FortranCompiler        = ifort
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


MakeOptions            = -k -s
MakeDependOptions      = 
StaticArchiverOptions  = rc
DynamicArchiverOptions = -shared
JavaArchiverOptions    =
JniOptions             =
RmiOptions             = -d $(TargetDir) -classpath $(ClassPath) \
			 -sourcepath $(SourceDir)
FortranOptions         =
F90Options             =  -I$(SourceDir) -module $(TargetDir) 
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
$(PackageTargetDir)/%.o : $(PackageSourceDir)/%.c
	$(Print) $@
	@$(CCompiler) $(COptions) -c $(IncludePath) $< -o $@

%.o : $(PackageSourceDir)/%.c
	$(MAKE) $(MakeOptions) $(PackageTargetDir)/$@

# .cpp -> .o
$(PackageTargetDir)/%.o : $(PackageSourceDir)/%.cpp
	$(Print) $@
	$(CppCompiler) $(CppOptions) -c $(IncludePath) $< -o $@

%.o : $(PackageSourceDir)/%.cpp
	$(MAKE) $(MakeOptions) $(PackageTargetDir)/$@

# .f -> .o
$(PackageTargetDir)/%.o : $(PackageSourceDir)/%.f
	$(Print) $@
	@$(FortranCompiler) $(FortranOptions) -c $< -o $@

%.o : $(PackageSourceDir)/%.f
	$(MAKE) $(MakeOptions) $(PackageTargetDir)/$@

# .F90 -> .o
$(PackageTargetDir)/%.o : $(PackageSourceDir)/%.F90
	$(Print) $@
	@$(F90Compiler) $(F90Options) -c $< -o $@

%.o : $(PackageSourceDir)/%.F90
	$(MAKE) $(MakeOptions) $(PackageTargetDir)/$@

# .java -> .class
$(PackageTargetDir)/%.class : $(PackageSourceDir)/%.java
	$(Print) $@
	@$(JavaCompiler) $(JavaCompilerOptions) $< 

%.class : $(PackageSourceDir)/%.java
	@$(MAKE) $(MakeOptions) $(PackageTargetDir)/$@

# .class -> .h
$(PackageSourceDir)/%.h : $(PackageTargetDir)/%.class
	$(Print) $@
	$(JniCompiler) $(JniOptions) $(JavaPackageName).$*

%.h : %.class
	$(MAKE) $(MakeOptions) $(PackageSourceDir)/$@

#.y -> .h
$(PackageSourceDir)/%.h : $(PackageSourceDir)/%.y
	$(Print) $@
	$(Print) @$(Yacc) -d $?
	@$(Yacc) -d  $?
	@$(Move) y.tab.h $*.h
	@$(Delete) y.tab.c

#.y -> .c
$(PackageSourceDir)/%.c : $(PackageSourceDir)/%.y
	$(Print) $@
	$(Print) @$(Yacc) -d $?
	@$(Yacc) -d  $?
	@$(Move) y.tab.c $*.c
	@$(Delete) y.tab.h

# .l -> .c
$(PackageSourceDir)/%.c : $(PackageSourceDir)/%.l
	$(Print) $@
	$(Print) @$(Lex) -o$@ $?
	@$(Lex) -o$@ $?

# .o -> .a
$(LibDir)/%.a : $(ObjectFiles)
	$(Print) $@
	@$(StaticArchiver) $(StaticArchiverOptions) $@ $(ObjectFiles)

%.a : $(ObjectFiles)
	$(MAKE) $(MakeOptions) $(LibDir)/$@

# .o -> .so
$(LibDir)/%.so : $(ObjectFiles)
	$(Print) $@	
	$(DynamicArchiver) $(ObjectFiles) $(DynamicArchiverOptions) -o $@

%.so : $(ObjectFiles)
	$(MAKE) $(MakeOptions) $(LibDir)/$@

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

# Executable
$(Executable) : $(ObjectFiles)
	$(Print) $@
	$(Linker) $(LinkOptions) $(LibDirs) $(LibList) $(ObjectFiles) -o $@

# Anything else is just copied from source to target
$(PackageTargetDir)/% : $(PackageSourceDir)/%
	$(Print) $@
	$(Copy) $< $@

# make (or make build)
build : $(PackageListLoop)
	$(Print) Done build.

_all : _buildall

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


# make clean
clean : $(PackageListLoop)
	$(Print) Done clean.	

_cleanall :
	$(Delete) $(PackageTargetDir)/* \
		  $(JarFile) \
	          $(SharedLibrary) \
	          $(StaticLibrary) \
	          $(Executable) \
	          $(DependencyFile)


# make depend
depend : $(PackageListLoop)
	$(Print) Done dependencies.

_dependall : $(DependencyFile)

$(DependencyFile) : $(DerivedSource)
	$(Print) $@ 
	$(Print) $(DerivedSource)
	@cd $(PackageSourceDir)
	$(CppCompiler) $(IncludePath) -MM $(CppFiles)  > Make.temp
	$(CCompiler)   $(IncludePath) -MM $(CFiles) $(DerivedCFiles)   >> Make.temp
	cat Make.temp | sed 's/^[a-zA-Z0-9]/$$\(DEV_ROOT\)\/obj\/&/g' >> $(DependencyFile)
	cat Make.temp | sed 's/^[a-zA-Z0-9]/$$\(DEV_ROOT\)\/MPIobj\/&/g' >> $(DependencyFile)
	#$(DEV_ROOT)/scripts/sfmakedepend -I $(DEV_ROOT)/src -d $(DEV_ROOT)/obj -f ./Make.ftemp -h *.F90
	#$(DEV_ROOT)/scripts/sfmakedepend -I $(DEV_ROOT)/src -d $(DEV_ROOT)/MPIobj -f ./Make.ftemp -h *.F90
	#cat Make.ftemp >> $(DependencyFile)
	$(Delete) Make.temp Make.ftemp

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

# Execute
_runexe :
	$(Executable) $(RunParameters)

_runjava :
	$(JavaExecute) $(JavaRunOptions) $(JavaMainClass) $(RunParameters)

run : _runjava


ifdef $(DependencyFile)
-include $(DependencyFile)
endif
