Source    = \
	definitions.F90 \
	fError.c \
	PropertyMap.cpp \
	oopseMPI.F90 \
	simError.c \
	status.F90 \
	StringUtils.cpp \
	vector_class.F90 \
        StringTokenizer.cpp \
	Utility.cpp \
	Trim.cpp

ifneq "$(words $(LIBOBJS))" "0"
   Source  = $(Source) getopt.c getopt1.c
endif

RmiSource =

Main      =

DEV_ROOT=../..
include $(DEV_ROOT)/make/Makefile
