Package = UseTheForce

Source    = \
	doForces.F90 \
	doForces_interface.F90 \
        ForceField.cpp \
        ForceFieldFactory.cpp \
        DUFF.cpp \
	EAM_FF.cpp \
	EADM_FF.cpp \
	SHAPES_FF.cpp \
	SC_FF.cpp \
	MnM_FF.cpp \
	Amber_FF.cpp \
	CLAYFF.cpp 

ParallelSource = doForces.F90

Modules = \
	doForces \
	notifyCutoffs

RmiSource =

Main      =

DEV_ROOT=../..
include $(DEV_ROOT)/make/Makefile

