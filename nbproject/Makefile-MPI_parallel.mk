#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=mpic++
CXX=mpic++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=MPI_parallel
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/CorrelatedGaussianObservables.o \
	${OBJECTDIR}/MixingModel.o \
	${OBJECTDIR}/histo.o \
	${OBJECTDIR}/main.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=`root-config --cflags` 
CXXFLAGS=`root-config --cflags` 

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L../HEPfit/BAT_parallel/lib -Wl,-rpath,'../HEPfit/BAT_parallel/lib'

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ddbar

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ddbar: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	mpic++ -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ddbar ${OBJECTFILES} ${LDLIBSOPTIONS} -L/usr/lib -lgsl -lgslcblas -lm `root-config --libs` -lMinuit -L../HEPfit/BAT_parallel/lib -lBAT -lBATmodels

${OBJECTDIR}/CorrelatedGaussianObservables.o: CorrelatedGaussianObservables.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I../HEPfit/BAT_parallel/include -I../HEPfit/gslpp/src -I/usr/include/root -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CorrelatedGaussianObservables.o CorrelatedGaussianObservables.cpp

${OBJECTDIR}/MixingModel.o: MixingModel.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I../HEPfit/BAT_parallel/include -I../HEPfit/gslpp/src -I/usr/include/root -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MixingModel.o MixingModel.cpp

${OBJECTDIR}/histo.o: histo.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I../HEPfit/BAT_parallel/include -I../HEPfit/gslpp/src -I/usr/include/root -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/histo.o histo.cpp

${OBJECTDIR}/main.o: main.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -I../HEPfit/BAT_parallel/include -I../HEPfit/gslpp/src -I/usr/include/root -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
