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
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_CONF=Debug
CND_DISTDIR=dist

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/1360937237/clique.o \
	${OBJECTDIR}/_ext/1360937237/vectormgm.o \
	${OBJECTDIR}/_ext/1360937237/strUtils.o \
	${OBJECTDIR}/_ext/1360937237/memory.o \
	${OBJECTDIR}/_ext/1360937237/conflict_discover.o \
	${OBJECTDIR}/eclq.o \
	${OBJECTDIR}/_ext/1360937237/vint_set.o \
	${OBJECTDIR}/_ext/1360937237/cgraph.o \
	${OBJECTDIR}/_ext/1360937237/clique_enum.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-lm

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-Debug.mk dist/Debug/GNU-Linux-x86/enum_clq

dist/Debug/GNU-Linux-x86/enum_clq: ${OBJECTFILES}
	${MKDIR} -p dist/Debug/GNU-Linux-x86
	${LINK.c} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/enum_clq ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/1360937237/clique.o: ../src/clique.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -I../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/clique.o ../src/clique.c

${OBJECTDIR}/_ext/1360937237/vectormgm.o: ../src/vectormgm.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -I../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/vectormgm.o ../src/vectormgm.c

${OBJECTDIR}/_ext/1360937237/strUtils.o: ../src/strUtils.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -I../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/strUtils.o ../src/strUtils.c

${OBJECTDIR}/_ext/1360937237/memory.o: ../src/memory.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -I../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/memory.o ../src/memory.c

${OBJECTDIR}/_ext/1360937237/conflict_discover.o: ../src/conflict_discover.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -I../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/conflict_discover.o ../src/conflict_discover.c

${OBJECTDIR}/eclq.o: eclq.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.c) -g -I../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/eclq.o eclq.c

${OBJECTDIR}/_ext/1360937237/vint_set.o: ../src/vint_set.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -I../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/vint_set.o ../src/vint_set.c

${OBJECTDIR}/_ext/1360937237/cgraph.o: ../src/cgraph.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -I../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/cgraph.o ../src/cgraph.c

${OBJECTDIR}/_ext/1360937237/clique_enum.o: ../src/clique_enum.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -I../include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/clique_enum.o ../src/clique_enum.c

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} dist/Debug/GNU-Linux-x86/enum_clq

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
