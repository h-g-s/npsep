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
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/1360937237/clique.o \
	${OBJECTDIR}/_ext/1360937237/clique_separation.o \
	${OBJECTDIR}/_ext/1360937237/vectormgm.o \
	${OBJECTDIR}/_ext/1360937237/strUtils.o \
	${OBJECTDIR}/_ext/1360937237/memory.o \
	${OBJECTDIR}/_ext/1360937237/lp_callbacks.o \
	${OBJECTDIR}/_ext/1360937237/conflict_discover.o \
	${OBJECTDIR}/_ext/1113355286/cbc_hooks.o \
	${OBJECTDIR}/_ext/1360937237/clique_extender.o \
	${OBJECTDIR}/_ext/1360937237/vint_set.o \
	${OBJECTDIR}/_ext/1360937237/grasp.o \
	${OBJECTDIR}/_ext/1360937237/cgraph.o \
	${OBJECTDIR}/cbcroot.o \
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
LDLIBSOPTIONS=-L/usr/lib/coin -L/usr/lib/coin/ThirdParty -lm -lgfortran -lcoinlapack -lcoinblas -lrt -lCoinUtils -lOsi -lClp -lOsiClp -lCgl -lCbc -lglpk -lpthread

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/cbcr

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/cbcr: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/cbcr ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/1360937237/clique.o: ../src/clique.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DDEBUG -I../include -I../cbc/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/clique.o ../src/clique.c

${OBJECTDIR}/_ext/1360937237/clique_separation.o: ../src/clique_separation.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DDEBUG -I../include -I../cbc/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/clique_separation.o ../src/clique_separation.c

${OBJECTDIR}/_ext/1360937237/vectormgm.o: ../src/vectormgm.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DDEBUG -I../include -I../cbc/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/vectormgm.o ../src/vectormgm.c

${OBJECTDIR}/_ext/1360937237/strUtils.o: ../src/strUtils.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DDEBUG -I../include -I../cbc/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/strUtils.o ../src/strUtils.c

${OBJECTDIR}/_ext/1360937237/memory.o: ../src/memory.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DDEBUG -I../include -I../cbc/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/memory.o ../src/memory.c

${OBJECTDIR}/_ext/1360937237/lp_callbacks.o: ../src/lp_callbacks.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DDEBUG -I../include -I../cbc/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/lp_callbacks.o ../src/lp_callbacks.c

${OBJECTDIR}/_ext/1360937237/conflict_discover.o: ../src/conflict_discover.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DDEBUG -I../include -I../cbc/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/conflict_discover.o ../src/conflict_discover.c

${OBJECTDIR}/_ext/1113355286/cbc_hooks.o: ../cbc/src/cbc_hooks.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1113355286
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I../include -I../cbc/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1113355286/cbc_hooks.o ../cbc/src/cbc_hooks.cpp

${OBJECTDIR}/_ext/1360937237/clique_extender.o: ../src/clique_extender.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DDEBUG -I../include -I../cbc/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/clique_extender.o ../src/clique_extender.c

${OBJECTDIR}/_ext/1360937237/vint_set.o: ../src/vint_set.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DDEBUG -I../include -I../cbc/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/vint_set.o ../src/vint_set.c

${OBJECTDIR}/_ext/1360937237/grasp.o: ../src/grasp.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DDEBUG -I../include -I../cbc/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/grasp.o ../src/grasp.c

${OBJECTDIR}/_ext/1360937237/cgraph.o: ../src/cgraph.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DDEBUG -I../include -I../cbc/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/cgraph.o ../src/cgraph.c

${OBJECTDIR}/cbcroot.o: cbcroot.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -DDEBUG -I../include -I../cbc/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/cbcroot.o cbcroot.cpp

${OBJECTDIR}/_ext/1360937237/clique_enum.o: ../src/clique_enum.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1360937237
	${RM} $@.d
	$(COMPILE.c) -g -Wall -DDEBUG -I../include -I../cbc/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1360937237/clique_enum.o ../src/clique_enum.c

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/cbcr

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
