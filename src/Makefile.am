INCLUDES = \
	-I/usr/include/coin/

AM_CFLAGS= -fopenmp
AM_CXXFLAGS= -fopenmp

bin_PROGRAMS=cbcr writeclqw cptest bk


cptest_SOURCES= cgraph.c\
				vint_set.c\
				memory.c\
				clique.c\
				vectormgm.c\
				strUtils.c\
				constraint_propagation.cpp\
				problem.cpp\
				build_cgraph.cpp\
				preprocess.cpp\
				cp_test.cpp

writeclqw_SOURCES=clique.c\
	memory.c\
	strUtils.c\
	vectormgm.c\
	vint_set.c\
	cgraph.c \
	writeclqw.cpp\
	problem.cpp\
	build_cgraph.cpp \
	constraint_propagation.cpp

cbcr_SOURCES=cgraph.c \
				 clique_extender.c \
				 grasp.c \
				 memory.c \
				 vectormgm.c \
				 clique.c \
				 clique_separation.c \
				 node_heap.c \
				 vint_set.c \
				 vint_queue.c \
				 clique_enum.c \
				 strUtils.c \
				 BKGraph.cpp \
				 BKVertex.cpp \
				 bron_kerbosch.cpp \
				 problem.cpp\
				 build_cgraph.cpp \
				 cbcroot.cpp \
				 CglEClique.cpp \
				 spaths.c \
				 oddhs.c \
				 macros.h


bk_SOURCES=bkRun.cpp\
			  cgraph.c \
			  bron_kerbosch.cpp\
			  vint_set.c\
			  memory.c\
			  clique.c\
			  strUtils.c\
			  vectormgm.c\
			  BKVertex.cpp\
			  BKGraph.cpp

cbcr_LDADD=-lz -ldl -lrt -lm  ${cbc_LIBS}

writeclqw_LDADD=-lz -ldl -lrt -lm  ${cbc_LIBS}

cptest_LDADD=-lz -ldl -lrt -lm  ${cbc_LIBS}

bk=-lz -ldl -lrt -lm  ${cbc_LIBS}
