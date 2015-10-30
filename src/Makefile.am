INCLUDES = \
	-I$(top_srcdir)/include/ -I$(top_srcdir)/CglEClique/  -I/usr/include/coin/

AM_CFLAGS= -fopenmp
AM_CXXFLAGS= -fopenmp

bin_PROGRAMS=cbcr writeclqw cptest


cptest_SOURCES= cgraph.c\
				vint_set.c\
				memory.c\
				clique.c\
				vectormgm.c\
				strUtils.c\
				constraint_propagation.cpp\
				osi_cgraph.cpp\
				preprocess.cpp\
				problem.cpp\
				cp_test.cpp

writeclqw_SOURCES=clique.c\
	memory.c\
	strUtils.c\
	vectormgm.c\
	vint_set.c\
	writeclqw.cpp\
	osi_cgraph.cpp \
	cgraph.c \
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
				 osi_cgraph.cpp \
				 cbcroot.cpp \
				 spaths.c \
				 oddhs.c \
				 ../CglEClique/CglEClique.cpp \
				 ../CglEClique/CglEClique.hpp \
				 clique_separation.h \
				 cgraph.h \
				 clique_extender.h \
				 oddhs.h \
				 clique.h \
				 osi_cgraph.h	\
				 bron_kerbosch.h \
				 spaths.h \
				 memory.h \
				 strUtils.h \
				 clique_enum.h \
				 node_heap.h \
				 vectormgm.h \
				 vint_set.h \
				 vint_queue.h \
				 grasp.h \
				 BKGraph.hpp \
				 BKVertex.hpp \
				 macros.h


#bk_SOURCES=bkRun.cpp\
#			  cgraph.c \
#			  bron_kerbosch.cpp\
#			  vint_set.c\
#			  memory.c\
#			  clique.c\
#			  strUtils.c\
#			  vectormgm.c\
#			  BKVertex.cpp\
#			  BKGraph.cpp

cbcr_LDADD=-lz -ldl -lrt -lm  ${cbc_LIBS}

writeclqw_LDADD=-lz -ldl -lrt -lm  ${cbc_LIBS}

cptest_LDADD=-lz -ldl -lrt -lm  ${cbc_LIBS}
