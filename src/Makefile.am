INCLUDES = \
	-I$(top_srcdir)/include/ -I$(top_srcdir)/CglEClique/  -I/usr/include/coin/

#grasp_CFLAGS=-fpic -static
#grasp_LDFLAGS=-static

# Possible directives:
#    -DDEBUG_CONF -DCONSIDER_BINARY_COMPLEMENT

AM_CFLAGS= -fopenmp
AM_CXXFLAGS= -fopenmp

bin_PROGRAMS=cbcr writeclqw cptest
#bin_PROGRAMS=grasp cbcr writeclqw bk


cptest_SOURCES= cgraph.c\
				vint_set.c\
				memory.c\
				clique.c\
				vectormgm.c\
				strUtils.c\
				constraint_propagation.cpp\
				osi_cgraph.cpp\
				cp_test.cpp


#grasp_SOURCES=cgraph.c\
#	clique.c\
#	conflict_discover.c\
#	grasp.c\
#	memory.c\
#	strUtils.c\
#	vectormgm.c\
#	clique_extender.c\
#	lp_callbacks.c\
#	grasp_run.c\
#	vint_set.c 

writeclqw_SOURCES=clique.c\
	conflict_discover.c\
	lp_callbacks.c\
	memory.c\
	strUtils.c\
	vectormgm.c\
	vint_set.c\
	vnode_degree.c\
	writeclqw.cpp\
	cbc_hooks.cpp\
	osi_cgraph.cpp \
	node_heap.c \
	cgraph.c \
	constraint_propagation.cpp
	

#dijkstra_SOURCES=dijkstra.c\
#					  spaths.c\
#					  node_heap.c\
#					  memory.c\
#					  strUtils.c

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
				 conflict_discover.c \
				 strUtils.c \
				 vnode_degree.c \
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
				 conflict_discover.h \
				 cpx_cgraph.h \
				 clique_extender.h \
				 digraph.h \
				 oddhs.h \
				 clique.h \
				 osi_cgraph.h	\
				 bron_kerbosch.h \
				 lp_callbacks.h \
				 spaths.h \
				 memory.h \
				 strUtils.h \
				 clique_enum.h \
				 node_heap.h \
				 vectormgm.h \
				 vint_set.h \
				 vint_queue.h \
				 grasp.h \
				 vnode_degree.h \
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

#writeclqw_LDADD=-lOsiClp -lOsi -lClp -lCoinUtils -lcoinmumps -lcoinlapack -lcoinblas -lgfortran -lz -ldl -lrt -lm
#writeclqw_LDADD=-lz -ldl -lrt -lm `pkg-config --libs cbc`

