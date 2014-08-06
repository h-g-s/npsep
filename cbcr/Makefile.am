INCLUDES = \
 -I$(top_srcdir)/include/ \
 -I../include/ \
 -I$(top_srcdir)/src/ \
 -I$(top_srcdir)/cbc/include \
 -I../cbc/include \
 -I$(top_srcdir)/cbc/src

bin_PROGRAMS = cbcr

cbcr_LDADD = -lm -lpthread -lrt -ldl -lcoinglpk -lcoinblas \
						-lcoinlapack -lgfortran -lCoinUtils -lOsi -lClp \
						-lOsiClp -lCgl -lCbc

cbcr_SOURCES =  cbcroot.cpp \
							  ../cbc/src/cbc_hooks.cpp \
							  ../src/vnode_degree.c \
							  ../src/cgraph.c \
							  ../src/clique.c \
							  ../src/clique_enum.c \
							  ../src/clique_extender.c \
							  ../src/conflict_discover.c \
							  ../src/grasp.c \
							  ../src/lp_callbacks.c \
							  ../src/node_heap.c \
							  ../src/clique_separation.c \
							  ../src/memory.c \
							  ../src/strUtils.c \
							  ../src/vectormgm.c \
							  ../src/vint_set.c 

