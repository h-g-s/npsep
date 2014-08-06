INCLUDES = \
 -I$(top_srcdir)/include/ \
 -I../../src/ \
 -I../../include/ 

bin_PROGRAMS = cbcsep

cbcsep_LDADD = -L/usr/lib/coin/ -L/usr/lib/coin/ThirdParty/ \
					-L/usr/local/lib/coin/ -L/usr/local/lib/coin/ThirdParty/ \
               -lm -lgfortran -lcoinlapack -lcoinblas -lrt \
					-lCoinUtils -lOsi -lClp -lOsiClp -lCgl -lCbc \
					-lglpk -lpthread

cbcsep_SOURCES = cbc_hooks.cpp \
 					main.cpp \
 					mpIndex.cpp \
 					MSepEvent.cpp \
 					GreedyCliqueExtender.cpp \
 					MIPSep.cpp \
 					mpName.cpp \
 					../../src/cgraph.c \
 					../../src/clique.c \
 					../../src/conflict_discover.c \
 					../../src/grasp.c \
 					../../src/memory.c \
 					../../src/strUtils.c \
 					../../src/vectormgm.c \
 					../../src/vint_set.c

