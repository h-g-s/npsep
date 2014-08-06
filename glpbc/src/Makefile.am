INCLUDES = \
 -I$(top_srcdir)/include/ \
 -I../../src/ \
 -I../../include/ 

bin_PROGRAMS = glpbc

glpbc_LDADD = -lm -lglpk

glpbc_SOURCES =  glp_hooks.c \
 					  glpbc.c \
 					  ../../src/cgraph.c \
 					  ../../src/clique.c \
 					  ../../src/conflict_discover.c \
 					  ../../src/grasp.c \
 					  ../../src/memory.c \
 					  ../../src/strUtils.c \
 					  ../../src/vectormgm.c \
 					  ../../src/vint_set.c \
 					  ../../src/clique_extender.c \
					  ../../src/lp_callbacks.c

