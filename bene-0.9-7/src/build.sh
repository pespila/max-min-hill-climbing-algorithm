#!/bin/bash

if [ ${1:-UNIX} == "WIN" ]; then
    CC="i586-mingw32msvc-gcc"
    EXT=".exe"
else
    CC=gcc
fi


# CFLAGS="-Wextra -g -ansi -pedantic"
# CFLAGS="-Wall -O3 -g -pg"
CFLAGS="-Wall  -O3"

$CC $CFLAGS -c -o files.o files.c
$CC $CFLAGS -c -o varpar.o varpar.c
$CC $CFLAGS -c -o gopt.o gopt.c

$CC $CFLAGS -o get_local_scores$EXT files.c reg.c ilogi.c ls_XIC.c ls_NML.c ls_BDe.c ls_LOO.c get_local_scores.c gopt.o -lm
$CC $CFLAGS -o split_local_scores$EXT split_local_scores.c files.o
$CC $CFLAGS -o reverse_local_scores$EXT reverse_local_scores.c files.o
$CC $CFLAGS -o get_best_parents$EXT get_best_parents.c files.o
$CC $CFLAGS -o get_best_sinks$EXT get_best_sinks.c files.o
$CC $CFLAGS -o get_best_order$EXT get_best_order.c
$CC $CFLAGS -o get_best_net$EXT get_best_net.c files.o varpar.o
$CC $CFLAGS -o score_net$EXT score_net.c files.o varpar.o
$CC $CFLAGS -o score_nets$EXT score_nets.c files.o varpar.o
$CC $CFLAGS -o net2parents$EXT net2parents.c 
$CC $CFLAGS -o parents2arcs$EXT parents2arcs.c
$CC $CFLAGS -o arcs2dot$EXT arcs2dot.c
$CC $CFLAGS -o ban_arc$EXT ban_arc.c files.o
$CC $CFLAGS -o force_arc$EXT force_arc.c files.o
$CC $CFLAGS -o gen_prior_file$EXT gen_prior_file.c gopt.o -lm
