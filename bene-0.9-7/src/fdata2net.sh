#!/bin/bash

binpath=`dirname $0`

if [ $# -ne 3 ]; then
    echo Usage: data2net.sh datafile score resultdir 1>&2
    exit 1
fi

datafile=$1;  shift
score=$1;  shift
rdir=$1;      shift

mkdir -p $rdir

nof_vars=`$binpath/L1_get_local_scores.py $datafile $score ${rdir}/res`
$binpath/split_local_scores   $nof_vars ${rdir}
$binpath/reverse_local_scores $nof_vars ${rdir}
$binpath/get_best_parents     $nof_vars ${rdir}
$binpath/get_best_sinks       $nof_vars ${rdir} ${rdir}/sinks
$binpath/get_best_order       $nof_vars ${rdir}/sinks ${rdir}/ord
$binpath/get_best_net         $nof_vars ${rdir} ${rdir}/ord ${rdir}/net
$binpath/score_net            ${rdir}/net ${rdir}
