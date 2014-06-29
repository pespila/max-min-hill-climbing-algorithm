#!/bin/bash

d=`dirname $0`

if [ $# -eq 1 ]; then
	ess=$1
else
    ess=1
fi

score=`$d/../src/data2net.sh $d/iris.vd $d/iris.idt $ess $d/resdir`
echo Score : $score
echo Arcs : 
$d/../src/net2parents $d/resdir/net - | $d/../src/parents2arcs - -

if which dot > /dev/null; then
    $d/../src/net2parents $d/resdir/net - \
	| $d/../src/parents2arcs - - \
	| $d/../src/arcs2dot $d/iris.vd - - \
	| dot -Tps -o $d/resdir/iris.ps
    echo See $d/resdir/iris.ps for a postscript picture of the net.
fi
