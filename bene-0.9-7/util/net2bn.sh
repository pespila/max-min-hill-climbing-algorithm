#!/bin/bash

if [ $# -ne 1 ]; then
    echo Usage: $0 netfile
    exit 1
fi
td=`dirname $0`

wc -l < $1
$td/../src/net2parents $1 - | $td/../src/parents2arcs - -
