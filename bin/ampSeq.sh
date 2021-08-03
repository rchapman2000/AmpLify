#!/usr/bin/env bash

help() {
    echo ""
    echo " AmpSeq Pipeline"
    echo "Usage: ampseq <module> [module options]"
    echo ""
    echo "Modules:"
    echo ""
    echo "  - preProcess"
    echo "  - calcDepth"
    echo "  - generateConsensus"
    echo ""
    echo "To view this message:"
    echo "  ampseq -h"
    echo ""
    echo "To view module help messages:"
    echo "  ampseq <module> -h"
    echo "";}

SCRIPT_DIR=$(dirname ${BASH_SOURCE})

if [ "$1" = preProcess ]; then
    echo $( dirname ${BASH_SOURCE})
    python3 "$SCRIPT_DIR"/preProcessing.py ${@:2}
elif [ "$1" = calcDepth ]; then
    python3 "$SCRIPT_DIR"/calcDepth.py ${@:2}
elif [ "$1" = generateConsensus ]; then
    python3 "$SCRIPT_DIR"/generateConsensus.py ${@:2}
elif [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    help
else 
    echo "Module "$1" does not exist"
    help
fi