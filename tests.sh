#!/bin/bash
clear
INITIAL_DIR="$(pwd)"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd $SCRIPT_DIR/tests
for f in $(find . -name \*.py)
do
    f="${f##*/}"
    python3 -m unittest -v $f
done
cd $INITIAL_DIR
