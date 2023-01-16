#!/bin/bash
clear
INITIAL_DIR="$(pwd)"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd $SCRIPT_DIR/src/test
for f in $(find . -name \*.py)
do
    echo ''
    echo ''
    echo ''
    f="${f##*/}"
    $SCRIPT_DIR/bin/python/bin/python3 -m unittest -v $f
    echo ''
    echo ''
    echo ''
done
cd $INITIAL_DIR
read
