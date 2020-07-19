#!/bin/bash

DIR=$(dirname $(readlink -f $0))
CAE=$DIR/src/cae

if [ -f $CAE ]; then
    echo "Running binary."
else
    echo "Binary does not exist. Running source code."
    CAE=$CAE.py
fi

exec $CAE $*
