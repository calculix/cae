#!/bin/bash

DIR="$(dirname $(readlink -f $0))"

if [ -z "$1" ]
then
    # Calling ./cae.sh
    exec "$DIR/src/cae"
else
    # Calling ./cae.sh -inp default.inp
    exec "$DIR/src/cae $1 $2"
fi
