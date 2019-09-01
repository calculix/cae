#!/bin/bash

DIR="$(dirname $(readlink -f $0))"

if [ -z "$1" ]
then
    exec "$DIR/src/cae"
else
    # ./cae.sh -inp default.inp
    # ./cae.sh -inp=default.inp
    exec "$DIR/src/cae $1 $2"
fi
