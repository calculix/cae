#!/bin/bash
# exec $(dirname $(readlink -f $0))/src/cae.py $*
# python3 $(dirname $(readlink -f $0))/src/cae.py $*
$(dirname $(readlink -f $0))/bin/python/bin/python3 $(dirname $(readlink -f $0))/src/cae.py $*
read