#!/bin/bash
exec $(dirname $(readlink -f $0))/src/cae.py $*
