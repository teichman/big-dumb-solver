#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: ./solve.sh YAML_CONFIG"
    exit 0
fi

CONFIG=$1
NAME=$(basename ${CONFIG%.yaml})
CPP=generated/$NAME.cpp
EXE=bin/$NAME

mkdir -p bin
mkdir -p generated

./codegen.py $CONFIG -o $CPP
echo "Compiling"
g++ -std=c++14 -O3 -I $(brew --prefix eigen)/include/eigen3/ -I. $CPP bds.cpp -o $EXE
echo "Optimizing"
./$EXE
