#!/bin/bash
rm -rf ./output
mkdir ./output

BENCH_DIR=/home/ga83dit/Documents/project/working/cookbook
CASE=tpv5
cp -r $BENCH_DIR/$CASE/* ./

BASE_DIR=/home/ga83dit/Documents/project/working/SeisSol
echo "$BASE_DIR/Maple/" > ./DGPATH
