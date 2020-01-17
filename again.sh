#!/bin/bash
rm -rf ./output
mkdir ./output

BENCH_DIR=/home/ga83dit/Documents/project/working/cookbook
#CASE=tpv5-giant
CASE=tpv5
cp -r $BENCH_DIR/$CASE/* ./


#TEST_DIR=/import/home/ga83dit/Documents/project/working/benchmarks/convergence_
#TEST=elastic
#cp -r $TEST_DIR$TEST/* ./


BASE_DIR=/home/ga83dit/Documents/project/working/SeisSol
echo "$BASE_DIR/Maple/" > ./DGPATH
