#!/bin/bash

for f in ./SNAP_output_part1/*.out ./SNAP_output_part1_bu/*.out ./SNAP_output_part2/*.out ./SNAP_output_part3/*.out ./SNAP_output_part3_bu/*.out ./SNAP_output_part4/*.out ./SNAP_output_part4_bu/*.out ./SNAP_output_part5/*.out ./SNAP_output_part5_bu/*.out ./SNAP_output_part6/*.out ./SNAP_output_part6_bu/*.out ./SNAP_output_part7/*.out ./SNAP_output_part7_bu/*.out ./SNAP_output_part8/*.out
do
	python snap-scores-mutOut.py $f
done
