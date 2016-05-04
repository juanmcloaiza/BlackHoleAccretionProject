#!/bin/bash
for snap in '010' '050' '100' '200' '300';
do
	./particles2bins.x /scratch/jcarmona/No_racc_res_001M/Output/snapshot_${snap} ./res_001M_t${snap}.txt
done
