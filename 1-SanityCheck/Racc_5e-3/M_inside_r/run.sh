#!/bin/bash
for res in '100k' '250k' '500k' '750k' '001M';
do
	./particles2bins.x /scratch/jcarmona/Racc_5e-3_res_${res}/Output/snapshot_000 ./res_${res}_t000.txt
done
