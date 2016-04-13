#!/bin/bash
for res in '100k' '250k' '500k' '750k' '001M';
do
	./particles2bins.x /scratch/jcarmona/No_racc_res_${res}/Output/snapshot_500 ./res_${res}_t500.txt
done
