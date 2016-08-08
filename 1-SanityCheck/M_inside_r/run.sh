#!/bin/bash
#for res in '100k' '250k' '500k' '750k' '001M' '002M' '004M' '008M' '010M';
#do
	#./particles2bins.x /scratch/jcarmona/Racc_1e-4_res_${res}/Output/snapshot_000 ./res_${res}_t000.txt
#	ls /scratch/jcarmona/Racc_1e-4_res_${res}/Output/snapshot_* | tail -n 1
#done

for res in '001M' #'002M' '004M' '008M' '010M';
do
	./particles2bins.x /scratch/jcarmona/Racc_1e-4_res_${res}/Output/snapshot_500 ./res_${res}_t500.txt
done
