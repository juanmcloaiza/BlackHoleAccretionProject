#!/bin/bash

for RACC in '1e-3' '5e-3' '1e-4'
do
	for RES in '250k' #'100k' '500k' '750k' '001M' '004M'
	do
		cp /scratch/jcarmona/Racc_${RACC}_res_${RES}/Output/Accretion.txt ./Accretion/Racc_${RACC}_res_${RES}.txt
	done
done
