#!/bin/bash
cd $PBS_O_WORKDIR
qsub RESTART.pbs
date +"%m/%d/%Y %H:%M:%S $HOSTNAME" >> launched_autorestarts.txt
