#!/bin/bash

for i in $(seq 1 1 4)
do
	 echo $i
	 screen -dmS cont$i bash -c "sh ../run_matlab_with_lasso.sh;exec bash"
	 sleep 15
done
