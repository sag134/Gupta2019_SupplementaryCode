#!/bin/bash
for i in $(seq 1 1 10)
do
	 echo $i
	 screen -dmS PA$i bash -c "sh ../run_matlab_with_lasso.sh;exec bash"
	 sleep 30
done
