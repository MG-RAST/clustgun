#!/bin/bash

# experimental..

ln -s $1 iteration1.fas
for i in {1..10}
do
	#echo iteration${i}.fas  iteration`expr $i + 1`.fas
	./clustgun iteration${i}.fas | tee iteration${i}.log
	cat iteration${i}.fas.consensus iteration1.fas  > iteration`expr $i + 1`.fas
done
