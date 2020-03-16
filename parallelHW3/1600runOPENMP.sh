#!/bin/bash

i=1

while [ $i -lt 4 ]

do
	start=`date +%s.%N`


	echo "Running " $i " threads"
    echo "Running " $i " threads" >> logOMP1600
    export OMP_NUM_THREADS=$i
	./ompCODE >>  logOMP1600

	end=`date +%s.%N`

	elapsed="$(bc <<<"$end-$start")"

	echo "Runtime: " $elapsed


	if [ $i = 1 ]

		then

			i=2

		else

			i=$((i*2))

	fi

   

done

#    start=`date +%s.%N`

# #    /opt/zeus2/bin/mpirun -np $i ./mpi_epic -itrun 100 epic.nc >& log



#    mpiexec -np $i ./partC >& log 

#    end=`date +%s`

#    elapsed="$(bc <<<"$end-$start")"

#    echo "Runtime: " $elapsed

 