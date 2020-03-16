#!/bin/bash

i=1

while [ $i -lt 4 ]

do
	start=`date +%s.%N`


	echo "Running " $i " threads"
    echo "Running " $i " threads" >> partBlogOMP1600
    export OMP_NUM_THREADS=$i
	./ompCODE >>  partBlogOMP1600

	end=`date +%s.%N`

	elapsed="$(bc <<<"$end-$start")"

	echo "Runtime: " $elapsed


	i=$((i*2))

	  

done

#    start=`date +%s.%N`

# #    /opt/zeus2/bin/mpirun -np $i ./mpi_epic -itrun 100 epic.nc >& log



#    mpiexec -np $i ./partC >& log 

#    end=`date +%s`

#    elapsed="$(bc <<<"$end-$start")"

#    echo "Runtime: " $elapsed

 