#!/bin/bash

i=1

while [ $i -lt 4 ]

do
	start=`date +%s.%N`


	echo "Running " $i " procs"

    echo "Running " $i " procs" >> logMPI16


	mpirun -n $i ./partC >> logMPI16

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

 