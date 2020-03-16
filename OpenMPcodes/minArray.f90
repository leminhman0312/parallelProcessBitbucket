program findMinArr

    implicit none
    
    !! DECLARE BASIC VARIABLES FOR THIS PROGRAM 
    INTEGER,PARAMETER :: MAX = 1e7

    DOUBLE PRECISION,dimension(MAX) :: x
    DOUBLE PRECISION,DIMENSION(10) :: my_min
    DOUBLE PRECISION :: rmin


    INTEGER :: num_threads
    INTEGER :: i, n 
    INTEGER :: id,start, stop

    !! DECLARE OPENMP FUNCTIONS
    INTEGER,EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS



    !$OMP PARALLEL PRIVATE(i,id,start,stop,num_threads,n)

    num_threads = OMP_GET_NUM_THREADS()
    n = MAX/num_threads

    !Find start
    start = id*n + 1

    !Find stop 
    if (id <> (num_threads-1)) then 
        stop = start + n 
    else
        stop = MAX
    end if

    my_min(id+1) = x(start)

    do i = start+1,stop
        if (x(i) < my_min(id+1)) then 
            my_min(id+1) = x(i)
        end if
    end do

    !$OMP END PARALLEL



    rmin = my_min(1)

    DO i = 2, num_threads
        IF ( rmin < my_min(i) ) THEN
        rmin = my_min(i)
        END IF
    END DO


    print *, "min = ", rmin

end program findMinArr