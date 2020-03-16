PROGRAM findmax
    use omp_lib 
    IMPLICIT NONE
    INTEGER,DIMENSION(10) :: x_arr 
    INTEGER :: xmax,i
    INTEGER :: num_threads 
    DOUBLE PRECISION :: tstart, tend, t_elapse 


    num_threads = 16
    call OMP_SET_NUM_THREADS(num_threads)


    x_arr = (/1,4,9,2,15,20,4,34,45,5/)


    !! FIND MAX 

    xmax = 0

    tstart = OMP_GET_WTIME()
    !$OMP PARALLEL DO REDUCTION(max:xmax)
    DO i = 1,10
        IF(x_arr(i) > xmax) THEN 
            xmax = x_arr(i)
        END IF 
    END DO 
    !$OMP END PARALLEL DO 

    tend = OMP_GET_WTIME()
    t_elapse = tend - tstart 


    WRITE(*,*) "MAX = ", xmax
    WRITE(*,*) ""
    WRITE(*,*) "Time taken :",t_elapse  

END PROGRAM findmax 
