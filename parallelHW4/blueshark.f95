PROGRAM hw4
    USE omp_lib 
    IMPLICIT NONE
    INTEGER*8:: n, iterHOTPO,iterlast, NMAX , iterNMAX, highN, MAX_HOTPO, n_start, NINIT     
    INTEGER :: num_threads 
    DOUBLE PRECISION :: tstart, tend, t_elapse
    INTEGER*8,DIMENSION(:),ALLOCATABLE :: iter_arr, max_arr, n_arr  
    INTEGER*8 :: GLOBAL_MAX 

    NMAX = 3 
    MAX_HOTPO = 1e8
    GLOBAL_MAX = 0
    NINIT = 2 

    ALLOCATE(iter_arr(NMAX-1))
    ALLOCATE(max_arr(NMAX-1)) 
    ALLOCATE(n_arr(NMAX-1)) 

    num_threads = 4
    call OMP_SET_NUM_THREADS(num_threads)


    tstart = OMP_GET_WTIME()

    !$OMP PARALLEL DO SCHEDULE(STATIC,1) DEFAULT(PRIVATE)&
    !$OMP& SHARED(MAX_HOTPO, iter_arr, max_arr,n_arr) REDUCTION(max:GLOBAL_MAX)
    DO iterNMAX = 2,3 
        n = iterNMAX   
        highN = 0 
        n_start = n 

        DO iterHOTPO = 1,MAX_HOTPO

            !! FIND MAX N  
            IF (n > highN) THEN 
                highN = n
            END IF


            CALL hotpo(n)
            !WRITE(*,*) "ni =        ",  n


            !! BREAK AT 1 
            IF (n == 1) THEN
                !iterlast = iterHOTPO
                !iter_arr(iterNMAX-1) = iterlast 
                !max_arr(iterNMAX-1)  = highN
                !n_arr(iterNMAX-1) = n_start 
                EXIT 
            END IF 
        END DO
        

        !! CALCULATE GLOBAL MAX 
        IF (highN > GLOBAL_MAX) THEN 
            GLOBAL_MAX = highN
        END IF



    END DO
    !$OMP END PARALLEL DO


    tend = OMP_GET_WTIME() 





    t_elapse = tend - tstart 


    WRITE(*,*) ""
    WRITE(*,*) "RUNNING WITH ", num_threads, " CORES"
    WRITE(*,*) "PROBLEM SIZE ", NMAX 
    WRITE(*,*) "Time taken :",t_elapse  




    !DO iterNMAX = 2,NMAX

        !WRITE(*,*) n_arr(iterNMAX-1), iter_arr(iterNMAX-1), max_arr(iterNMAX-1) 
    
    !END DO 


    WRITE(*,*) "GLOBAL MAX = ", GLOBAL_MAX

END PROGRAM hw4



SUBROUTINE hotpo(ndum)
    IMPLICIT NONE
    INTEGER*8, INTENT(INOUT) :: ndum

    IF (mod(ndum,2) == 0) THEN
        ndum = ndum / 2 
    ELSE
        ndum = (3 * ndum ) + 1
    END IF

END SUBROUTINE hotpo 
