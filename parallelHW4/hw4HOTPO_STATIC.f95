PROGRAM hw4
    USE omp_lib 
    IMPLICIT NONE
    INTEGER*8:: n, iterHOTPO,iterNMAX, highN, MAX_HOTPO       
    INTEGER :: num_threads 
    DOUBLE PRECISION :: tstart, tend, t_elapse
    INTEGER*8 :: GLOBAL_MAX 
    INTEGER :: nchunk 

    MAX_HOTPO = 1e8
    GLOBAL_MAX = 0

    num_threads = 4
    call OMP_SET_NUM_THREADS(num_threads)

    nchunk = 1

    tstart = OMP_GET_WTIME()

    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(PRIVATE)&
    !$OMP& SHARED(MAX_HOTPO) REDUCTION(max:GLOBAL_MAX)
    DO iterNMAX = 2,20000
        n = iterNMAX   
        highN = 0 

        DO iterHOTPO = 1,MAX_HOTPO

            !! FIND MAX N  
            IF (n > highN) THEN 
                highN = n
            END IF


            CALL hotpo(n)

            !! BREAK AT 1 
            IF (n == 1) THEN
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
    WRITE(*,*) "STATIC, NCHUNK = ", nchunk
    WRITE(*,*) "RUNNING WITH ", num_threads, " CORES"
    WRITE(*,*) "Time taken :",t_elapse  
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
