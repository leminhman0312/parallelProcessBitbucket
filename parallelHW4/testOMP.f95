PROGRAM test

    use omp_lib 
    IMPLICIT NONE
    INTEGER*8 :: n,iterNMAX,j,highN, GLOBAL_MAX 
    INTEGER*8 :: maxHOTPO 

    INTEGER*8 :: nstart, nend 

    INTEGER :: num_threads 
    DOUBLE PRECISION :: tstart, tend, t_elapse
    
    
    GLOBAL_MAX = 0 
    maxHOTPO = 1e8

    num_threads = 4
    call OMP_SET_NUM_THREADS(num_threads)

    nstart = 2
    nend = 2000 

    tstart = OMP_GET_WTIME()




    !!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iterNMAX,j,highN)&
    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(PRIVATE)&
    !$OMP& SHARED(maxHOTPO,nstart,nend) REDUCTION(max:GLOBAL_MAX) 
    DO iterNMAX = nstart,nend
        n = iterNMAX
        highN = 0 

        DO j = 1,maxHOTPO 
        
            !! local max 
            IF(n > highN) THEN 
                highN = n 
            END IF 

            
            call hotpo(n)

            IF (n == 1) THEN 
                EXIT 
            END IF 


            !! global max 
            IF(highN > GLOBAL_MAX) THEN 
                GLOBAL_MAX = highN
            END IF  



        END DO 

    END DO
    !$OMP END PARALLEL DO 
    

    tend = OMP_GET_WTIME() 

    t_elapse = tend - tstart 


    WRITE(*,*) ""
    WRITE(*,*) "RUNNING WITH ", num_threads, " CORES"
    WRITE(*,*) "Time taken :",t_elapse  
    WRITE(*,*) "GLOBAL MAX = ", GLOBAL_MAX



END PROGRAM test 



SUBROUTINE hotpo(ndum)
    IMPLICIT NONE
    INTEGER*8, INTENT(INOUT) :: ndum

    IF (mod(ndum,2) == 0) THEN
        ndum = ndum / 2 
    ELSE
        ndum = (3 * ndum ) + 1
    END IF

END SUBROUTINE hotpo 
