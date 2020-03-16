
PROGRAM Compute_PI
    IMPLICIT NONE

    interface
        FUNCTION f(a)
        double precision a
        double precision f
        END FUNCTION
    end interface


    INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS     

    INTEGER           N, i
    INTEGER           id, num_threads
    DOUBLE PRECISION  w, x, sum
    DOUBLE PRECISION  pi, mypi


    N = 50000000         !! Number of intervals
    w = 1.0d0/N          !! width of each interval

    sum = 0.0d0

    !$OMP    PARALLEL PRIVATE(i, id, num_threads, x, mypi)

    num_threads = omp_get_num_threads()
    id = omp_get_thread_num()

    mypi = 0.0d0;

    DO i = id,   N-1,   num_threads  !! Parallel DO LOOP
        x = w * (i + 0.5d0)
        mypi = mypi + w*f(x)
    END DO


    !$OMP CRITICAL   !!Executed by ONE thread at a time 
    pi = pi + mypi
    !$OMP END CRITICAL


    !$OMP    END PARALLEL

    PRINT *, "Pi = ", pi

END PROGRAM Compute_PI



FUNCTION f(a)
    IMPLICIT NONE

    double precision a
    double precision f

    f = 2.d0 / SQRT(1.d0 - a*a)
END
