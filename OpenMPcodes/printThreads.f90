program main
    implicit none
    integer :: nthreads, myid, maxthreads
    integer,external :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS, OMP_GET_MAX_THREADS


    !$OMP PARALLEL private(nthreads,myid)  !! Make these private

    myid = OMP_GET_THREAD_NUM()

    print*, "Hello I am thread ", myid

    ! Print stuff out at processor 0. 
    if (myid == 0) then
        nthreads = OMP_GET_NUM_THREADS()
        maxthreads = OMP_GET_MAX_THREADS()
        print*, "Number of threads = ", nthreads
        print*, "Maxthreads = ", nthreads

    end if

    !$OMP END PARALLEL 







end program main