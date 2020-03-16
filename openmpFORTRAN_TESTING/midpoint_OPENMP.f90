program midpoint_speedup


    use mpi
    implicit none

    integer :: ierr, rank, nprocs,n,i, isend, nsend, nreceive
    integer,dimension(MPI_STATUS_SIZE) :: status1
    double precision :: real_PI
    double precision mypi, pi, h, sum, x
    double precision startwtime, endwtime

    

    !nprocs = 1

    startwtime = 0.0
    real_PI = 3.141592653589793238462643

    ! OPEN UP THE FILE
    open (21, file = 'resultMPI.dat',action='write',position='append')



    !Start OPENMP
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)


    n = 0
    ! SET PROBLEM SIZE AT PROCESS 0
    if (rank == 0) then
        n =  600000000
        startwtime =  MPI_WTIME()
    end if


   


    !BROADCAST PROBLEM SIZE TO ALL
    call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)


 



    ! COMPUTE PARTIAL SUM
    h = 1.0/(1.0*n)
    sum = 0.0
    do i = rank+1,n,nprocs
        x = h*((1.0)*i-0.5)
        sum = sum + 4.0/(1.0+x*x)
    end do
    mypi = h*sum

    !COMBINE THE LOCAL SUM, ADD THEM AT PROCESS 0
    call MPI_REDUCE(mypi, pi, 1, MPI_DOUBLE_PRECISION,&
        MPI_SUM, 0, MPI_COMM_WORLD, ierr)


    !AT PROCESS 0, PRINT OUT RESULTS

    write(21,*) " "
    write(21,*) "----------------------------------------------------------------------"

    if (rank == 0) then
        endwtime =  MPI_WTIME()
        write (6,200) pi, abs(pi-real_PI), nprocs
        write(21,200) pi, abs(pi-real_PI),nprocs
    200 format(' ','pi is approximately ',f16.12,',Error is ',f16.12,' nprocs = ', I2)
        write (6,300) endwtime-startwtime
        write (21,300) endwtime-startwtime
    300 format(' ','wall click time = ', f16.12)

    end if

    !Close OPENMP
    call MPI_FINALIZE(ierr)

    ! Close the file
    close(21)



end program midpoint_speedup




