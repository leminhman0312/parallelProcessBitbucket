program midpoint_speedup


    ! include 'mpif.h'
    
    use mpi
    implicit none
    integer :: ierr, rank, nprocs,n,i
    integer,dimension(MPI_STATUS_SIZE) :: status1
    double precision startwtime, endwtime
    


    !nprocs = 1

    startwtime = 0.0


    !Start OPENMP
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)


   
    ! BROADCAST PROBLEM SIZE TO ALL
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


    if (rank == 0) then
        endwtime =  MPI_WTIME()
        write (6,300) endwtime-startwtime
        300 format(' ','wall click time = ', f16.12)

    end if

    
    

    ! Close OPENMP
    call MPI_FINALIZE(ierr)



end program midpoint_speedup




