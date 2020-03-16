program test

    use mpi 
    implicit none



    !! DECLARE MPI VARIABLES
    integer :: rank,nprocs,ierr

    integer :: nblocks = 2, ncol = 10

    integer, dimension(:),ALLOCATABLE :: x_local, x_global 

    !Start MPI
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

    allocate(x_local(nblocks))
    allocate(x_global(ncol))

    x_local(:) = 2

    call MPI_ALLgather(x_local,nblocks,MPI_INTEGER,x_global,ncol,MPI_INTEGER,MPI_COMM_WORLD,ierr)







    ! End MPI

    call MPI_FINALIZE(ierr)


end program test
