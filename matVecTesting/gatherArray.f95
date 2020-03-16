program gathervector

    use mpi
    implicit none
    integer :: i, j, ierr,nprocs,rank, counter

    integer :: x_local

    integer,dimension(:),allocatable :: x_global


    !Start OPENMP
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

    allocate(x_global(nprocs))


    ! Each processor has a different data 
    if (rank == 0) then 
        x_local = 1
        print*, "I am proc ", rank, " and I have local data = ", x_local

        x_global(1) = x_local

    else
        x_local = 100
        print*, "I am proc ", rank, " and I have local data = ", x_local

        

    end if


    call MPI_GATHER(x_local,1,MPI_INT,x_global,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
   


   

    !PRINT OUT GLOBAL at P0

    print*, ""

    if (rank == 0 ) then
        do i = 1,nprocs
            print*, x_global(i)
        end do
    end if


    !END MPI

    call MPI_Finalize(ierr)






end program gathervector