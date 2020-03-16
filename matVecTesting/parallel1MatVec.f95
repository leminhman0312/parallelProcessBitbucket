program matvec

    !! SOLVING Ax = y
    
    use mpi
    implicit none
    
    !Declare MPI Variables
    integer :: ierr,rank,nprocs
    double precision startwtime, endwtime



    !Declare other data
    integer :: i,j
    real,dimension(3) :: x_global
    real :: x_local


    integer,dimension(3) :: n_global
    integer :: n_local



    real, dimension(9) :: a_local  ! spread a out into 1 single vector
    real,dimension(3) :: y_local



    a_local(1) = 1
    a_local(2) = 2
    a_local(3) = 3

    a_local(4) = 4
    a_local(5) = 5
    a_local(6) = 6

    a_local(7) = 7
    a_local(8) = 8
    a_local(9) = 9





    !Start MPI
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)


    if (rank == 0) then
        startwtime =  MPI_WTIME() 
        x_local = 10
        n_local = 3
        print*,"Rank ", rank, "has local data = ", x_local
    else
        x_local = 10+rank
        n_local = 3
        print*,"Rank ", rank, "has local data = ", x_local
    end if



   
    ! OR, USE MPI_ALLGATHER

    call MPI_ALLGATHER(x_local,1,MPI_REAL,x_global,1,MPI_REAL,MPI_COMM_WORLD,ierr)


    ! call MPI_ALLGATHER(n_local,1,MPI_INT,n_global,1,MPI_INT,MPI_COMM_WORLD,ierr)




    do i = 1,3
        y_local(i) = 0.0
        do j = 1,3
            y_local(i) = y_local(i) + a_local((i-1)*n_local+j)*x_global(j)
        end do
    end do 


    print*, ""

    if (rank == 0) then
        endwtime =  MPI_WTIME()
        do i = 1,3
            write(6,300) y_local(i)
            300 format(F6.2)
        end do
        write(6,400) (endwtime-startwtime)
        400 format('', 'Time taken', F16.12)
    end if

    !End MPI
    call MPI_Finalize(ierr)




end program matvec