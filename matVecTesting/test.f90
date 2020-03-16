
program test

    use mpi
    implicit none
    integer, dimension(8) :: global      ! only root has this
    integer, dimension(2) :: local       ! everyone has this
    integer, parameter    :: root = 0
    integer :: rank, comsize
    integer :: i, ierr

    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, comsize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    if (rank == root) then
        global = [ (i, i=1,8) ]
    endif

    call MPI_Scatter(global, 2, MPI_INTEGER, &    ! send everyone 2 ints from global
                    local,  2, MPI_INTEGER, &    ! each proc recieves 2 into
                    root,                   &    ! sending process is root,
                    MPI_COMM_WORLD, ierr)        ! all procs in COMM_WORLD participate



    if (rank == root) then 
        print*, "THIS IS ROOT"
        do i = 1,8
            print*, global(i)
        end do
    else
        print*, "THIS IS NOT ROOT, RANK = ", rank
        do i = 1,8
            print*, global(i)
        end do
    end if


    call MPI_FINALIZE(ierr)


end program test