program send_recv

    use mpi
    implicit none


    integer :: ierr, rank,nprocs,n0,n1

    integer, dimension(MPI_STATUS_SIZE) :: status1

    call MPI_INIT(ierr)

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)

    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)



    ! rank = 0
    
    if (rank == 0) then
        n0 = 5
        print*, "Rank = ", Rank, " data = ", n0
        call MPI_SEND(n0,1,MPI_INT,1,1,MPI_COMM_WORLD,ierr)   
        print*, "Rank ", Rank, "sending to 1"
    end if


    ! rank = 1

    if (rank == 1) then 
        call MPI_RECV(n1,1,MPI_INT,0,1,MPI_COMM_WORLD,status1,ierr)
        print*, "Rank", Rank, "has data from 0, data1 = ", n1

    end if


    call MPI_FINALIZE(ierr)

end program send_recv