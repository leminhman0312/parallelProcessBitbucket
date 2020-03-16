program pingpong


    ! include 'mpif.h'
    
    use mpi
    implicit none
    integer :: ierr, rank, nprocs,i,j
    integer,dimension(MPI_STATUS_SIZE) :: status1
    double precision startwtime, endwtime
    real,dimension(10001):: data
    integer,dimension(21):: narr
    integer :: left, right, current
    real :: timediff

    real,dimension(21):: timearr


    

    startwtime = 0.0


    !Start OPENMP
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)




    ! Populate sizearr
    narr(1) = 1
    do i = 2,21
        narr(i)=narr(i-1)+500
    end do




    do i = 1,21
        startwtime = MPI_WTIME()
        do j = 1,100
            if (rank == 0) then 
                current = rank
                left = current-1
                right = current+1
                ! Send from 0 to 1 PING
                call MPI_SEND(data,narr(i),MPI_REAL,right,1,MPI_COMM_WORLD,ierr)
                ! print*, "Rank ", current, "send ", narr(i), "# of data to ", right
            
                ! 0 receive from 1 PONG
                call MPI_RECV(data,narr(i),MPI_REAL,right,2,MPI_COMM_WORLD,status1,ierr)
                ! print*, "Rank ", current, "receive ", narr(i), "# of data ", rank


            else if (rank == 1) then 

                current = rank
                left = current-1
                right = current+1

                !1 receive from 0 PING
                call MPI_RECV(data,narr(i),MPI_REAL,left,1,MPI_COMM_WORLD,status1,ierr)
                ! print*, "Rank ", rank, "receive ", narr(i), "# of data from ", left


                ! 1 send back to 0 PONG
                call MPI_SEND(data,narr(i),MPI_REAL,left,2,MPI_COMM_WORLD,ierr)
                ! print*, "Rank ", rank, "send ", narr(i), "# of data to ", left

            end if
        end do

       
        endwtime = MPI_WTIME()
        timediff = endwtime-startwtime
        timearr(i) = timediff

        


    end do





    open (21, file = 'result.dat',action='write')
   


    !AT PROCESS 0, PRINT OUT RESULTS

    if (rank ==  0) then 
        do i = 1,21

            write(21,300) narr(i), timearr(i) 
            300 format(' ', I5 ,F16.12)
        end do
    end if
  


    close(21)
    

    ! Close OPENMP
    call MPI_FINALIZE(ierr)


    !from new mint



end program pingpong




