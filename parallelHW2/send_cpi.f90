program midpoint_speedup
    
    include 'mpif.f'

    integer :: ierr, rank, nprocs,n,i, isend, nsend, nreceive
    integer,dimension(MPI_STATUS_SIZE) :: status1
    double precision :: real_PI
    double precision :: mypi, pi, h, sum, x
    double precision startwtime, endwtime
    double precision :: mypi0, mypi_recv


    !nprocs = 1

    startwtime = 0.0
    real_PI = 3.141592653589793238462643

  


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


    ! SEND CALL
    if (rank == 0) then
        do i = 1,nprocs-1
            call MPI_SEND(n,1,MPI_INT,i,1,MPI_COMM_WORLD,ierr)
        end do
    ! RECEIVE CALL
    else
        call MPI_RECV(n,1,MPI_INT,0,1,MPI_COMM_WORLD,status1,ierr)
    end if


    ! COMPUTE PARTIAL SUM
    h = 1.0/(1.0*n)
    sum = 0.0
    do i = rank+1,n,nprocs
        x = h*((1.0)*i-0.5)
        sum = sum + 4.0/(1.0+x*x)
    end do
    mypi = h*sum


    if (rank == 0) then 
        pi = mypi
    end if


    
        
    if (rank == 0) then

        write(6,99) rank, mypi
        99 format('Initial, process ', I2, ' has', f16.12)

        do i = 1,nprocs-1
            ! receive my pi
            call MPI_RECV(mypi_recv,1,MPI_DOUBLE_PRECISION,i,2,MPI_COMM_WORLD,status1,ierr)

            write(6,100) rank, mypi_recv,i
            100 format('Process', I2, ' received', f16.12, ' from process ', I2)

            pi = pi + mypi_recv
            
            
        end do

    else
           
        

        ! send mypi
        call MPI_SEND(mypi,1,MPI_DOUBLE_PRECISION,0,2,MPI_COMM_WORLD,ierr)
        write(6,101) rank, mypi
        101 format('Process ',I2,' sent ', f16.12, ' to 0')
    end if






    


    



    if (rank == 0) then
        endwtime =  MPI_WTIME()
        write (6,200) pi, abs(pi-real_PI), nprocs
    200 format(' ','pi is approximately ',f16.12,',Error is ',f16.12,' nprocs = ', I2)
        write (6,300) endwtime-startwtime
    300 format(' ','wall click time = ', f16.12)

    end if

   
   

    ! !Close OPENMP
    call MPI_FINALIZE(ierr)



 


end program midpoint_speedup




