      program cpi

      include 'mpif.h'


      integer n, myid, numprocs, i, ierr
      double precision PI25DT
      double precision mypi, pi, h, sum, x
      double precision startwtime, endwtime

      startwtime = 0.0
      PI25DT = 3.141592653589793238462643

C     Initialize MPI
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)

      n = 0
C     Process zero sets problem size and starts timer
      if (myid.EQ.0) then
	   n=600000000
	   startwtime =  MPI_WTIME()
      end if

C     Broadcast problem size to all processors
      call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

C     Compute local sum
      h   = 1.0 / (1.0*n)
      sum = 0.0
      do 100 i=myid+1,n,numprocs
            x = h * ((1.0)*i - 0.5)
            sum = sum + 4.0/(1.0+x*x)
 100  continue
      mypi = h * sum

C     Combine local sums
      call MPI_REDUCE(mypi, pi, 1, MPI_DOUBLE_PRECISION,
     &                      MPI_SUM, 0, MPI_COMM_WORLD, ierr)

C     Process zero prints results
            if (myid.EQ.0) then
		endwtime =  MPI_WTIME()
                write (6,200) pi, abs(pi-PI25DT)
 200            format(' ','pi is approximately ',f16.12,', 
     &                      Error is ',f16.12)
                write (6,300) endwtime-startwtime
 300            format(' ','wall click time = ', f16.12)
            end if
    
      call MPI_FINALIZE(ierr)
      end
