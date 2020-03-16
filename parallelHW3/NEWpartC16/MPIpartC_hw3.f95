program partC
    use mpi
    ! include 'mpif.h'
    implicit none

    !! DECLARE PROBLEM VARIABLES

    integer :: i,j,iTrue,k
    DOUBLE PRECISION,dimension(:,:),ALLOCATABLE:: s_local, s_global, b_local, b_global, g_local, g_global
    DOUBLE PRECISION :: q = 0.15
    DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: x_local, y_local, y_global, x_global
    integer :: ncol = 1600

    DOUBLE PRECISION :: ytempt

    integer :: timecount 

    DOUBLE PRECISION :: ymin, yminlocal

    DOUBLE PRECISION :: ymax, ymaxlocal

    DOUBLE PRECISION :: tstart,tend,timeelapse





    !! DECLARE MPI VARIABLES
    integer :: rank,nprocs,ierr



    !! DECLARE OTHER VARIABLES
    integer :: IHIGH,ILOW
    integer :: nblocks
    integer :: jDia, jSuper, jSub



    !Start MPI
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)





    nblocks = FLOOR(REAL(ncol/nprocs))

    ALLOCATE(s_local(nblocks,ncol))
    ALLOCATE(b_local(nblocks,ncol))
    ALLOCATE(g_local(nblocks,ncol))
    ALLOCATE(x_local(nblocks))
    ALLOCATE(y_local(nblocks))
    ALLOCATE(y_global(ncol))
    ALLOCATE(x_global(ncol))



    ! If nblocks is constant => broadcast nblocks to other processes
    ! call MPI_BCAST(nblocks,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)



    ILOW  = (rank*nblocks)+1
    IHIGH = (rank+1)*nblocks

    do i = 1, nblocks

        !! --------- SETTING UP S LOCAL -------------------------------  !!
        iTrue =i+(ILOW-1) 
        jDia = iTrue
        jSuper = jDia + 1
        jSub = jDia - 1
        
        s_local(i,:) = 0.0
        s_local(i,jDia) = 0.0
    
        !! Special cases
        if (jDia /= 1) then 
            s_local(i,jSub) = 0.5
        end if 

        if (jDia /= ncol) then 
            s_local(i,jSuper) = 0.5
        end if 


        if (iTrue == ncol) then
            !! More special points
            s_local(i,1) = 0.5
        end if


        if (iTrue == ncol-1) then
            !! More special points 
            s_local(i,ncol) = 1.0
        end if


        !!----------- SETTING UP Y LOCAL ----------------------------- !!
        y_local(:) = 0.0
    end do

    do i=1,nblocks

        !!------------ SETTING UP X LOCAL------------------------------ !!

        x_local(i) = 1.00/ncol
    end do





    tstart = MPI_WTIME()    
    do timecount = 1,1000

        call MPI_Allgather(x_local, nblocks,MPI_DOUBLE_PRECISION,x_global,nblocks,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)


        do i = 1,nblocks
            y_local(i) = 0.0

            do j = 1,ncol
                y_local(i) = y_local(i) + (1-q)*s_local(i,j)*x_global(j)
                
            end do
            y_local(i) = y_local(i) + (q/DBLE(ncol)) 

        end do

        !   UPDATE X 
        do k = 1,nblocks
            x_local(k) = y_local(k)
        end do
    end do
    tend = MPI_WTIME()

    ! tend = MPI_WTIME()
        



    !!------------ CALCULATE MAX  and MIN
    ymaxlocal = maxval(y_local) 
    call MPI_Reduce(ymaxlocal,ymax,1,MPI_DOUBLE_PRECISION,MPI_MAX, 0, MPI_COMM_WORLD,ierr)


    
    yminlocal = minval(y_local) 
    call MPI_Reduce(yminlocal,ymin,1,MPI_DOUBLE_PRECISION,MPI_MIN, 0, MPI_COMM_WORLD,ierr)





    !! GATHER Y LOCAL TO Y GLOBAL

    ! call MPI_Allgather(y_local, nblocks,MPI_DOUBLE_PRECISION,y_global,nblocks,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

    ! if (rank == 0) then 
    !     ! print*, y_global
    !     print*, "Ymax = ", ymax
    !     print*, "Ymin = ", ymin
    ! end if



    if (ncol == 16) then 

        if (rank == 0 ) then 


            WRITE(6,*) " Y = "
            do i = 1,ncol
                write(6,20) y_global(i)
                20 format('  ',F7.5)
            end do

            write(6,*) "YMAX = "
            WRITE(6,25) ymax
            25 format('  ',F7.5)


            write(6,*) " " 

            write(6,*) "YMIN =  "
            WRITE(6,26) ymin
            26 format('  ',F7.5)

            write(6,*) " Time Elapsed ="
            write(6,30) (tend-tstart)
            30 format(F16.12)



                

        end if
	 
        else 

        if (rank == 0 ) then 

            write(6,*) "YMAX = "
            WRITE(6,31) ymax
            31 format('  ',F7.5)


            write(6,*) " " 

            write(6,*) "YMIN =  "
            WRITE(6,32) ymin
            32 format('  ',F7.5)
            
            
            write(6,*) " Time Elapsed ="
            write(6,33) (tend-tstart)
            33 format(F16.12)


        end if
    end if
    



    


    DEALLOCATE(s_local)
    DEALLOCATE(b_local)
    DEALLOCATE(g_local)


    !Stop MPI
    call MPI_FINALIZE(ierr)






end program partC
