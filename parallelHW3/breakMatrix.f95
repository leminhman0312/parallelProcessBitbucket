program decomp
    use mpi
    implicit none
    
    !! PROBLEM VARIABLES
    real, dimension(:,:),ALLOCATABLE :: main_matrix 
    integer :: i,j,ncol
    integer :: IHIGH, ILOW
    integer :: nblocks,ilocal,jglobal,k
    integer :: totalnum
    

    real,dimension(:,:),ALLOCATABLE :: a_local
    
    real,dimension(:),ALLOCATABLE :: x_global, x_local
    real,dimension(:),ALLOCATABLE :: a_local_1d

    real,DIMENSION(:),ALLOCATABLE ::y_local,y_global





    !! MPI VARIABLES
    integer :: rank, nprocs, ierr

  

   

    !! *****************************************  START OF MPI *********************************************************** !!
    !Start MPI
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

    

    ncol = 10

    allocate(main_matrix(ncol,ncol))
    allocate(x_global(ncol))
    allocate(y_global(ncol))


    ! Set up matrix
    do j = 1,ncol
        do i = 1,5
            main_matrix(i,j) = 55
        end do
    end do



    do j = 1,ncol
        do i = 6,ncol
            main_matrix(i,j) = 99
        end do
    end do


    if (rank == 0 ) then
        do i = 1,ncol
            print*, main_matrix(i,:)
        end do
    end if


    
   

    nblocks = FLOOR(REAL(ncol/nprocs))
    ALLOCATE(a_local(nblocks,ncol))
    ALLOCATE(x_local(nblocks))
    ALLOCATE(y_local(nblocks))


    totalnum = nblocks*ncol

    allocate(a_local_1d(totalnum))
   
   

     

    ! If nblocks is constant => broadcast nblocks to other processes
    call MPI_BCAST(nblocks,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    ! Then, we can also broadcast a_local => constant 2d array
    call MPI_BCAST(a_local,totalnum,MPI_REAL,0,MPI_COMM_WORLD,ierr) 
    call MPI_BCAST(a_local_1d,totalnum,MPI_REAL,0,MPI_COMM_WORLD,ierr) 



    ILOW  = (rank*nblocks)+1
    IHIGH = (rank+1)*nblocks

    do i = 1,nblocks
        a_local(i,:) = main_matrix(i+(ILOW-1),:)
    end do





    
    ! if (rank == 0 ) then
    !     !Calculate IHIGH, ILOW 
    !     ILOW  = (rank*nblocks)+1
    !     IHIGH = (rank+1)*nblocks

    !     ! print*, "I am process ", rank, " nblocks = ", nblocks, " ILOW = ", ILOW, " IHIGH = ", IHIGH


    !     ! Set up a_local chunk
    !     do ilocal = ILOW,IHIGH
    !         do i = 1,nblocks
    !             a_local(i,:) = main_matrix(ilocal,:) 
    !         end do
    !     end do

    !     !Set up X local
        
    !     do i = 1,nblocks
    !         x_local(i) = i
    !     end do

    !     ! print*, x_local


    !     ! Convert A_local into 1d
       
    !     ! a_local = transpose(a_local)
    !     ! a_local_1d = reshape(a_local, [totalnum] ) 



   
    ! else
    !     !Calculate IHIGH, ILOW other process
    !     ILOW  = (rank*nblocks)+1
    !     IHIGH = (rank+1)*nblocks     
        
    !     ! print*, "I am process ", rank, " nblocks = ", nblocks," ILOW = ", ILOW, " IHIGH = ", IHIGH

    !     ! Set up a_local chunk other process
    !     do ilocal = ILOW,IHIGH
    !         do i = 1,nblocks
    !             a_local(i,:) = main_matrix(ilocal,:) 
    !         end do
    !     end do


    !     ! Setup X_local other process


    !     do i = 1,nblocks
    !         x_local(i) = i
    !     end do
        



    !     ! ! ! Convert A_local into 1d other process
    !     ! a_local = transpose(a_local)
    !     ! a_local_1d = reshape(a_local, [totalnum] ) 




    ! end if



    

    ! All process needs whole X => use All gather so all proces has 1 big X_global = 12345 12345


    ! call MPI_ALLGATHER(x_local,nblocks,MPI_REAL,x_global,nblocks,MPI_REAL,MPI_COMM_WORLD,ierr)

    ! call MPI_ALLGATHER(x_local,1,MPI_REAL,x_global,1,MPI_REAL,MPI_COMM_WORLD,ierr)





    ! if (rank == 0 ) then
    !     print* ,"Rank = ", rank 
    !     print*, x_local
    ! else
    !     print* ,"Rank = ", rank
    !     print*, x_local
    ! end if




    

    do i = 1,nblocks
        y_local(i) = 0.0
        do j = 1,ncol
            ! y_local(i) = y_local(i) + a_local_1d((i-1)*ncol+j)*x_global(j)            !!  USE 1D A LOCAL
            y_local(i) = y_local(i) + a_local(i,j)*x_global(j)                          !!  USE 2D A LOCAL
        end do
    end do

    
    

    call MPI_ALLGATHER(y_local,nblocks,MPI_REAL,y_global,nblocks,MPI_REAL,MPI_COMM_WORLD,ierr)

    print*, ""

    if (rank == 0) then
        print*, y_global
    end if


























    




    

     !Deallocate array

    DEALLOCATE(a_local)
    DEALLOCATE(main_matrix)

    



    !Stop MPI
    call MPI_FINALIZE(ierr)

    !! *****************************************  END OF MPI *********************************************************** !!


   

end program decomp