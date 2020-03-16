program partA
    use omp_lib
    implicit none

    ! PROBLEM VARIABLES
    real*8,dimension(:,:),allocatable :: s2d,bmat,gmat
    real*8,dimension(:),ALLOCATABLE :: gmat1d,x_vector,y_vector
    integer :: ncol = 1600
    integer :: nshape 
    integer :: i,j
    real*8 :: q
    integer :: timecount
    DOUBLE PRECISION :: tstart,tend,timeelapse
    CHARACTER :: np
    REAL*8 :: sumi




    !OPEN MP VARIABLES

    integer :: num_threads


    !----------------  SETTING UP S-----------------------------!
    interface 
        subroutine formS(nsize,smat)
            implicit none
            integer,INTENT(IN) :: nsize
            real*8,INTENT(OUT),DIMENSION(:,:),ALLOCATABLE :: smat

        end subroutine formS
    end interface

    num_threads = 4
    call omp_set_num_threads(num_threads)


    !CALL get_environment_variable("OMP_NUM_THREADS", np)
    !num_threads = ICHAR(np)
    !read(np, *) num_threads
    !call omp_set_num_threads(num_threads)
    
    !num_threads = omp_get_num_threads()    

    !print*, "Running with ", num_threads, " threads ", np


    allocate(s2d(ncol,ncol))
    call formS(ncol,s2d)    

    !----------------  END SETTING UP S------------------------!




    !----------------  SETTING UP G----------------------------!

    q = 0.15

    allocate(bmat(ncol,ncol))

    ! First, setup Bmatrix, each element = 1/numPage
    bmat(:,:) = 1.00/ncol



    ! Now, forming G
    allocate(gmat(ncol,ncol))
    ! gmat = (1.00-q)*s2d + q*bmat 



    do i = 1,ncol
        do j = 1,ncol
            gmat(i,j) = (1.00-q)*s2d(i,j) + q*bmat(i,j)
        end do
    end do

    

    ! Convert G into 1d
    nshape = ncol**2
    allocate(gmat1d(nshape))
    gmat = transpose(gmat)
    gmat1d = reshape(gmat, [nshape] ) 

    !----------------  END SETTING UP G------------------------!



    !----------------  SETTING UP X----------------------------!


    allocate(x_vector(ncol))
    x_vector(:) = 1.00/ncol

    allocate(y_vector(ncol))  !set up y_vector as well
    y_vector(:) = 0.0   ! zero out   

    !----------------  END SETTING UP X------------------------!




    !----------------- MAIN MATVEC ---------------------------------------------------!

    tstart =  omp_get_wtime()
    
    do timecount = 1,1000
	!$OMP PARALLEL DO private(sumi)
        do i = 1,ncol
 	    sumi = 0.0
            do j = 1,ncol
                sumi = sumi + gmat1d((i-1)*ncol+j)*x_vector(j)
            end do
	    y_vector(i) = sumi
        end do  
	!$OMP END PARALLEL DO
        
        !$OMP DO
        do i = 1,ncol
            x_vector(i) = y_vector(i)  !UPDATE X
        end do
        !$OMP END DO
    end do



   
    tend =  omp_get_wtime()
    
    timeelapse = tend-tstart

    !!!$ timeelapse = timeelapse/num_threads
    
    !----------------- END MAIN MATVEC -------------------------------------------------!


    if (ncol == 16) then 

        write(6,*) "--------------------------------------------------------------------"


        WRITE(6,*) " Y = "
        do i = 1,ncol
            write(6,20) y_vector(i)
            20 format('  ',F7.5)
        end do

        write(6,*) " "

        write(6,*) "YMAX = "
        WRITE(6,25) MAXVAL(y_vector)
        25 format('  ',F7.5)


        write(6,*) " "

        write(6,*) "YMIN =  "
        WRITE(6,26) MINVAL(y_vector)
        26 format('  ',F7.5)



        write(6,*) " "

        write(6,*) " Time Elapsed ="
        write(6,30) timeelapse
        30 format(F16.12)

        write(6,*) "--------------------------------------------------------------------"
    else

!    if (ncol == 1600) then
        write(6,*) "--------------------------------------------------------------------"

        write(6,*) "YMAX = "
        WRITE(6,31) MAXVAL(y_vector)
        31 format('  ',F7.5)


        write(6,*) " "

        write(6,*) "YMIN =  "
        WRITE(6,32) MINVAL(y_vector)
        32 format('  ',F7.5)

        write(6,*) " "

        write(6,*) " Time Elapsed ="
        write(6,33) timeelapse
        33 format(F16.12)

        write(6,*) "--------------------------------------------------------------------"
    end if



    














    !----------------- END MAIN MATVEC -------------------------!















    !! DEALLOCATE 
    deallocate(gmat)
    deallocate(bmat)
    deallocate(s2d)
    
end program partA 



subroutine formS(nsize,smat)
    implicit none
    integer ::  icol
    integer,INTENT(IN) :: nsize
    
    real*8,INTENT(OUT),DIMENSION(:,:),ALLOCATABLE :: smat


    ALLOCATE(smat(nsize,nsize))

    !Super diagonal
    do icol = 1,nsize-1 !loop row
       
        smat(icol,icol+1) = 0.5
      
    end do

    !Sub diagonal

    do icol = 2,nsize
        smat(icol,icol-1) = 0.5
    end do
   
    smat(nsize,1) = 0.5
    smat(nsize-1,nsize) = 1

end subroutine






