program partB
    use omp_lib
    implicit none

    ! PROBLEM VARIABLES
    real,dimension(:,:),allocatable :: s2d
    real,dimension(:),ALLOCATABLE :: x_vector,y_vector, Sx_vector, Bx_vector
    integer :: ncol = 10
    real :: q= 0.15
    integer :: timecount
    integer :: nshape
    DOUBLE PRECISION :: tstart, tend, timetaken 

    ! CSR STUFF

    REAL,DIMENSION(:),ALLOCATABLE :: a_vec
    INTEGER,DIMENSION(:),ALLOCATABLE :: j_vec
    INTEGER,DIMENSION(:),ALLOCATABLE:: i_vec
    INTEGER,DIMENSION(:),ALLOCATABLE :: nsizezero
    INTEGER :: nsize, counter
    INTEGER :: i,j,k








    !----------------  SETTING UP S-----------------------------!
    interface 
        subroutine formS(nsize,smat)
            implicit none
            integer,INTENT(IN) :: nsize
            real,INTENT(OUT),DIMENSION(:,:),ALLOCATABLE :: smat

        end subroutine formS
    end interface  
    
   

    allocate(s2d(ncol,ncol))
    call formS(ncol,s2d)  






    !----------------  FORMING CSR ----------------------------!

    nsize = count (s2d/=0)      ! FIND # of NON ZEROS in S
    ALLOCATE(nsizezero(ncol))
    ALLOCATE(a_vec(nsize))      ! A has length # of non zero
    ALLOCATE(j_vec(nsize))      ! J has length # of non zero
    ALLOCATE(i_vec(ncol+1))     ! I has length NCOL+1

    counter = 1


    ! CALCULATE A, and jA, get non zero in each row
    do i = 1,ncol
        nsizezero(i) = count(s2d(i,:)/=0)
        do j = 1,ncol
            if (s2d(i,j) /= 0) then 
                a_vec(counter) = s2d(i,j)
                j_vec(counter) = j
                counter = counter + 1
            end if
        end do
    end do

    i_vec(1) = 1
    do i = 2,ncol+1
        i_vec(i) = i_vec(i-1) + nsizezero(i-1)
    end do



    !----------------  SETTING UP X----------------------------!


    allocate(x_vector(ncol))
    x_vector(:) = 1.00/ncol

    allocate(y_vector(ncol))  !set up y_vector as well
    y_vector(:) = 0.0   ! zero out   

    allocate(Bx_vector(ncol))  !Bx vector
    ! Base on property of X adds up to 1
    Bx_vector(:) = (1.00/ncol)







    ! !-----------------  MAIN MULTIPLY ------------------------------  !

    ALLOCATE(Sx_vector(ncol))       ! to store (1-q)Sx
    call cpu_time(tstart)
    do timecount = 1,1000
        ! Get Sx using CSR. Correct
        do i = 1,ncol
            Sx_vector(i) = 0.0
            !$OMP CRITICAL 
            do k = i_vec(i),i_vec(i+1)-1
                Sx_vector(i) = Sx_vector(i)+a_vec(k)*x_vector(j_vec(k))
            end do
            !$OMP END CRITICAL
            y_vector(i) = (1-q)*Sx_vector(i) + q/ncol
        end do

        !$OMP DO 
        ! Update X 
        do i = 1,ncol
            x_vector(i) = y_vector(i)
        end do
        !$OMP END DO
        
    end do

    call cpu_time(tend)


    timetaken = tend-tstart


    !$ timetaken = timeelapse/num_threads



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
        write(6,30) timetaken
        30 format(F16.12)

        write(6,*) "--------------------------------------------------------------------"
    end if

    if (ncol == 1600) then
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
        write(6,33) timetaken
        33 format(F16.12)

        write(6,*) "--------------------------------------------------------------------"
    end if

            







end program partB



subroutine formS(nsize,smat)
    implicit none
    integer ::  icol
    integer,INTENT(IN) :: nsize
    
    real,INTENT(OUT),DIMENSION(:,:),ALLOCATABLE :: smat


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






