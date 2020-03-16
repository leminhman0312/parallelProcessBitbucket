program midpoint
    implicit none
    ! ******************** INTERFACE FUNCTION **************!
    interface
        function f1 ( x )
            real*8 :: f1
            real , intent ( in ) :: x
        end function f1
    end interface


    interface
        function f2 ( x )
            real*8 :: f2
            real , intent ( in ) :: x
        end function f2
    end interface
    ! *************************DECLARE VARIABLES **********!
    ! Declaring basic variables
    real :: lower , upper , delx , n
    integer :: i , j
    real*8 :: sum
    integer :: status
    integer :: numberofN
    real :: stopT , startT
    real :: timediff
    real*8 :: analytical
    ! Declaring arrays for storage
    real*8 , dimension (:) , allocatable :: narr
    real*8 , dimension (:) , allocatable :: resultarr
    real*8 , dimension (:) , allocatable :: analyticarr
    real*8 , dimension (:) , allocatable :: errarr
    real*8 , dimension (:) , allocatable :: timearr
    real*8 , dimension (:) , allocatable :: xmid



    ! *************************PUT NUMBERS IN **************!
    ! DEFINE THE UPPER AND LOWER BOUND OF THE INTEGRAL
    lower = 0.0
    upper = 1.0


    ! DEFINE THE ANALYTICAL RESULT , WE KNOW IT IS PI
    analytical = 4.0* atan (1.0)


    ! ALLOCATE STORAGE ARRAYS
    numberofN = 6.00
    allocate ( narr ( numberofN ))
    allocate ( resultarr ( numberofN ))
    allocate ( analyticarr ( numberofN ))
    allocate ( errarr ( numberofN ))
    allocate ( timearr ( numberofN ))


    ! POPULATE ANALYTICAL
    do j = 1 , numberofN
        analyticarr ( j ) = analytical
    end do

    ! *************************MAIN SOLVER ******************!

    n = 2
    sum = 0.0
    do while(n < 200000)
        do j = 1,numberofN
            call cpu_time(startT)
            ! Calculate dx
            delx = (upper-lower)/n
            allocate(xmid(int(n)))
            do i = 1,int(n)
                !populate xmid
                xmid(i) = lower + (0.5*delx)+(i-1)*delx
                ! midpoint sum
                sum = sum + delx*(f2(xmid(i)))
            end do
            call cpu_time(stopT)
            timediff = stopT-startT

            !Store into arrays
            timearr(j) = timediff
            narr(j) = n
            resultarr(j) = sum
            !calulate error
            errarr ( j ) = abs ( resultarr ( j ) - analyticarr ( j ))
            n = n *10
            deallocate ( xmid , stat = status )
            sum = 0.0
        end do
    end do



    ! ******************** WRITE TABLE TO FILE ************!

    ! OPEN UP FILE
    open (21, file = 'result.dat',action='write')

    !Write to file
    write (21 ,*) ""
    write (21 , " (6x , a1 ,10 x , a9 ,16x , a10 ,10x , a5 ,15x , a7 ) " ) &
        "N " , "NUMERICAL " , "ANALYTICAL " , " ERROR " , "TIME ( s ) "

    write(21,*) ''

    do i = 1,6
        write (21 , " ( T2 , I6 , T15 , F14 .12 , T40 , F14 .12 , T60 , F9 .6 , T80 , F14 .12) ") &
        int (narr (i)) , resultarr (i), analyticarr(i), errarr(i), timearr(i)
    end do

    write(21,*)''


    !DEALLOCATING THE ARRAYS
    deallocate ( narr )
    deallocate ( resultarr )
    deallocate ( errarr )
    deallocate ( analyticarr)

    ! CLOSE THE FILE
    close(21)




end program midpoint


    ! *************************FUNCTION PROTOTYPES *************!

    function f1(x)
        implicit none
        real*8 :: f1
        real,intent(in) :: x
        f1 = 4.00/(1+x**2)

    end function f1


    function f2(x)
        implicit none
        real*8 :: f2
        real,intent(in) :: x
        f2 = 4.00*sqrt(1-x**2)

    end function f2
