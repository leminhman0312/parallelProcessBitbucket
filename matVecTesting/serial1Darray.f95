program serial1Darray

    implicit none
    integer :: i, j 
    real, dimension(9) :: a  ! spread a out into 1 single vector
    real, dimension(3) :: x,y

    a(1) = 1
    a(2) = 2
    a(3) = 3

    a(4) = 4
    a(5) = 5
    a(6) = 6

    a(7) = 7
    a(8) = 8
    a(9) = 9


    x(1) = 10
    x(2) = 11
    x(3) = 12



    ! matrix multiplication 


    do i = 1,3
      
        y(i) = 0
        do j = 1,3
    
            y(i) = y(i) + a((i-1)*3+j)*x(j)
            

            ! print*, "Y(i) = ", y(i), "at i = ", i , "j = ", j 
        end do
    end do




    print*, ""
    do i = 1,3
        print*, y(i)
    end do









end program serial1Darray