program serialmatvec2d


    implicit none
    integer :: i, j 
    real, dimension(3,3) :: a
    real, dimension(3) :: x,y

    a(1,1) = 1
    a(1,2) = 2
    a(1,3) = 3

    a(2,1) = 4
    a(2,2) = 5
    a(2,3) = 6

    a(3,1) = 7
    a(3,2) = 8
    a(3,3) = 9


    print*, "THIS IS THE MATRIX A"
    do i = 1,3
            print*, a(i,:)
    end do


    x(1) = 10
    x(2) = 11
    x(3) = 12


    ! matrix multiplication 


    do i = 1,3
        y(i) = 0
        do j = 1,3
            y(i) = y(i) + a(i,j)*x(j)
        end do
    end do



    print*, ""
    do i = 1,3
        print*, y(i)
    end do






end program serialmatvec2d