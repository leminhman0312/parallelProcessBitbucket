program tri1d
    implicit none

    integer :: i,j, N=5, indDiag, indSuper,indSub
    real,DIMENSION(:),ALLOCATABLE :: s_arr, b_arr


    ALLOCATE(s_arr(N**2))
    ALLOCATE(b_arr(N**2))



    do i = 1,N**2
        b_arr(i) = 100
    end do


    do i = 2,N
        indDiag = i+(i-1)*N
        indSuper = indDiag + 1 
        indSub = indDiag - 1
        s_arr(indDiag) = 0.0        ! main diagonal
        s_arr(indSuper) = 0.5       ! super diagonal
        s_arr(indSub)  = 0.5        ! sub diagonal
    end do


    !First element of Supper Diagonal = 0.5

    s_arr(2) = 0.5

    ! Last element of Supper diagonal = 1 

    s_arr(N**2-N) = 1

    ! Last column, first element is 0.5

    s_arr(N**2-N+1) = 0.5





    ! print*, ""

    do i = 1,N**2
        print*, i, s_arr(i)
    end do



    DEALLOCATE(s_arr)


end program tri1d