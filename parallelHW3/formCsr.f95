PROGRAM test_pack_2
    implicit none
    real,dimension(:),ALLOCATABLE :: amat
    INTEGER,DIMENSION(:),ALLOCATABLE :: jamat
    real :: mymat(5,5)
    INTEGER,DIMENSION(5) :: nsizezero
    INTEGER,DIMENSION(6) :: iamat
    INTEGER :: i,j
    integer :: nsize
    INTEGER :: counter
    
    

    mymat(:,:) = 0
    mymat(1,2) = 1
    
    mymat(2,1) = 0.5
    mymat(2,3) = 0.5
    mymat(2,4) = 0.3

    mymat(3,1) = 0.5
    mymat(3,4) = 0.3

    mymat(4,3) = 0.5
    mymat(4,5) = 1
    mymat(5,4) = 0.3
    
    nsize = count (mymat/=0)

    allocate(amat(nsize))
    allocate(jamat(nsize))

    counter = 1





    do i = 1,5
        nsizezero(i) = count(mymat(i,:)/=0)

        do j = 1,5
            if (mymat(i,j) /= 0) then 
                amat(counter) = mymat(i,j)
                jamat(counter) = j
                counter = counter + 1
            end if

        end do
    end do



    iamat(1) = 1
    do i = 2,6
        iamat(i) = iamat(i-1) + nsizezero(i-1)
    end do



    

    print*, "FORM CSR", size(mymat)

    








END PROGRAM