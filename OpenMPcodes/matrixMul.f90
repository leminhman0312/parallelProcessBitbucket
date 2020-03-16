program  matrix
    use omp_lib
    implicit none
    
    INTEGER :: num_threads = 8
    INTEGER, PARAMETER :: m = 100,n = 100, o = 100
    INTEGER :: i,j,k
    real*8,DIMENSION(1:m,1:n) :: a
    real*8,DIMENSION(1:n,1:o) :: b
    real*8,DIMENSION(1:m,1:n) :: c
    real*8 :: t1, t2, ep

    a = 1.0
    b = 1.0 
    c = 0.0



    !$ call omp_set_num_threads(num_threads)

    !start
    call cpu_time(t1)


    !$omp parallel do 
    do i = 1,m 
        do j = 1,o
            do k = 1,n
                c(i,j) = c(i,j)+a(i,k)*b(k,j)
            end do
        end do
    end do
    !$omp end parallel do


    !stop
    call cpu_time(t2)
    ep = t2-t1

    !$ ep = ep/num_threads
    !$ print*, "Parallel Mode On" 


    print*, " Time taken in seconds: ", ep
    print*, " A Sample value : ", c(10,10)


end program  matrix