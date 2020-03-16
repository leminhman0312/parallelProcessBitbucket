program main

    !! Program to set N a private var
    !! N will not be used in the parallel loop
    
    use omp_lib
    implicit none 
    integer :: N !This is a shared variable
    integer :: num_threads = 4

    !$ call omp_set_num_threads(num_threads)
        
    N = 1001
    print*, "Before parallel section: N = ", N 

    !$OMP PARALLEL PRIVATE(N)
    N = N + 1
    print*, "Inside parallel section: N = ", N 
    !$OMP END PARALLEL


    print*, "After parallel section: N = ", N

end program main 
