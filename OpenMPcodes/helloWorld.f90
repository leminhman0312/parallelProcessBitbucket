program hello
    use omp_lib 
    implicit none
    integer :: num_threads = 4
    integer :: thread_num =   0 
    !$ call omp_set_num_threads(num_threads)
    print*, "Hello"
    print*, "Number of threads used = ", num_threads

    !$omp parallel 
        !$omp critical  
            !$ thread_num = omp_get_thread_num()
                print*, "Hello from ", thread_num 
        !$omp end critical 
    !$omp end parallel 
    

end program hello 
