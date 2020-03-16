program main
	use omp_lib
	implicit none 
	integer :: N		!This is a shared variable
	integer :: num_threads = 4

	!$ call omp_set_num_threads(num_threads)
		
	N = 1001
	print*, "Before parallel section: N = ", N 

	!$OMP PARALLEL
	N = N + 1
	print*, "Inside parallel section: N = ", N 
	!$OMP END PARALLEL


	print*, "After parallel section: N = ", N

end program main 
