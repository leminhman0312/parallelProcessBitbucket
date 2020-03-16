!! Name: Max Le 
!! Paralle Process Final Project Spring 2019 

!! SOLVING 1D EULER
!! GRID AS FOLLOW
!!  |---------- |--------|------------|
!!  |   i-1     |    i   |     i + 1  |
!!  |---------- |--------|------------|
!!             i-1/2     i+1/2 



!! ASSUME POSITIVE FLUX TO THE RIGHT, NEGATIVE FLUX TO THE LEFT 

!! EACH EDGE (HALF) EXPERIENCES 2 FLUXES: POSITIVE AND NEGATIVE


!! GRID NODE 1 and NMAX are BC.  CALCULATIONS OCCUR 2 to NMAX-2

program SHOCK


    use mpi

    !SOLVING THE SHOCKTUBE PROBLEM 
    implicit none
    integer :: IMAX = 1001
    integer :: I,J,N
    real*8  :: dx,dt,cfl = 0.67, gamma = 1.4
    real*8,dimension(:),allocatable :: a 

    !! FLUX VECTORS
    real*8,dimension(3):: fiph, fimh
    

    real*8,dimension(:,:),allocatable :: fGLOBALp, fGLOBALm 

    !! STATE VECTOR 
    real*8,dimension(:,:),allocatable :: state_vector 


    !! Mach vector 
    real*8,dimension(:),allocatable   :: mach_array 

    !! OTHERS
    real*8 :: t, tmax 
    real*8, dimension(:,:,:),allocatable :: U 
    real*8,dimension(:),allocatable :: x_arr_local 
    real*8 :: xstart, xend, contactpoint, factor   

    integer :: iTrue, iStart!! real index  

    !! MPI VARIABLES 

    integer :: rank, nprocs, ierr 
    integer :: IHIGH, ILOW 
    integer :: nblocks
    real*8 :: min_dt, global_dt
    integer :: leftREQUEST, rightREQUEST, tempREQUEST 


    !! LOCAL VAR 
    integer :: localLOW, localHIGH
    integer :: ghostPoints = 1 
    integer :: maxLocal
    
    !! TIME VARIABLES 
    DOUBLE PRECISION :: tstart, tend, telapse 


    !*********************  INTERFACE   ***********************!
    interface
        subroutine timestep(w_arr,a_arr,delx,Nsize,cfldum,dtdum,factordum)
            implicit none 
            integer,intent(in) :: Nsize 
            real*8, intent(in),DIMENSION(:), ALLOCATABLE   :: a_arr 
            real*8, intent(in),DIMENSION(:,:), ALLOCATABLE :: w_arr
            real*8, intent(in) :: delx, cfldum 
            real*8,intent(out) :: dtdum, factordum
        end subroutine timestep 
    end interface
    
    interface
        subroutine vanleer(fp, fm, w_arr, Nsize, g)
            implicit none
            !! INTENT VARIABLES 
            real*8,intent(in),DIMENSION(:,:), ALLOCATABLE :: w_arr
            real*8,intent(in) :: g
            integer,intent(in) :: Nsize
            real*8,intent(out),DIMENSION(3):: fp, fm 
        end subroutine 
    end interface 






    !!***************** MPI START  *************************** !
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

    !************* INITIALIZING *************************** !
    PRINT*, "INITIALIZING"
    tmax = 2.00
    xstart = -10.00
    xend = 10.00
    contactpoint = 0.00
    dx = (xend-xstart)/DBLE((IMAX-1))

    

    print*, "SETTING NBLOCKS, ILOW, IHIGH" 

    nblocks = FLOOR(REAL(IMAX/nprocs))

    iStart = (nblocks*rank) + ghostPoints 
   
    localLOW = 1+ghostPoints 
    localHIGH = nblocks+ghostPoints 
    maxLocal = localHIGH + ghostPoints 


    ILOW = (rank*nblocks)+1 
    IHIGH = (rank+1)*nblocks

    print*, rank, localLOW, localHIGH 

    print*, "FINISHED SETTING NBLOCKS" 


    print*, "ALLOCATING LOCAL ARRAYS"
    ALLOCATE(a(maxLocal))

    ALLOCATE(fGLOBALm(3,maxLocal))

    ALLOCATE(fGLOBALp(3,maxLocal))

    ALLOCATE(U(3,maxLocal,2))

    ALLOCATE(state_vector(3,maxLocal))

    ALLOCATE(x_arr_local(IMAX+2*ghostPoints))

    ALLOCATE(mach_array(maxLocal))

    print*, "SETTING INITIAL CONDITIONS" 



    !if (rank == 0) then 
        do i = 1,IMAX+2*ghostPoints 
            x_arr_local(i) = xstart + (i-1-ghostPoints)*dx  
        end do
    !end if 
  


    !! SETTING INITIAL RHO, U, P 
    do i = localLOW-ghostPoints ,localHIGH+ghostPoints
        iTrue = i+iStart- ghostPoints 

        ! print*, x_arr_local(iTrue), i , iTrue  

        if (x_arr_local(iTrue) .LE. contactpoint) then
            state_vector(1,i) = 1.00                  !! RHO LEFT 
            state_vector(2,i) = 0.0                   !! U LEFT 
            state_vector(3,i) = 1.00/(gamma-1)        !! P LEFT
        else 
            state_vector(1,i) = 0.125                 !! RHO RIGHT 
            state_vector(2,i) = 0.0                   !! U RIGHT 
            state_vector(3,i) = 0.125/(gamma-1)       !! P RIGHT
        end if
    end do


    
    !! AT TIME = 0 
    do i = localLOW-ghostPoints,localHIGH+ghostPoints 
        U(1,i,1) = state_vector(1,i)                                                    !U(1) = W(1) = RHO
        U(2,i,1) = state_vector(1,i)*state_vector(2,i)                                  !U(2) = W(2) = RHO*U 
        U(3,i,1) = state_vector(3,i)/(gamma-1)  &
            + 0.5*state_vector(1,i)*state_vector(2,i)*state_vector(2,i)                 !ENERGY FORMULA 
        a(i) = sqrt((gamma*state_vector(3,i))/state_vector(1,i))
        mach_array(i) = state_vector(2,i)/(a(i))
    end do


    print*, "FINISHED INITIALIZING" 

    t = 0 


    PRINT*, "SOLVING" 

    tstart = MPI_WTIME()
    do  while (t .LE. tmax)

            !EXIT  

            min_dt = 100 

            do i = localLOW, localHIGH 


                call timestep(state_vector,a,dx,i,cfl,dt,factor)
               
        


                if (dt < min_dt) then 
                    min_dt = dt 
                end if 
                

            end do 

            call MPI_ALLREDUCE(min_dt, global_dt, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD,ierr)
          

            !! COMPUTE FLUXES ON THE EDGE
            do i = localLOW-ghostPoints , localHIGH + ghostPoints  

                call vanleer(fiph, fimh, state_vector, i, gamma)

                fGLOBALm(:,i) = fimh(:) 

                fGLOBALp(:,i) = fiph(:) 
                
            
            end do
            
               


            !!!FINITE DIFFERENCE EQUATION (OLD)
            do i = localLOW,localHIGH 
                do j=1,3
                    U(j,i,2) = U(j,i,1)- (global_dt/dx)*(fGLOBALp(j,i)-fGLOBALp(j,i-1))- &
                            (global_dt/dx)*(fGLOBALm(j,i+1)-fGLOBALm(j,i))        
                end do
            end do

             



            do j = 1,3
                if (rank > 0) then 
                    call MPI_IRECV(U(j,localLOW-1,2),1,MPI_DOUBLE_PRECISION, rank-1,1,MPI_COMM_WORLD,rightREQUEST,ierr)
                    call MPI_Isend(U(j,localLOW ,2),1,MPI_DOUBLE_PRECISION, rank-1,0,MPI_COMM_WORLD,tempREQUEST,ierr)
                end if


                if (rank < nprocs-1) then 
                    call MPI_Irecv(U(j,localHIGH+1,2),1, MPI_DOUBLE_PRECISION, rank+1,0,MPI_COMM_WORLD,leftREQUEST,ierr)
                    call MPI_Isend(U(j,localHIGH,2),1,MPI_DOUBLE_PRECISION,rank+1,1,MPI_COMM_WORLD,tempREQUEST,ierr)
                end if 

                if (rank .NE. 0) then 
                    call MPI_Wait(rightREQUEST,MPI_STATUS_IGNORE,ierr)
                end if 
    
    
                if (rank .NE. nprocs-1) then 
                    call MPI_Wait(leftREQUEST,MPI_STATUS_IGNORE,ierr)
                end if 
            end do 




            





            !! UPDATE U 
            do i = localLOW-ghostPoints,localHIGH+ghostPoints   
                do j = 1,3
                    U(j,i,1) = U(j,i,2)
                end do
            end do   
            !!PULL VARIABLES 
            do i = localLOW-ghostPoints ,localHIGH + ghostPoints  
                state_vector(1,i)   = U(1,i,2)                                     ! PULL RHO
                state_vector(2,i)   = U(2,i,2)/U(1,i,2)                            ! PULL U = RHOU / RHO 
                state_vector(3,i)   = (U(3,i,2)-0.5*U(1,i,2)*(state_vector(2,i)**2))*(gamma-1)  ! PULL PRESSURE 
                a(i)     = sqrt((gamma*state_vector(3,i))/(state_vector(1,i)))                ! UPDATE SOUND SPEED 
            end do   


            !!! BOUNDARY CONDITIONS 

            if (rank == 0) then 
                state_vector(:,1) = state_vector(:,2)
                a(1) = a(2) 
            end if  


            if (rank == nprocs-1) then 
                state_vector(:,localHIGH+1) = state_vector(:,localHIGH) 
                a(localHIGH+1) = a(localHIGH)
            end if  

            t = t + global_dt

            !EXIT  
        end do 
        print*, "FINISHED SOLVING" 

        tend = MPI_WTIME()
        
    
    print*, "WRITE TO FILE" 

    if (nprocs > 1) then
        if (rank == 0) then
            open(unit=  42, file = 'result.txt', action = 'write') 
            do i = localLOW,localHIGH 
                iTrue = i+iStart- ghostPoints 
                write(42,*) x_arr_local(iTrue), state_vector(1,i), state_vector(2,i), state_vector(3,i),state_vector(2,i)/a(i)
            end do
            close(42)
            write(*,*) (tend-tstart)
        end if
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
        
        if (nprocs > 2) then 
            do N = 1,nprocs-2
                if (rank == N) then
                    open(unit=  42, file = 'result.txt', status='old',position='append',action = 'write') 
                    do i = localLOW,localHIGH 
                        iTrue = i+iStart- ghostPoints 
                        write(42,*) x_arr_local(iTrue), state_vector(1,i), state_vector(2,i), state_vector(3,i),&
                                state_vector(2,i)/a(i)
                    end do
                    close(42)
                    write(*,*) (tend-tstart) 
                end if 
                call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
            end do
        end if
        
        if (rank == nprocs-1) then
            open(unit=  42, file = 'result.txt', status='old',position='append',action = 'write') 
            do i = localLOW,localHIGH 
                iTrue = i+iStart- ghostPoints 
                write(42,*) x_arr_local(iTrue), state_vector(1,i), state_vector(2,i), state_vector(3,i),state_vector(2,i)/a(i)
            end do
            close(42)
            write(*,*) (tend-tstart)            
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr) 


    else  !! FOR 1 core case 
        open(unit=  42, file = 'result.txt', action = 'write') 
        do i = localLOW,localHIGH 
            write(42,*) x_arr_local(i), state_vector(1,i), state_vector(2,i), state_vector(3,i),state_vector(2,i)/a(i)
        end do 
        close(42)
        write(*,*) (tend-tstart) 
    end if 

    call MPI_FINALIZE(ierr)
    print*, "DEALLOCATE ARRAYS "

   
    !!!! DEALLOCATE ARRAY 
    DEALLOCATE(a)
    DEALLOCATE(fGLOBALm)
    DEALLOCATE(fGLOBALp)
    DEALLOCATE(U)
    DEALLOCATE(state_vector)
    DEALLOCATE(x_arr_local)    
    DEALLOCATE(mach_array)


end program SHOCK




subroutine timestep(w_arr,a_arr,dxdum,Nsize,cfldum,dtdum,factordum)
    implicit none 
    integer :: i 
    real*8, intent(in),dimension(:,:),ALLOCATABLE  :: w_arr
    real*8, intent(in),dimension(:),ALLOCATABLE    :: a_arr
    real*8, intent(in)  :: cfldum,dxdum 
    integer,intent(in)  :: Nsize 
    real*8, intent(out) :: dtdum
    real*8,intent(out) :: factordum


    real*8 :: max_wave_speed = 0.0

    real*8,dimension(:),allocatable :: local_max 

    allocate(local_max(Nsize))

    do i = 2,Nsize-1
        local_max(i) = ABS(w_arr(2,i)) + a_arr(i) 
    end do 

    max_wave_speed = MAXVAL(local_max)

    dtdum = (cfldum*dxdum)/max_wave_speed

    factordum = max_wave_speed*(dtdum/dxdum) 


end subroutine timestep 







subroutine vanleer(fp, fm, w_arr,Nsize, g)

    implicit none


    !! INTENT VARIABLES 
    real*8,intent(in),DIMENSION(:,:), ALLOCATABLE :: w_arr
    real*8,intent(in) :: g
    integer,intent(in) :: Nsize
    real*8,intent(out),DIMENSION(3) :: fp, fm 




    !! LOCAL VARIABLES
    integer :: i  
    real*8 :: machNumber 
    real*8 :: faplus = 0.0, fbplus = 0.0, faminus = 0.0, fbminus = 0.0 

    real*8 :: ai, Ei  


    i = Nsize 

    !! Calculate ai, pi 

    ai = sqrt((g*w_arr(3,i))/(w_arr(1,i)))

    !! Build variables 
    machNumber = w_arr(2,i)/ai                       !! mach = u/ai 

    if (abs(machNumber) .LE. 1) then 
        faplus = 0.25*w_arr(1,i)*ai*(machNumber+1)**2 
        fbplus = (g-1.0)*machNumber*ai + 2*ai 

        faminus = -0.25*w_arr(1,i)*ai*(machNumber-1)**2 
        fbminus = -(g-1.0)*machNumber*ai - 2*ai 

        !! Building up Van Leer flux vector 
        fp(1) = faplus 
        fp(2) = faplus*fbplus*(1.0/g)
        fp(3) = faplus*fbplus**2*(1.0/(2.0*(g**2-1)))


        fm(1) = faminus 
        fm(2) = faminus*fbminus*(1.0/g)
        fm(3) = faminus*fbminus**2*(1.0/(2.0*(g**2-1)))

    else 
        Ei = w_arr(3,i)/(g-1.0) + 0.5*w_arr(1,i)*w_arr(2,i)**2
        if (w_arr(2,i) > 0) then 
            fp(1) = w_arr(1,i)*w_arr(2,i)
            fp(2) = w_arr(1,i)*w_arr(2,i)**2 + w_arr(3,i) 
            fp(3) = (Ei+w_arr(3,i))*w_arr(2,i)

            fm(1) = 0.0 
            fm(2) = 0.0
            fm(3) = 0.0


        else 

            fm(1) = w_arr(1,i)*w_arr(2,i)
            fm(2) = w_arr(1,i)*w_arr(2,i)**2 + w_arr(3,i) 
            fm(3) = (Ei+w_arr(3,i))*w_arr(2,i)

            fp(1) = 0.0 
            fp(2) = 0.0 
            fp(3) = 0.0
        
        end if
         !! Reset to zero for next run 
        machNumber = 0.0
        faplus = 0.0
        fbplus = 0.0 
        faminus = 0.0 
        fbminus = 0.0 
    end if 

   



end subroutine 



!subroutine output


    !if (nprocs > 1) then
        !if (rank == 0) then
            !open(unit=  42, file = 'result.txt', action = 'write') 
            !do i = localLOW,localHIGH 
                !iTrue = i+iStart- ghostPoints 
                !write(42,*) x_arr_local(iTrue), state_vector(1,i), state_vector(2,i), state_vector(3,i),state_vector(2,i)/a(i)
            !end do
            !close(42)
        !end if
        
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
        
        !if (nprocs > 2) then 
            !do N = 1,nprocs-2
                !if (rank == N) then
                    !open(unit=  42, file = 'result.txt', status='old',position='append',action = 'write') 
                    !do i = localLOW,localHIGH 
                        !iTrue = i+iStart- ghostPoints 
                        !write(42,*) x_arr_local(iTrue), state_vector(1,i), state_vector(2,i), state_vector(3,i),&
                                !state_vector(2,i)/a(i)
                    !end do
                    !close(42)
                !end if 
                !call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
            !end do 
        !end if
        
        !if (rank == nprocs-1) then
            !open(unit=  42, file = 'result.txt', status='old',position='append',action = 'write') 
            !do i = localLOW,localHIGH 
                !iTrue = i+iStart- ghostPoints 
                !write(42,*) x_arr_local(iTrue), state_vector(1,i), state_vector(2,i), state_vector(3,i),state_vector(2,i)/a(i)
            !end do
            !close(42)            
        !end if
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr) 


    !else  !! FOR 1 core case 
        !open(unit=  42, file = 'result.txt', action = 'write') 
        !do i = localLOW,localHIGH 
            !write(42,*) x_arr_local(i), state_vector(1,i), state_vector(2,i), state_vector(3,i),state_vector(2,i)/a(i)
        !end do 
        !close(42)
    !end if 


!end subroutine output 















            !do j = 1,3
                !if (rank >  0) then
                    !call MPI_ISEND(U(j,ILOW,2),1,MPI_REAL,rank-1,0,MPI_COMM_WORLD,ierr) 
                !end if
                
                !if (rank < nprocs -1 ) then 
                    !call MPI_IRECV(U(j,IHIGH+1,2),1,MPI_REAL,rank+1,0,MPI_COMM_WORLD,ierr)
                !end if 

                !if (rank < nprocs-1) then 

                    !call MPI_ISEND(U(j,IHIGH,2),1,MPI_REAL,rank+1,1,MPI_COMM_WORLD,ierr) 
                !end if 

                !if (rank > 0) then 

                    !call MPI_IRECV(U(j,ILOW-1,2),1,MPI_REAL,rank+1,1,MPI_COMM_WORLD,ierr)
                !end if    

            !end do

