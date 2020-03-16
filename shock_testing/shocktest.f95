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


    ! use mpi

    !SOLVING THE SHOCKTUBE PROBLEM 
    implicit none
    integer :: IMAX = 10 
    integer :: I,J,M,N
    real*8  :: dx,dt,cfl = 0.67, gamma = 1.4
    real*8,dimension(:),allocatable :: a, H

    !! FLUX VECTORS
    real*8,dimension(:,:),allocatable :: fiph, fimh


    !! STATE VECTOR 
    real*8,dimension(:,:),allocatable :: state_vector


    !! Mach vector 
    real*8,dimension(:),allocatable   :: mach_array 

    !! OTHERS
    real*8 :: t, tmax 
    real*8, dimension(:,:,:),allocatable :: U 
    real*8,dimension(:),allocatable :: x_arr 
    real*8 :: xstart, xend, contactpoint, factor   

    !! FILES 
    character(len=50) :: filename 
    integer :: count


    !! MPI VARIABLES 

    integer :: rank, nprocs, ierr 
    integer :: IHIGH, ILOW 
    integer :: nblocks 





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
        subroutine steger(fph_arr, fmh_arr, w_arr, a_arr, h_arr, Nsize, g) 
            implicit none
            integer,intent(in) :: Nsize 
            real*8,intent(in)  :: g
            real*8, intent(in),DIMENSION(:,:), ALLOCATABLE :: w_arr      
            real*8, intent(out),DIMENSION(:,:), ALLOCATABLE :: fph_arr, fmh_arr  
            real*8, intent(in),DIMENSION(:), ALLOCATABLE   :: a_arr, h_arr
        end subroutine steger
    end interface 

    interface
        subroutine vanleer(fp, fm, w_arr, a_arr, h_arr, Nsize, g,mach_arr_dum)
            implicit none
            !! INTENT VARIABLES 
            real*8,intent(in),DIMENSION(:,:), ALLOCATABLE :: w_arr
            real*8,intent(in),DIMENSION(:), ALLOCATABLE   :: a_arr, h_arr 
            real*8,intent(in) :: g
            integer,intent(in) :: Nsize
            real*8,intent(out),DIMENSION(:,:),ALLOCATABLE :: fp, fm 
            real*8,intent(out),DIMENSION(:),ALLOCATABLE :: mach_arr_dum  
        end subroutine 
    end interface 




    !**************  ALLOCATE ARRAYS  ************************ !

    ALLOCATE(a(IMAX))

    ALLOCATE(H(IMAX))


    ALLOCATE(fiph(3,IMAX))

    ALLOCATE(fimh(3,IMAX))

    ALLOCATE(U(3,IMAX,2))

    ALLOCATE(state_vector(3,IMAX))

    ALLOCATE(x_arr(IMAX))

    ALLOCATE(mach_array(IMAX))


    !************* INITIALIZING *************************** !

    PRINT*, "INITIALIZING"
    tmax = 2.00
    xstart = -10.00
    xend = 10.00
    contactpoint = 0.0 
    x_arr(1) = xstart
    x_arr(IMAX) = xend 
    dx = (xend-xstart)/DBLE((IMAX-1))


    do i = 2,IMAX-1
        x_arr(i) = x_arr(i-1) + dx 
    end do 



    !! SETTING INITIAL RHO, U, P 
    do i = 1,IMAX
        if (x_arr(i) .LE. contactpoint) then
            state_vector(1,i) = 1.00        !! RHO LEFT 
            state_vector(2,i) = 0.0         !! U LEFT 
            state_vector(3,i) = 1.00/(gamma-1)      !! P LEFT
        else 
            state_vector(1,i) = 0.125        !! RHO RIGHT 
            state_vector(2,i) = 0.0         !! U RIGHT 
            state_vector(3,i) = 0.125/(gamma-1)        !! P RIGHT 
        end if
    end do



    ! AT TIME = 0 
    do i = 2,IMAX-1
        U(1,i,1) = state_vector(1,i)                                                    !U(1) = W(1) = RHO
        U(2,i,1) = state_vector(1,i)*state_vector(2,i)                                  !U(2) = W(2) = RHO*U 
        U(3,i,1) = state_vector(3,i)/(gamma-1)  &
            + 0.5*state_vector(1,i)*state_vector(2,i)*state_vector(2,i)                 !ENERGY FORMULA 
        a(i) = sqrt((gamma*state_vector(3,i))/state_vector(1,i))
        H(i) = (U(3,i,1) + state_vector(3,i))/state_vector(1,i)
        mach_array(i) = state_vector(2,i)/(a(i))
    end do

    !! ZERO GRADIENT BC
    H(1) = H(2) 
    a(1) = a(2) 
    H(IMAX) = H(IMAX-1)
    a(IMAX) = a(IMAX-1)

    mach_array(1) = mach_array(2)
    mach_array(IMAX) = mach_array(IMAX-1)


    t = 0.0


    !!***************** SOLVING *************************** !
  

    count = 0 
    
    PRINT*, "SOLVING"
    do  while (t .LE. tmax)


        write (filename, "(A6,I0.4,A4)") "result", count, ".txt"

        !print*, filename

        

        !!! CALCULATE TIMESTEP BASE ON CFL   
        call timestep(state_vector,a,dx,IMAX,cfl,dt,factor)

        !print*, dt , factor 

        !!! CALCULATE FLUXES
        call vanleer(fiph, fimh, state_vector, a, H, IMAX, gamma,mach_array)    


        
        !!!FINITE DIFFERENCE EQUATION
        do i = 2,IMAX-1
            do j=1,3
                U(j,i,2) = U(j,i,1)- (dt/dx)*(fiph(j,i)-fiph(j,i-1))- (dt/dx)*(fimh(j,i+1)-fimh(j,i))        
            end do
        end do   
        
       !!PULL VARIABLES 
        do i = 2,IMAX-1
            do j=1,3
                state_vector(1,i)   = U(1,i,2)                                     ! PULL RHO
                state_vector(2,i)   = U(2,i,2)/U(1,i,2)                            ! PULL U = RHOU / RHO 
                state_vector(3,i)   = (U(3,i,2)-0.5*U(1,i,2)*(state_vector(2,i)**2))*(gamma-1)  ! PULL PRESSURE 
                H(i)     = (U(3,i,1)+state_vector(3,i))/state_vector(1,i)                     ! UPDATE ENTHALPY
                a(i)     = sqrt((gamma*state_vector(3,i))/(state_vector(1,i)))                ! UPDATE SOUND SPEED 
            end do
        end do   


        !!! BOUNDARY CONDITIONS 
        do j = 1,3 
            state_vector(j,1) = state_vector(j,2) 
            state_vector(j,IMAX) = state_vector(j,IMAX-1)
        end do


        !!! UPDATE ENTIRE U 

        do i = 2,IMAX-1
            do j = 1,3
                U(j,i,1) = U(j,i,2)
            end do
        end do   

        !! UPDATE TIME 
        t = t + dt 
        count = count + 1

    end do  


    print*, "WRITING TO FILE" 

    !print*, "WRITE TO FILE" 
    open(unit=  42, file = 'result.txt', action = 'write')

    do i = 1,IMAX 
        write(42,*) x_arr(i), state_vector(1,i), state_vector(2,i), state_vector(3,i), mach_array(i) 
    end do 


    print*, "DEALLOCATE ARRAYS "



    !!!! CLOSE THE FILE
    close(42)

    !! DEALLOCATE ARRAY 
    DEALLOCATE(a)
    DEALLOCATE(H)
    DEALLOCATE(fiph)
    DEALLOCATE(fimh)
    DEALLOCATE(U)
    DEALLOCATE(state_vector)

    print*, "PROGRAM FINISHED" 


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







subroutine vanleer(fp, fm, w_arr, a_arr, h_arr, Nsize, g,mach_arr_dum)

    implicit none


    !! INTENT VARIABLES 
    !real*8,intent(out),DIMENSION(:,:), ALLOCATABLE :: fph_arr, fmh_arr
    real*8,intent(in),DIMENSION(:,:), ALLOCATABLE :: w_arr
    real*8,intent(in),DIMENSION(:), ALLOCATABLE   :: a_arr, h_arr 
    real*8,intent(in) :: g
    integer,intent(in) :: Nsize
    real*8,intent(out),DIMENSION(:,:),ALLOCATABLE :: fp, fm 
    real*8,intent(out),DIMENSION(:),ALLOCATABLE :: mach_arr_dum  



    !! LOCAL VARIABLES
    integer :: i, j 
    real*8  :: machNumber 
    real*8 :: faplus = 0.0, fbplus = 0.0, faminus = 0.0, fbminus = 0.0 

    real*8 :: ai, Ei  

    ALLOCATE(fp(3,Nsize))
    ALLOCATE(fm(3,Nsize))

    ALLOCATE(mach_arr_dum(Nsize))

    !ALLOCATE(fph_arr(3,Nsize))
    !ALLOCATE(fmh_arr(3,Nsize))





    !DO i = 2,Nsize-1

    DO i = 1,Nsize
        !! Calculate ai, pi 

    
        ai = sqrt((g*w_arr(3,i))/(w_arr(1,i)))



        !! Build variables 
        machNumber = w_arr(2,i)/ai                       !! mach = u/ai 

        mach_arr_dum(i) = machNumber 
        


        if (abs(machNumber) .LE. 1) then  

            faplus = 0.25*w_arr(1,i)*ai*(machNumber+1)**2 
            fbplus = (g-1.0)*machNumber*ai + 2*ai 

            faminus = -0.25*w_arr(1,i)*ai*(machNumber-1)**2 
            fbminus = -(g-1.0)*machNumber*ai - 2*ai 

            !! Building up Van Leer flux vector 
            fp(1,i) = faplus 
            fp(2,i) = faplus*fbplus*(1.0/g)
            fp(3,i) = faplus*fbplus**2*(1.0/(2.0*(g**2-1)))


            fm(1,i) = faminus 
            fm(2,i) = faminus*fbminus*(1.0/g)
            fm(3,i) = faminus*fbminus**2*(1.0/(2.0*(g**2-1)))

        else
            Ei = w_arr(3,i)/(g-1.0) + 0.5*w_arr(1,i)*w_arr(2,i)**2
            if (w_arr(2,i) > 0) then 
                fp(1,i) = w_arr(1,i)*w_arr(2,i)
                fp(2,i) = w_arr(1,i)*w_arr(2,i)**2 + w_arr(3,i) 
                fp(3,i) = (Ei+w_arr(3,i))*w_arr(2,i)

                fm(1,i) = 0.0 
                fm(2,i) = 0.0
                fm(3,i) = 0.0


            else 

                fm(1,i) = w_arr(1,i)*w_arr(2,i)
                fm(2,i) = w_arr(1,i)*w_arr(2,i)**2 + w_arr(3,i) 
                fm(3,i) = (Ei+w_arr(3,i))*w_arr(2,i)

                fp(1,i) = 0.0 
                fp(2,i) = 0.0 
                fp(3,i) = 0.0
            end if 
        end if 

        !! Reset to zero for next run 
        machNumber = 0.0
        faplus = 0.0
        fbplus = 0.0 
        faminus = 0.0 
        fbminus = 0.0 

    END DO 


    !DO i = 3,Nsize-2 
        !do j = 1,3
            !fph_arr(j,i) = 0.0
            !fmh_arr(j,i) = 0.0

            !fph_arr(j,i) = fp(j,i)+ fm(j,i+1)
            !fmh_arr(j,i) = fp(j,i-1)+fm(j,i)
        !end do
    !end do



end subroutine 















