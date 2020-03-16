!! Name: Max Le
!! Parallel Process Final Project Spring 2019

!! SOLVING 1D EULER
!! GRID AS FOLLOW
!!  |---------- |--------|------------|
!!  |   i-1     |    i   |     i + 1  |
!!  |---------- |--------|------------|
!!             i-1/2     i+1/2

!! ASSUME POSITIVE FLUX TO THE RIGHT, NEGATIVE FLUX TO THE LEFT

!! EACH EDGE (HALF) EXPERIENCES 2 FLUXES: POSITIVE AND NEGATIVE

!! GRID NODE 1 and NMAX are BC.  CALCULATIONS OCCUR 2 to NMAX-2



!!**********************  MODULE FOR ALL VARIABLES **************************** !! 

module globalvar

    integer :: IMAX = 10000
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
    real*8 :: t, tmax , ai
    real*8, dimension(:,:,:),allocatable :: U
    real*8,dimension(:),allocatable :: x_arr_local
    real*8 :: xstart, xend, contactpoint

    integer :: iTrue, iStart!! real index

    !! MPI VARIABLES

    integer :: rank, nprocs, ierr
    integer :: nblocks, rem
    real*8 :: min_dt, global_dt
    integer :: leftREQUEST, rightREQUEST, tempREQUEST

    integer :: sendLEFT,sendRIGHT, recvLEFT, recvRIGHT 

    !! LOCAL VAR
    integer :: localLOW, localHIGH
    integer :: ghostPoints = 2
    integer :: maxLocal, minLocal



    !! SEND/RECV MATRIX
    DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: s1mat, r1mat, s2mat, r2mat

    !! FLAG 

    integer :: flagdum,npointsdum  

    !*********************  INTERFACE   ***********************!
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


end module globalvar 



!!**********************  END MODULE  **************************************** !! 













!!**************************************************************************** !! 
!!**********************  MAIN PROGRAM *************************************** !! 
!!**************************************************************************** !! 

program SHOCK
    use globalvar 
    include 'mpif.h' 


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

    rem = IMAX - nblocks*nprocs


    !! starting REAL index on each chunk
    iStart = (nblocks*rank) + ghostPoints       

    localLOW = 1+ghostPoints
    !localLOW = ghostPoints
    localHIGH = nblocks+ghostPoints
    maxLocal = localHIGH + ghostPoints
    minLocal = localLOW-ghostPoints

    ! adding the remaining point
    if(rem /= 0) then
        if(rank == nprocs-1) then
            localHIGH = localHIGH + 1
            maxlocal  = maxlocal + 1
        end if
    end if


123 FORMAT("I am rank: ",i2," and ILO=",i4," and IHI=",i4,&
    " minlocal=", &i4, " maxlocal=",i4)
    WRITE(*,123) rank,localLOW,localHIGH,minlocal, maxlocal


    print*, "FINISHED SETTING NBLOCKS"


    print*, "ALLOCATING LOCAL ARRAYS"
    ALLOCATE(a(maxLocal))

    ALLOCATE(fGLOBALm(3,maxLocal))

    ALLOCATE(fGLOBALp(3,maxLocal))

    ALLOCATE(U(3,maxLocal,2))

    ALLOCATE(state_vector(3,maxLocal))

    ALLOCATE(x_arr_local(IMAX+2*ghostPoints))

    ALLOCATE(mach_array(maxLocal))

    ALLOCATE(s1mat(3,ghostPoints))
    ALLOCATE(s2mat(3,ghostPoints))
    ALLOCATE(r1mat(3,ghostPoints))
    ALLOCATE(r2mat(3,ghostPoints))

    print*, "SETTING INITIAL CONDITIONS"



    do i = 1,IMAX+2*ghostPoints
        x_arr_local(i) = xstart + (i-1-ghostPoints)*dx
    end do



    !! SETTING INITIAL RHO, U, P
    do i = minLocal , maxLocal
        iTrue = i+iStart- ghostPoints

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
    do i = minLocal,maxLocal
        U(1,i,1) = state_vector(1,i)                                                    
        U(2,i,1) = state_vector(1,i)*state_vector(2,i)                                  
        U(3,i,1) = state_vector(3,i)/(gamma-1)  &
            + 0.5*state_vector(1,i)*state_vector(2,i)*state_vector(2,i)                 
        a(i) = sqrt((gamma*state_vector(3,i))/state_vector(1,i))
        mach_array(i) = state_vector(2,i)/(a(i))
    end do


    flagdum = 0
    print*, "FINISHED INITIALIZING"



    PRINT*, "SOLVING"
    t = 0
    do  while (t .LE. tmax)

            !!! GET GLOBAL DT 
            min_dt = 100

            do i = localLOW, localHIGH

                ai = sqrt(gamma*state_vector(3,i)/state_vector(1,i))
                dt = cfl*dx/(abs(state_vector(2,i)) + ai)

                if (dt < min_dt) then
                    min_dt = dt
                end if

            end do

            call MPI_ALLREDUCE(min_dt, global_dt, 1, MPI_DOUBLE_PRECISION, &
                MPI_MIN, MPI_COMM_WORLD,ierr)



            !!! CALCULATE FLUX AND UPDATE U 
            call update(flagdum)  

            !! UPDATE U
            do i = localLOW-ghostPoints,localHIGH+ghostPoints
                do j = 1,3
                    U(j,i,1) = U(j,i,2)
                end do
            end do
            

            !!PULL VARIABLES
            do i = minLocal , maxLocal
                 ! PULL RHO
                state_vector(1,i)   = U(1,i,2)     
                ! PULL U = RHOU / RHO                               
                state_vector(2,i)   = U(2,i,2)/U(1,i,2)            
                ! PULL PRESSURE                
                state_vector(3,i)   = (U(3,i,2)-0.5*U(1,i,2)*&
                    (state_vector(2,i)**2))*(gamma-1)  
                ! UPDATE SOUND SPEED
                a(i)     = sqrt((gamma*state_vector(3,i))/(state_vector(1,i)))               
            end do

            !! UPDATE TIME STEP 
            t = t + global_dt
            print*, "t = ", t 

    end do
    
    print*, "FINISHED SOLVING"

    call output 

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
    DEALLOCATE(s1mat)
    DEALLOCATE(s2mat)
    DEALLOCATE(r1mat)
    DEALLOCATE(r2mat)



end program SHOCK



!!**************************************************************************** !! 
!!**********************  FINISHED MAIN PROGRAM ****************************** !! 
!!**************************************************************************** !! 














!!**********************  SUBROUTINE FOR FLUX  ******************************* !! 


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




!!**********************  END SUBROUTINE FOR FLUX  *************************** !! 
























!!**********************  SUBROUTINE FOR SEND/RECV *************************** !! 


subroutine bc_recv 

    use globalvar
    include 'mpif.h' 

    if (rank > 0) then
        call MPI_Wait(sendLEFT,MPI_STATUS_IGNORE,ierr)
        call MPI_WAIT(recvLEFT,MPI_STATUS_IGNORE,ierr) 
    end if


    if (rank < nprocs-1) then
        call MPI_Wait(sendRIGHT,MPI_STATUS_IGNORE,ierr)
        call MPI_Wait(recvRIGHT,MPI_STATUS_IGNORE,ierr) 

    end if

    !! SETUP RECV MATRIX
    do j = 1,ghostPoints
        U(:,localLOW-(ghostPoints-j+1),2) = r1mat(:,j)
        U(:,localHIGH+j,2) = r2mat(:,j)
        
    end do

end subroutine bc_recv  

subroutine bc_send
    use globalvar
    include 'mpif.h' 
    

    !! SEND/RECV
    if (rank > 0 ) then
        call MPI_Irecv(r1mat,ghostPoints*3,MPI_DOUBLE_PRECISION,rank-1,0,&
            MPI_COMM_WORLD,recvLEFT,ierr)

    end if


    if (rank <nprocs-1) then
        call MPI_Irecv(r2mat,ghostPoints*3,MPI_DOUBLE_PRECISION,rank+1,1,&
            MPI_COMM_WORLD,recvRIGHT,ierr)

    end if


    !! SETUP SEND MATRIX
    do j = 1,ghostPoints
        s1mat(:,j) = U(:,localHIGH-(ghostPoints-j),2)
        s2mat(:,j) = U(:,localLOW+(j-1),2)

    end do



    if (rank > 0 ) then
        call MPI_Isend(s2mat,ghostPoints*3,MPI_DOUBLE_PRECISION,rank-1,1,&
            MPI_COMM_WORLD,sendLEFT,ierr)
    end if


    if (rank <nprocs-1) then
        call MPI_Isend(s1mat,ghostPoints*3,MPI_DOUBLE_PRECISION,rank+1,0,&
            MPI_COMM_WORLD,sendRIGHT,ierr)
    end if



end subroutine bc_send





!!**********************  END SUBROUTINE FOR SEND/RECV ************************ !! 






























!!**********************  SUBROUTINE FOR OUTPUT  ***************************** !! 


subroutine output 
    use globalvar 
    include 'mpif.h' 


    print*, "WRITE TO FILE"

    if (nprocs > 1) then
        if (rank == 0) then
            open(unit=  42, file = 'result.txt', action = 'write')
            do i = localLOW,localHIGH
                iTrue = i+iStart- ghostPoints
                write(42,*) x_arr_local(iTrue), state_vector(1,i), &
                state_vector(2,i), state_vector(3,i),state_vector(2,i)/a(i)
            end do
            close(42)
        end if

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if (nprocs > 2) then
            do N = 1,nprocs-2
                if (rank == N) then
                    open(unit=  42, file = 'result.txt', status='old',&
                        position='append',action = 'write')
                    do i = localLOW,localHIGH
                        iTrue = i+iStart- ghostPoints
                        write(42,*) x_arr_local(iTrue), state_vector(1,i), &
                            state_vector(2,i), state_vector(3,i),&
                            state_vector(2,i)/a(i)
                    end do
                    close(42)
                end if
                call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            end do
        end if

        if (rank == nprocs-1) then
            open(unit=  42, file = 'result.txt', status='old',&
                position='append',action = 'write')
            do i = localLOW,localHIGH
                iTrue = i+iStart- ghostPoints
                write(42,*) x_arr_local(iTrue), state_vector(1,i), &
                    state_vector(2,i), state_vector(3,i),state_vector(2,i)/a(i)
            end do
            close(42)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)


    else  !! FOR 1 core case
        open(unit=  42, file = 'result.txt', action = 'write')
        do i = localLOW,localHIGH
            write(42,*) x_arr_local(i), state_vector(1,i), state_vector(2,i), &
                state_vector(3,i),state_vector(2,i)/a(i)
        end do
        close(42)
    end if

end subroutine output 




!!**********************  END SUBROUTINE FOR OUTPUT  ************************** !! 








































!!**********************  SUBROUTINE FOR UPDATE  ***************************** !! 


subroutine update(temp) 
    use globalvar 
    include 'mpif.h' 

    integer, intent(inout) :: temp 


    !! FLUX CALCULATIONS 
    do i = minLocal, maxLocal      
        call vanleer(fiph, fimh, state_vector, i, gamma)
        fGLOBALm(:,i) = fimh(:)
        fGLOBALp(:,i) = fiph(:)
    end do


    ! FINITE DIFFERENCE 
    do i = localLOW,localHIGH 
        do j=1,3
            U(j,i,2) = U(j,i,1)- (global_dt/dx)*(fGLOBALp(j,i)-fGLOBALp(j,i-1))- &
                (global_dt/dx)*(fGLOBALm(j,i+1)-fGLOBALm(j,i))
        end do
    end do


    if (temp == 0) then
        print*, "DOING SEND/RECV" 
        call bc_send 
        call bc_recv
        temp = ghostPoints 
    end if 

    temp = temp -1 
    call boundary

end subroutine update 


!!**********************  END SUBROUTINE FOR UPDATE  ************************** !! 













!!**********************  SUBROUTINE FOR BOUNDARY  *************************** !! 

subroutine boundary 
    use globalvar 
    include 'mpif.h' 

    !!! BOUNDARY CONDITIONS
    if (rank == 0) then
        do j = 1,ghostPoints
            U(:,localLOW-j,2) = U(:,localLOW-j+1,2)
        end do
    end if

    if (rank == nprocs-1) then
        do j = 1,ghostPoints
            U(:,localHIGH+j,2) = U(:,localHIGH+j-1,2)
        end do
    end if
end subroutine boundary 

!!**********************  END SUBROUTINE FOR BOUNDARY  *********************** !! 
