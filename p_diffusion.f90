!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program diffusion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none

  ! Dimensions of the system
  real(8) :: L1 = 10.d0, L2=10.d0, L3=20.d0
  ! Grid size
  integer :: n1 = 48, n2 = 48, n3 = 96
  ! time step for time integration
  real(8) :: dt = 2.d-3
  ! number of time steps
  integer :: nstep = 10
  ! Radius of diffusion channel
  real(8) :: rad_diff = 0.7d0
  ! Radius of starting concentration
  real(8) :: rad_conc = 0.6d0

  real :: start, finish

  real(8), dimension(:,:,:) , allocatable :: diffusivity, conc, dconc

  real(8), dimension(:,:,:), allocatable    :: aux1, aux2


  integer :: status, i1, i2, i3, ipol, istep
  real(8) :: f1conc, f2conc, f3conc, f1diff, f2diff, f3diff, fac, ss
  real(8) :: x1,x2,x3, rr, r2mean
  logical :: direct_to_reciprocal

  fac = L1*L2*L3/(n1*n2*n3)
  !
  ! Allocate arrays
  !
  ! diffusivity is the diffusion coefficient in direct space
  allocate(diffusivity(0:n1-1,0:n2-1,0:n3-1), stat=status)
  call errore('diffusion','error allocating array diffusivity',status)

  ! conc is the concentration in direct space
  allocate(conc(0:n1-1,0:n2-1,0:n3-1), stat=status)
  call errore('diffusion','error allocating array conc',status)

  ! dconc is the time derivative of the concentration in direct space
  allocate(dconc(0:n1-1,0:n2-1,0:n3-1), stat=status)
  call errore('diffusion','error allocating array dconc',status)

  ! aux1 is an auxiliary function in direct space 
  allocate(aux1(0:n1-1,0:n2-1,0:n3-1), stat=status)
  call errore('diffusion','error allocating array aux1',status)

  ! aux2 is an auxiliary function in direct space 
  allocate(aux2(0:n1-1,0:n2-1,0:n3-1), stat=status)
  call errore('diffusion','error allocating array aux2',status)
  ! 
  ! Define the diffusivity inside the system and 
  ! the starting concentration
  !
  ! ss is to integrate (and normalize) the concentration
  !
  ss = 0.d0

  do i3=0,n3-1
     x3 = L3*real(i3,8)/n3
     f3diff = exp(-((x3-0.5d0*L3)/rad_diff)**2)
     f3conc = exp(-((x3-0.5d0*L3)/rad_conc)**2)
     do i2=0,n2-1
        x2 = L2*real(i2,8)/n2
        f2diff = exp(-((x2-0.5d0*L2)/rad_diff)**2)
        f2conc = exp(-((x2-0.5d0*L2)/rad_conc)**2)
        do i1=0,n1-1
           x1 = L1*real(i1,8)/n1
           f1diff = exp(-((x1-0.5d0*L1)/rad_diff)**2)
           f1conc = exp(-((x1-0.5d0*L1)/rad_conc)**2)

           diffusivity(i1,i2,i3) = max(f1diff*f2diff,f2diff*f3diff)
           conc(i1,i2,i3) = f1conc*f2conc*f3conc
           ss = ss+conc(i1,i2,i3)
        enddo
     enddo
  enddo

  call plot_data_2d('diffusivity',n1, n2, n3, 1, diffusivity) 
  call plot_data_2d('diffusivity',n1, n2, n3, 2, diffusivity) 
  call plot_data_2d('diffusivity',n1, n2, n3, 3, diffusivity) 

  ! Now normalize the concentration
  ss = 1.d0/(ss*fac)
  conc(:,:,:) = conc(:,:,:)*ss


  ! Now everything is defined: system size, diffusivity inside the system, and
  ! the starting concentration
  !
  ! Start the dynamics
  !
  call cpu_time(start)
  do istep=1,nstep
     dconc(:,:,:) = 0.d0
     do ipol=1,3
        call derivative(n1,n2,n3,L1,L2,L3,ipol,conc,aux1)
        aux1(:,:,:) = aux1(:,:,:)*diffusivity(:,:,:)
        
        call derivative(n1,n2,n3,L1,L2,L3,ipol,aux1,aux2)
  ! summing up contributions from the three spatial directions      
        dconc(:,:,:) = dconc(:,:,:) + aux2(:,:,:)
     enddo

     conc(:,:,:) = conc(:,:,:) + dt*dconc(:,:,:)
    
     if(mod(istep,30).eq.1) then
        ! Check the normalization of conc
        ss = 0.d0
        r2mean = 0.d0
        do i3=0,n3-1
           x3 = L3*real(i3,8)/n3
           do i2=0,n2-1
              x2 = L2*real(i2,8)/n2
              do i1=0,n1-1
                 x1 = L1*real(i1,8)/n1
                 rr = (x1-0.5d0*L1)**2 + (x2-0.5d0*L2)**2 + (x3-0.5d0*L3)**2 
                 ss = ss+conc(i1,i2,i3)
                 r2mean = r2mean + conc(i1,i2,i3)*rr
              enddo
           enddo
        enddo
        ss = ss*fac
        r2mean = r2mean*fac

        call cpu_time(finish)
        write(6,*) istep, r2mean, ss ,  ' elapsed time per iteration ', (finish-start)/istep
        call plot_data_2d('concentration',n1, n2, n3, 2, conc) 
        call plot_data_1d('1d_conc',n1, n2, n3, 3, conc) 

     endif

  enddo

  deallocate(aux1)
  deallocate(aux2)
  deallocate(diffusivity)
  deallocate(conc)
  deallocate(dconc)

end program diffusion
  
