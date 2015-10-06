!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine derivative(n1,n2,n3,L1,L2,L3,ipol,data,deriv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !
  ! Calculate the derivative in direction ipol of the array 'data'
  !

  implicit none

  integer, intent(in) :: n1,n2,n3
  real(8), intent(in) :: L1,L2,L3
  integer, intent(in) :: ipol
  real(8), dimension(0:n1-1, 0:n2-1, 0:n3-1), intent(in) :: data 
  real(8), dimension(0:n1-1, 0:n2-1, 0:n3-1), intent(out) :: deriv

  real(8), parameter :: pi = 3.14159265358979323846d0
  
  logical :: use_3d_fft
  complex(8), dimension(:,:,:), allocatable :: aux
  complex(8), dimension(:), allocatable :: aux_1d
  real(8), dimension(:), allocatable :: aux_real_1d
  integer :: status
  real(8) :: G
  integer :: i,i1,i2,i3

  ! 
  ! Here, you can find to versions of calculating the derivative
  ! of a 3D-array in the direction ipol
  ! The first implementation, used if use_3d_fft=.true., will perform
  ! a full 3D FFT of the data, and then derive.
  ! 
  ! However, for a derivative in one direction (ipol) only, this is 
  ! not necessary, and a series of 1D FFTs are enough.
  ! This version is used if use_3d_fft=.false.

  use_3d_fft = .true.

  if(use_3d_fft) then

     allocate(aux(0:n1-1, 0:n2-1, 0:n3-1),stat=status)
     call errore('derivative','error allocating array aux',status)
     !
     ! First get the FFT of data
     !
     call fft_3d(n1,n2,n3,data,aux,.true.)
     !
     if(ipol.eq.1) then
        G = 2.d0*pi/L1
        do i1=0,n1-1
           i = i1
           if(i1.gt.n1/2) i = i-n1
           if(i1.eq.n1/2) i = 0
           aux(i1,:,:) = aux(i1,:,:) * cmplx(0.d0,G*i,8)
        enddo
     else if(ipol.eq.2) then
        G = 2.d0*pi/L2
        do i1=0,n2-1
           i = i1
           if(i1.gt.n2/2) i = i-n2
           if(i1.eq.n2/2) i = 0
           aux(:,i1,:) = aux(:,i1,:) * cmplx(0.d0,G*i,8)
        enddo
     else if(ipol.eq.3) then
        G = 2.d0*pi/L3
        do i1=0,n3-1
           i = i1
           if(i1.gt.n3/2) i = i-n3
           if(i1.eq.n3/2) i = 0
           aux(:,:,i1) = aux(:,:,i1) * cmplx(0.d0,G*i,8)
        enddo
     else
        call errore('derivative','wrong value for ipol',1)
     endif
     !
     ! Now go back to real space
     !
     call fft_3d(n1,n2,n3,deriv,aux,.false.)


     deallocate(aux)

  else ! use 1D-FFTs only
     if(ipol.eq.1) then
        G = 2.d0*pi/L1
        allocate(aux_real_1d(0:n1-1),stat=status)
        call errore('derivative','error allocating array aux_real_1d',status)
        allocate(aux_1d(0:n1-1),stat=status)
        call errore('derivative','error allocating array aux_1d',status)
        
        do i3=0,n3-1
           do i2=0,n2-1
              aux_real_1d(0:n1-1) = data(0:n1-1,i2,i3)
              call fft_1d(n1,aux_real_1d,aux_1d,.true.)
              do i1=0,n1-1
                 i = i1
                 if(i1.gt.n1/2) i = i-n1
                 if(i1.eq.n1/2) i = 0
                 aux_1d(i1) = aux_1d(i1) * cmplx(0.d0,G*i,8)
              enddo
              call fft_1d(n1,aux_real_1d,aux_1d,.false.)
              deriv(0:n1-1,i2,i3) = aux_real_1d(0:n1-1)
           enddo
        enddo
     else if(ipol.eq.2) then
        G = 2.d0*pi/L2
        allocate(aux_real_1d(0:n2-1),stat=status)
        call errore('derivative','error allocating array aux_real_1d',status)
        allocate(aux_1d(0:n2-1),stat=status)
        call errore('derivative','error allocating array aux_1d',status)
        
        do i3=0,n3-1
           do i1=0,n1-1
              aux_real_1d(0:n2-1) = data(i1,0:n2-1,i3)
              call fft_1d(n2,aux_real_1d,aux_1d,.true.)
              do i2=0,n2-1
                 i = i2
                 if(i2.gt.n2/2) i = i-n2
                 if(i2.eq.n2/2) i = 0
                 aux_1d(i2) = aux_1d(i2) * cmplx(0.d0,G*i,8)
              enddo
              call fft_1d(n2,aux_real_1d,aux_1d,.false.)
              deriv(i1,0:n2-1,i3) = aux_real_1d(0:n2-1)
           enddo
        enddo
     elseif(ipol.eq.3) then
        G = 2.d0*pi/L3
        allocate(aux_real_1d(0:n3-1),stat=status)
        call errore('derivative','error allocating array aux_real_1d',status)
        allocate(aux_1d(0:n3-1),stat=status)
        call errore('derivative','error allocating array aux_1d',status)
        
        do i2=0,n2-1
           do i1=0,n1-1
              aux_real_1d(0:n3-1) = data(i1,i2,0:n3-1)
              call fft_1d(n3,aux_real_1d,aux_1d,.true.)
              do i3=0,n3-1
                 i = i3
                 if(i3.gt.n3/2) i = i-n3
                 if(i3.eq.n3/2) i = 0
                 aux_1d(i3) = aux_1d(i3) * cmplx(0.d0,G*i,8)
              enddo
              call fft_1d(n3,aux_real_1d,aux_1d,.false.)
              deriv(i1,i2,0:n3-1) = aux_real_1d(0:n3-1)
           enddo
        enddo
     else
        call errore('derivative','wrong value of variable ipol',1)
     endif

     deallocate(aux_1d)
     deallocate(aux_real_1d)

  endif

end subroutine derivative
  
