!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_3d(n1,n2,n3, data_direct, data_rec, direct_to_reciprocal)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! This subroutine uses fftw to calculate 3-dimensional discrete FFTs.
  ! The data in direct space is assumed to be real-valued
  ! The data in reciprocal space is complex. 
  ! direct_to_reciprocal indicates in which direction the FFT is to be calculated
  ! 
  ! Note that for real data in direct space (like here), we have
  ! F(N-j) = conj(F(j)) where F is the array in reciprocal space.
  ! Here, we do not make use of this property.
  ! Also, we do not use the special (time-saving) routines of FFTW which
  ! allow one to save time and memory for such real-to-complex transforms.
  !
  ! f: array in direct space
  ! F: array in reciprocal space
  ! 
  ! F(k) = \sum_{l=0}^{N-1} exp(- 2 \pi I k*l/N) f(l)
  ! f(l) = 1/N \sum_{k=0}^{N-1} exp(+ 2 \pi I k*l/N) F(k)
  !
  ! See also: http://www.fftw.org/fftw3_doc/Calling-FFTW-from-Modern-Fortran.html#Calling-FFTW-from-Modern-Fortran
  !

  use, intrinsic :: iso_c_binding
  
  implicit none

  include 'fftw3.f03'
  
  integer, intent(in) :: n1,n2,n3
  real(C_DOUBLE), dimension(0:n1-1,0:n2-1,0:n3-1), intent(inout) :: data_direct
  complex(C_DOUBLE_COMPLEX), dimension(0:n1-1,0:n2-1,0:n3-1), intent(inout) :: data_rec
  logical, intent(in) :: direct_to_reciprocal

  type(C_PTR) :: plan  
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), pointer :: fftw_data
  type(C_PTR) :: p
  real(C_DOUBLE) :: fac

  ! Allocate aligned memory
  ! See details here:
  ! http://www.fftw.org/fftw3_doc/Allocating-aligned-memory-in-Fortran.html#Allocating-aligned-memory-in-Fortran

  p = fftw_alloc_complex(int(n1*n2*n3, C_SIZE_T))
  call c_f_pointer(p, fftw_data, [n1,n2,n3])
  !
  ! Now distinguish in which direction the FFT is performed
  !
  if(direct_to_reciprocal) then

     plan = fftw_plan_dft_3d(n3,n2,n1, fftw_data, fftw_data, FFTW_FORWARD,FFTW_ESTIMATE)
     fftw_data(:,:,:) = cmplx(data_direct(:,:,:),0.d0,C_DOUBLE_COMPLEX)
     
     call fftw_execute_dft(plan, fftw_data, fftw_data)

     data_rec(:,:,:) = fftw_data(:,:,:)

     call fftw_destroy_plan(plan)
     call fftw_free(p)

  else
     
     plan = fftw_plan_dft_3d(n3,n2,n1, fftw_data, fftw_data, FFTW_BACKWARD,FFTW_ESTIMATE)
     fftw_data(:,:,:) = data_rec(:,:,:)
     
     call fftw_execute_dft(plan, fftw_data, fftw_data)

     data_direct(:,:,:) = real(fftw_data(:,:,:), C_DOUBLE)
     
     fac = 1.d0/(n1*n2*n3)

     data_direct(:,:,:) = data_direct(:,:,:) * fac

     call fftw_destroy_plan(plan)
     call fftw_free(p)

  endif

end subroutine fft_3d
  



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_1d(n1, data_direct, data_rec, direct_to_reciprocal)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! This subroutine uses fftw to calculate 1-dimensional discrete FFTs.
  ! The data in direct space is assumed to be real-valued
  ! The data in reciprocal space is complex. 
  ! direct_to_reciprocal indicates in which direction the FFT is to be calculated
  ! 
  ! Note that for real data in direct space (like here), we have
  ! F(N-j) = conj(F(j)) where F is the array in reciprocal space.
  ! Here, we do not make use of this property.
  ! Also, we do not use the special (time-saving) routines of FFTW which
  ! allow one to save time and memory for such real-to-complex transforms.
  !
  ! f: array in direct space
  ! F: array in reciprocal space
  ! 
  ! F(k) = \sum_{l=0}^{N-1} exp(- 2 \pi I k*l/N) f(l)
  ! f(l) = 1/N \sum_{k=0}^{N-1} exp(+ 2 \pi I k*l/N) F(k)
  !
  ! See also: http://www.fftw.org/fftw3_doc/Calling-FFTW-from-Modern-Fortran.html#Calling-FFTW-from-Modern-Fortran
  !

  use, intrinsic :: iso_c_binding
  
  implicit none

  include 'fftw3.f03'
  
  integer, intent(in) :: n1
  real(C_DOUBLE), dimension(0:n1-1), intent(inout) :: data_direct
  complex(C_DOUBLE_COMPLEX), dimension(0:n1-1), intent(inout) :: data_rec
  logical, intent(in) :: direct_to_reciprocal

  type(C_PTR) :: plan  
  complex(C_DOUBLE_COMPLEX), dimension(:), pointer :: fftw_data
  type(C_PTR) :: p
  real(C_DOUBLE) :: fac

  ! Allocate aligned memory
  ! See details here:
  ! http://www.fftw.org/fftw3_doc/Allocating-aligned-memory-in-Fortran.html#Allocating-aligned-memory-in-Fortran

  p = fftw_alloc_complex(int(n1, C_SIZE_T))
  call c_f_pointer(p, fftw_data, [n1])
  !
  ! Now distinguish in which direction the FFT is performed
  !
  if(direct_to_reciprocal) then

     plan = fftw_plan_dft_1d(n1, fftw_data, fftw_data, FFTW_FORWARD,FFTW_ESTIMATE)
     fftw_data(:) = cmplx(data_direct(:),0.d0,C_DOUBLE_COMPLEX)
     
     call fftw_execute_dft(plan, fftw_data, fftw_data)

     data_rec(:) = fftw_data(:)
  else
     
     plan = fftw_plan_dft_1d(n1, fftw_data, fftw_data, FFTW_BACKWARD,FFTW_ESTIMATE)
     fftw_data(:) = data_rec(:)
     
     call fftw_execute_dft(plan, fftw_data, fftw_data)

     data_direct(:) = real(fftw_data(:), C_DOUBLE)
     
     fac = 1.d0/real(n1,8)

     data_direct(:) = data_direct(:) * fac

  endif

  call fftw_destroy_plan(plan)
  call fftw_free(p)


end subroutine fft_1d
  
