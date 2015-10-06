F90           = gfortran
LIBS          = -lfftw3		
F90FLAGS      = -I/usr/include -I${FFTW3_INC} -O3
LDFLAGS       =
POBJS         = p_errore.o p_plot_data_1d.o  p_plot_data_2d.o p_fft_wrapper.o p_derivative.o 

all : diffusion.x 

%.o : %.f90
	$(F90) $(F90FLAGS) -c $<

diffusion.x : $(POBJS) p_diffusion.o 
	$(F90) $(LDFLAGS) -o $@ p_diffusion.o $(POBJS) $(LIBS) 

clean :
	- /bin/rm -f *.x *.o *.mod *~ *.a *.F90 concentration* diffusivity* 1d_conc*

