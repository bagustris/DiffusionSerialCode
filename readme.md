#Diffusion Serial Code, run job with the right modules

<h3> Step by step: </h3>
1. Login to master/frontend on cluster as root
2. Download the latest fftw3
3. Make directory on /opt/fftw3
3. Install fftw3 as single precision with:
  - ./configure --prefix=/opt//fftw3 --disable-shared \
--enable-static --enable-single --enable-fortran
  - make
  - make install
4. Install fftw3 with double precision with:
make clean
  - ./configure --prefix=/opt/soft/libs/fftw3 --disable-shared \
--enable-static --enable-fortran
  - make
  - make install
5. Set the prefix in file fftw according to installation path (/opt/fftw3)
6. mv fftw to modulefiles (/usr/share/Module/modulefiles)
7. Compile job
  - make all
8. Load module >> module load fftw3
9. Make bash script for the jobb >> diffusion.sh
10. Submit job to the cluster >> qsub diffusion.sh
11. Check if there is any error:
  - qstat -f
  - more diffusion.sh.oXX
