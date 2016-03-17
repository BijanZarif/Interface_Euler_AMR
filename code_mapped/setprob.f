      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /param/ gammagas, gammaplas, gammawat
      common /param/ pinfgas,pinfplas,pinfwat
      common /param/ omegas,omeplas,omewat
      common /param/ rhog,rhop,rhow
      
      common /mappedparam/rsqr,rout,rinn
      
c     # Read parameter from python setrun.py and loads them into fortran
c     # Set EOS and mapped grid parameters and
c     # passed to the required subroutines in common blocks
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)

      read(7,*) gammagas
      read(7,*) gammaplas
      read(7,*) gammawat
      read(7,*) pinfgas
      read(7,*) pinfplas
      read(7,*) pinfwat
      read(7,*) rhog
      read(7,*) rhop
      read(7,*) rhow
      read(7,*) omegas
      read(7,*) omeplas
      read(7,*) omewat
      read(7,*) rsqr
      read(7,*) rout
      read(7,*) rinn
      

      return
      end

