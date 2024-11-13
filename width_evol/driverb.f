      PROGRAM grbwdth_Driver
      IMPLICIT NONE
      INTEGER NE
      PARAMETER (NE=100)
      DOUBLE PRECISION EAR(0:NE),PHOTAR(NE),param(5), PHOTER(NE)
      REAL *8 alpha,width,ep,epv
      INTEGER ifl
      INTEGER i,j
      REAL *8 emin,emax,nmin,nmax,temp,fp,Planck,aa,bb,beta,c

      emin = 4.0e+0 
      emax = 1.0e+7 
      DO i = 0,NE
         EAR(i)=emin*(emax/emin)**(i*1.0/NE)
      ENDDO

      DO i = 1,NE
         PHOTAR(i) = 0.d0
      ENDDO

      alpha = 0.815726
      width = 0.68510
      ep = 156.075
      epv = 100.0
      c = -3.15

      param(1) = alpha
      param(2) = ep
      param(3) = width
      param(4) = epv
      param(5) = c

      CALL beta_calc(param,beta)

      print *,"beta calculated =",beta,"width=",width

      END
