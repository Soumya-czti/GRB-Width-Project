      PROGRAM grbwdth_Driver
      IMPLICIT NONE
      INTEGER NE
      PARAMETER (NE=100)
      DOUBLE PRECISION EAR(0:NE),PHOTAR(NE),param(4), PHOTER(NE)
      REAL *8 alpha,width,ep,epv
      INTEGER ifl
      INTEGER i,j
      REAL *8 emin,emax,nmin,nmax,temp,fp,Planck,aa,bb,beta,e0

      emin = 4.0e+0 
      emax = 1.0e+7 
      DO i = 0,NE
         EAR(i)=emin*(emax/emin)**(i*1.0/NE)
      ENDDO

      DO i = 1,NE
         PHOTAR(i) = 0.d0
      ENDDO

      alpha = 0.83
      beta = -3.286
      e0 =  55.407
      ep = e0*(2+alpha)
      epv = 100.0

      param(1) = alpha
      param(2) = ep
      param(3) = beta
      param(4) = epv

      CALL width_calc(param,width)

      print *,"ep=",ep,"beta=",beta,"width calculated=",width

      END
