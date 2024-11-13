      PROGRAM grbwdth_Driver
      IMPLICIT NONE
      INTEGER NE
      PARAMETER (NE=100)
      DOUBLE PRECISION EAR(0:NE),PHOTAR(NE),param(5)
      REAL *8 alpha,width,ep,epv,beta
      INTEGER ifl
      INTEGER i,j
      REAL *8 emin,emax,nmin,nmax,temp,fp,Planck,aa,bb,c

      emin = 4.0e+0 
      emax = 1.0e+7 
      DO i = 0,NE
         EAR(i)=emin*(emax/emin)**(i*1.0/NE)
      ENDDO

      DO i = 1,NE
         PHOTAR(i) = 0.d0
      ENDDO

      alpha = -0.5
      width = 1.16403
      ep = 300.d0
      epv = 100.0
      c = -2.5

      param(1) = alpha
      param(2) = ep
      param(3) = width
      param(4) = epv
      param(5) = c

      Planck    =  6.6260755E-27

      DO i = NE,1,-1
          nmax = EAR(i)*1.6021765e-09/Planck
          nmin = EAR(i-1)*1.6021765e-09/Planck
            
          aa = DLOG(nmin)
          bb = DLOG(nmax)
          CALL bandn(EAR,NE,param,ifl,PHOTAR,beta)

      ENDDO

      fp = 0.0
      DO i =1,NE
           if((EAR(i)*PHOTAR(i)).GT.fp)fp = (EAR(i)**2*PHOTAR(i))
      ENDDO

      print*,"#beta=",beta
      DO i =1,NE
           print *,EAR(i),EAR(i)**2 *PHOTAR(i)/fp
      ENDDO

      END
