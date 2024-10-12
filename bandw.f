      subroutine band(ear, ne, param, ifl, photar, photer)
C
      IMPLICIT NONE
      INTEGER          ne, ifl, ie
      DOUBLE PRECISION ear(0:ne), param(3), photar(ne), photer(ne)
      DOUBLE PRECISION wd,e_p,bisect,lambert,bet,temp
      DOUBLE PRECISION nprime, alfa, beta, amb,e_piv,fac
      Double precision tem, temi, e_b,betap1,b1,b2,w1,w2 
      Double precision elo, vlo, ehi, vhi, foo, output
      PARAMETER        (e_piv=1.d+2)
      COMMON           /grbpar/wd,alfa,temp
      external         bisect,lambert,bet
C
C-->  Approximation to the Gamma-Ray Burst model by D. Band, et al.,
C     1993 (ApJ 413, 281) for XSPEC programs:
C     N(E) = n  * (E/E_p)**alfa * exp(-E/tem), if E < (alfa-beta)tem
C          = n' * (E/E_p)**beta              , if E > (alfa-beta)tem
C
C     The parameters in order are n, alfa, beta, and tem.
C     Suggested default values and increments are (note that
C     the normalization "n" is handled outside this function):
C         n    = 0.01        Increment = 1e-4
C         alfa = -1.0        Increment = 1e-2
C         beta = -2.0        Increment = 1e-2
C         tem  = 300.        Increment = 10.
C
C-->  Revisions:
C     Written and revised for BSAS by D. Band and S. Bansal (1/27/92; 
C     7/06/92; 3/02/93);  Adapted for TDAS by S. Bansal (4/25/94); 
C     Modified for XSPEC by H. Seifert (7/23/96);
C

c suppress a warning message from the compiler
      ie = ifl

c this model does not calculate errors
      DO ie = 1, ne
         photer(ie) = 0.0
      ENDDO

C-->  Initialize variables...
      alfa = param(1)
      wd = param(2)
      e_p  = param(3)
      
      fac = 2.0+alfa
      tem = e_p/fac
      
      b1 = -2.2d0
      b2 = -10.0d0
       
      w1=0d0
      w2=-1d0
      
      temp = bisect(w1,w2,lambert)
      beta = bisect(b1,b2,bet)
C
C-->  Define some secondary variables...
      temi = 1.d0/tem
C
      amb = alfa-beta
      if (amb .lt. temi) amb = temi
      e_b = amb*tem
C
      nprime = ((e_b/e_piv)**amb)*exp(-amb)
      betap1 = beta+1.d0

C
C-->  Get first data point and consider the cases
C     when E is >/< (alfa-beta)tem; finally, also   
C     consider also the special case of beta=-1.
      elo = ear(0)
      vlo = 0.0d0
      if (elo .le. e_b) then
         vlo = ((elo/e_piv)**alfa)*exp(-elo*temi)
      else
         if (beta .eq. -1.d0) then
            vlo = nprime*e_piv*LOG(elo)
         else
            vlo = nprime*e_piv*((elo/e_piv)**betap1)/betap1
         endif
      endif

C
C-->  Loop over all energy bins, considering again
C     the cases when E is >/< (alfa-beta)tem and, again,
C     also consider the special case of beta=-1.
C     Note that if both elo and ehi are >(alfa-beta)tem,
C     we can do the exact integral; otherwise, we do a
C     simple two-point approximation.
      do ie = 1, ne
         ehi = ear(ie)
         if (ehi .le. e_b) then
            vhi = ((ehi/e_piv)**alfa)*exp(-ehi*temi)
         else
            if (elo .gt. e_b) then
               if (beta .eq. -1.d0) then
                  vhi = nprime*e_piv*LOG(ehi)
               else
                  vhi = nprime*e_piv*((ehi/e_piv)**betap1)/betap1
               endif
            else
               vhi = nprime*(ehi/e_piv)**beta
               foo = nprime*e_piv*((ehi/e_piv)**betap1)/betap1
            endif
         endif
C-->     Fill output array... 
         if (elo .gt. e_b) then
            output = vhi-vlo  
         else
            output = (vhi+vlo)*(ehi-elo)/2.d0
         endif

         IF ( output .LT. 2.0e-38 ) THEN
            photar(ie) = 0.0
         ELSE IF ( output .GT. 2.0e38 ) THEN
            photar(ie) = 2.0e38
         ELSE
            photar(ie) = SNGL(output)
         ENDIF

C-->     Set up for next loop...
         if ((ehi .gt. e_b) .and. (elo .le. e_b)) then
            vlo = foo
         else
            vlo = vhi
         endif
         elo = ehi
      enddo
C
      return
      end
      
!=================================================================
       real*8 function bisect(a,b,f)
       implicit none
       
       double precision :: a,b,f,c,x1,x2
       integer:: i,n
       double precision wd,alp,temp,ep,beta
        COMMON /grbpar/wd,alp,temp
        external f

       n = 1000
       do i = 1,n
	  x1 = a
	  do while( x1 >= b)
              x2 = x1 - 0.1d0
              do while (x2 >= b)
                  if (f(x2)== 0d0) then
                      bisect = x2
                      exit
                  end if
                  if (f(x1) * f(x2) < 0.0d0) exit
                  x2 = x2 - 0.1d0
              end do
              if (f(x1) * f(x2) < 0.0d0) exit
              if (f(x1)== 0d0) then
                      bisect = x1
                      exit
              end if
              if (f(x2)== 0d0) then
                      bisect = x2
                      exit
              end if
              x1 = x1 - 0.1d0
          end do
         
	 if (f(x1) * f(x2) <= 0.0d0) exit
		!print*,'Error: Could not find valid initial values'
		wd = wd + ((-1d0)**n)*rand()/10d0
		if (i==n) then
		print*,'Error: Could not find valid initial values'
		stop
		end if
	 !end if
	 end do


        do i = 1,n
          c = (x1 + x2)/2.0d0
          if (f(x1)*f(c)<0.0d0) then
            x2 = c
          else if (f(x2)*f(c)<0.0d0) then
            x1 = c
          end if
          if (f(c) == 0.0d0 .or. abs((x1-x2)/x1)<1.0d-12) then
            bisect = c
            exit
          end if
          !if (i==n) then
            !print*,'Bisection did not converge ',a,b
          !end if
       end do

       

      end 

!=======================================================
       DOUBLE PRECISION FUNCTION Lambert(w)
       IMPLICIT NONE
       DOUBLE PRECISION ene,alfa,ep,temp,wd,w,x,beta
       COMMON /grbpar/wd,alfa,temp

       x = -1.0d0/(exp(1.0d0)*2.0d0**(1.0d0/(alfa+2.0d0)))

       Lambert = w*exp(w)-x
      

       END
!=======================================================
       double precision function bet(b)
       implicit none

       double precision:: wd,alfa,a,temp,d,b,ep,k,beta
       COMMON /grbpar/wd,alfa,temp

       a = alfa + 2.0d0
       k = dexp(1.0d0)*2.0d0**(1.0d0/(b + 2.0d0))
       d = k*(temp)*((alfa - b)**((alfa - b)/(b + 2.0d0)))

       bet = dlog10(-(a)**((alfa - b)/(b + 2.0d0))/d) - wd 
       
       end
