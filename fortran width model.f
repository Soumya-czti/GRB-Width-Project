       SUBROUTINE band(EAR,alfa,width,ebr,k,PHOTAR)
       IMPLICIT NONE
      
       DOUBLE PRECISION EAR,PHOTAR,bisect,lambert,temp

Cf2py  intent(in) EAR
Cf2py  intent(in) alp     
Cf2py  intent(in) wd
Cf2py  intent(in) ebr
Cf2py  intent(in) k
Cf2py  intent(out) PHOTAR   

       DOUBLE PRECISION beta,alp,ebr,fac,k,wd,alfa,width
       DOUBLE PRECISION b1,b2,w1,w2
       DOUBLE PRECISION bet
       COMMON           /grbpar/wd,alp,temp
       external         bisect,lambert,bet
       alp = alfa
       wd = width
       b1 = -2.2d0
       b2 = -10.0d0
       
       w1=0d0
       w2=-1d0
      
       temp = bisect(w1,w2,lambert)
       beta = bisect(b1,b2,bet)
       fac = ebr*(alp-beta)
         
       IF(EAR.LE.fac) THEN
             PHOTAR = k*((EAR/100.0d0)**alp)
             PHOTAR = PHOTAR*DEXP(-EAR/ebr)
       ELSE    
             PHOTAR = k*(((alp-beta)*(ebr/100.0d0))**(alp-beta))
             PHOTAR = PHOTAR*DEXP(-(alp-beta))
             PHOTAR = PHOTAR*(EAR/100.0d0)**beta
     
       ENDIF

       END  
       
       real*8 function bisect(a,b,f)
       implicit none
       
       double precision :: a,b,f,c,x1,x2
       integer:: i,n
       double precision wd,alp,temp,ep,beta
        COMMON /grbpar/wd,alp,temp
!        external f

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
       DOUBLE PRECISION ene,alp,ep,temp,wd,w,x,beta
       COMMON /grbpar/wd,alp,temp

       x = -1.0d0/(exp(1.0d0)*2.0d0**(1.0d0/(alp+2.0d0)))

       Lambert = w*exp(w)-x
      

       END
!=======================================================
       double precision function bet(b)
       implicit none

       double precision:: wd,alp,a,temp,d,b,ep,k,beta
       COMMON /grbpar/wd,alp,temp

       a = alp + 2.0d0
       k = dexp(1.0d0)*2.0d0**(1.0d0/(b + 2.0d0))
       d = k*(temp)*((alp - b)**((alp - b)/(b + 2.0d0)))

       bet = dlog10(-(a)**((alp - b)/(b + 2.0d0))/d) - wd 
       
       end
