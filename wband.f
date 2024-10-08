       SUBROUTINE band(EAR,NE,param,ifl,PHOTAR,PHOTER)
       IMPLICIT NONE

       INTEGER NE,IFL,ie
       DOUBLE PRECISION EAR(0:NE),PHOTAR(NE),PHOTER(NE)
       DOUBLE PRECISION PARAM(4)

       DOUBLE PRECISION temp,alp,beta,e_p,wd,e_b,epiv
       DOUBLE PRECISION nprime,tem,fac,alb,temi
       DOUBLE PRECISION  betap1,elo,vlo,ehi,vhi,foo
       DOUBLE PRECISION bisect,lambert,bet 
       DOUBLE PRECISION w1,w2,b1,b2
       DOUBLE PRECISION alp1,beta1,wd1
       save alp1, beta1, wd1
       EXTERNAL bisect,lambert,bet
       
       alp1 = -1.0
       wd1 = 1.0
       
       alp = param(1) 
       e_p = param(2)  !in keV
       wd = param(3)  !in log (wd=log(E2/E1))
       epiv = param(4) !in keV

       fac = 2.0+alp
       tem = e_p/fac
       if (tem .lt. 1.d0) tem = 1.d0
       temi = 1.0d0/tem
       
       if (alp == alp1 .and. wd == wd1) then
           beta = beta1
       else
 
           w1=0
           w2=-1
       
           temp = bisect(w1,w2,lambert)
       !t = temp
       
           b1 = -2.2
           b2 = -5.0
       
           beta = bisect(b1,b2,bet)
       
       end if
       
       alp1 = alp
       wd1 = wd
       beta1 = beta
       
       alb = alp-beta
       if (alb .lt. temi) alb = temi
       e_b = tem*alb

       nprime = ((e_b/epiv)**alb)*DEXP(-alb)
       betap1 = beta+1.d0
      
       elo = EAR(0)
       if (elo .le. e_b) then
         vlo = ((elo/epiv)**alp)*exp(-elo*temi)
       else
         if (beta .eq. -1.d0) then
            vlo = nprime*epiv*LOG(elo)
         else
            vlo = nprime*epiv*((elo/epiv)**betap1)/betap1
         endif
       endif
       
       DO ie = 1, ne
         ehi = EAR(ie)
         if (ehi .le. e_b) then
            vhi = ((ehi/epiv)**alp)*exp(-ehi*temi)
         else
            if (elo .gt. e_b) then
               if (beta .eq. -1.d0) then
                  vhi = nprime*epiv*LOG(ehi)
               else
                  vhi = nprime*epiv*((ehi/epiv)**betap1)/betap1
               endif
            else
               vhi = nprime*(ehi/epiv)**beta
               foo = nprime*epiv*((ehi/epiv)**betap1)/betap1
            endif
         endif

         if (elo .gt. e_b) then
            PHOTAR(ie) = vhi-vlo  
         else
            PHOTAR(ie) = (vhi+vlo)*(ehi-elo)/2.d0
         endif
         
         if ((ehi .gt. e_b) .and. (elo .le. e_b)) then
            vlo = foo
         else
            vlo = vhi
         endif
         elo = ehi
      enddo

       !DO j = 1,NE
       !   ene = EAR(j)
       !   IF(ene.LE.e_b) THEN
       !      PHOTAR(j) = ((ene/epiv)**alp)*DEXP(-ene/tem)
       !   ELSE
       !      PHOTAR(j) = nprime*(ene/epiv)**beta
       !   ENDIF
       !ENDDO

       RETURN
       END

!=======================================================

       real*8 function bisect(a,b,f)
       implicit none
       
       double precision width,alp,temp,ep,beta
       COMMON /grbpar/width,alp,temp,ep,beta
       real*8 :: a,b,f,c
       integer:: i,n

       external f

       n = 100

       if (f(a)*f(b)>= 0) then
          !print*,'Warning: Invalid initial guesses.Searching for valid guesses...'
	  a = -2.2d0
	  do while( a >= -10.0d0)
			b = a - 0.1d0
              do while (b >= -10.0d0)
                  if (f(a) * f(b) < 0.0d0) exit
                  b = b - 0.1d0
              end do
              if (f(a) * f(b) < 0.0d0) exit
              a = a - 0.1d0
         end do
         !if (f(a) * f(b) >= 0.0d0) then
          !    width = 0.5
          !    do while(width <= 2.0d0)
          !  	  if (f(a)*f(b)>= 0.0d0) then
	!			a = -2.2d0
	!			do while( a >= -10.0d0)
	!				b = -10.0d0
	!				do while (b < a)
	!					if (f(a) * f(b) < 0.0d0) exit
	 !					     b = b + 0.1d0
	!			        end do
	!					     if (f(a) * f(b) < 0.0d0) exit
	!					     a = a - 0.1d0
	!			end do
	!	endif
	!	width = width + 1.0d-3
	!	end do
	!end if
	if (f(a) * f(b) >= 0.0d0) then
		print*,'Error: Could not find valid initial values'
		stop
	end if
       end if


       do i = 1,n
          c = (a + b)/2.0d0
          if (f(c)*f(b)<0.0d0) then
            a = c
          else if (f(a)*f(c)<0.0d0) then
            b = c
          end if
          if (f(c) == 0.0d0 .or. abs(a-b)<1.0d-6) then
            exit
          end if
          if (i==n) then
            print*,'Bisection did not converge ',a,b
          end if
       end do

       bisect = c

       end 

!=======================================================
       DOUBLE PRECISION FUNCTION Lambert(w)
       IMPLICIT NONE
       DOUBLE PRECISION ene,alp,ep,temp,width,w,x,beta
       COMMON /grbpar/width,alp,temp,ep,beta

       x = -1.0d0/(exp(1.0d0)*2.0d0**(1/(alp+2.0d0)))

       Lambert = w*exp(w)-x
      


       RETURN
       END
!=======================================================
       double precision function bet(b)
       implicit none

       double precision:: width,alp,a,temp,d,b,ep,k,beta
       COMMON /grbpar/width,alp,temp,ep,beta

       a = alp + 2.0d0
       k = dexp(1.0d0)*2.0d0**(1.0d0/(b + 2.0d0))
       d = k*(temp)*((alp - b)**((alp - b)/(b + 2.0d0)))

       bet = dlog10((-a)**((alp - b)/(b + 2.0d0))/d) - width

       end function
