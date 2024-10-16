      program check
      implicit none
      
      real*8:: width,alpha,e_p,bta
      integer:: i
      
      open(11,file='wid_re.dat')
      open(10,file='beta_re.dat')
      
      do i = 1,15
      read(11,*) alpha,bta,e_p,width
      
      call band(alpha,width)
      write(10,*) 
      end do
      
      end program
      
      subroutine band(alp,width)

      IMPLICIT NONE
      INTEGER          ne, ifl, ie
      !DOUBLE PRECISION ear(0:ne), param(3), photar(ne), photer(ne)
      DOUBLE PRECISION wd,e_p,bisect,lambert,bet,temp,bta,bisect1
      DOUBLE PRECISION alfa, beta,alp,width,wd1,beta1,beta2
      Double precision tem, temi, e_b,betap1,b1,b2,w1,w2 
      Double precision elo, vlo, ehi, vhi, foo, output
      
      COMMON           /grbpar/wd,alfa,temp,bta
      external         bisect,lambert,bet,bisect1

      alfa = alp
      wd = width
      !e_p  = param(3)
      
      !fac = 2.0+alfa
      !tem = e_p/fac
      
      b1 = -2.2d0
      b2 = -50.0d0
       
      w1=0d0
      w2=-1d0
      
      temp = bisect1(w1,w2,lambert)
      print*,''
      beta = bisect(b1,b2,bet)
      !print*,beta
      
      !print*,''
      end
!=================================================================
real*8 function bisect(a,b,f)
implicit none
       
double precision :: a,b,f,c,x1,x2,a1,a2,h,x
integer:: i,n,j,k
double precision wd,alp,temp,ep,beta,wd1
COMMON /grbpar/wd,alp,temp
external f

n = 1000   
print*,'a',a
x = a
x1 = x
print*,x,x1,'i'
h = 0.01d0

do while (x1 .ge. b .and. x1 .le. a)
	a1 = x
	a2 = x1
	!print*,x1
	if (f(x1) == 0.0d0) then
		print*, "Root:",x1
		bisect = x1
		write(10,'(F10.3, 2F10.3)',advance='no') bisect
		
	elseif (f(a1)*f(a2) < 0.0d0) then
		do j = 1,10000
			!print*,a1,a2
			c = (a1+a2)/2.0d0
			!print*,c
			if (f(c)*f(a2) < 0.0d0) then
				a1 = c
			elseif (f(c)*f(a1) < 0.0d0) then
				a2 = c
			end if
			if (f(c) == 0.0d0 .or. abs((a2-a1)/a1) <= 1d-12) then
				print*,'Converged',c,f(c)
				bisect = c
				write(10,'(F10.3, 2F10.3)',advance='no') bisect
				exit
			end if
		end do
		!x = x1
	end if	
		x = x1
		x1 = x1 - h
end do

      end 

!=======================================================

real*8 function bisect1(a,b,f)
implicit none
       
double precision :: a,b,f,c,x1,x2,a1,a2,h,x
integer:: i,n,j,k
double precision wd,alp,temp,ep,beta,wd1
COMMON /grbpar/wd,alp,temp
external f

n = 1000   
!print*,'a',a
x = a
x1 = x
!print*,x,x1,'i'
h = 0.01d0

do while (x1 .ge. b .and. x1 .le. a)
	a1 = x
	a2 = x1
	!print*,x1
	if (f(x1) == 0.0d0) then
		!print*, "Root:",x1
		bisect1 = x1
		!write(10,'(F10.3, 2F10.3)',advance='no') bisect
		
	elseif (f(a1)*f(a2) < 0.0d0) then
		do j = 1,10000
			!print*,a1,a2
			c = (a1+a2)/2.0d0
			!print*,c
			if (f(c)*f(a2) < 0.0d0) then
				a1 = c
			elseif (f(c)*f(a1) < 0.0d0) then
				a2 = c
			end if
			if (f(c) == 0.0d0 .or. abs((a2-a1)/a1) <= 1d-12) then
				!print*,'Converged',c,f(c)
				bisect1 = c
				!write(10,'(F10.3, 2F10.3)',advance='no') bisect
				exit
			end if
		end do
		!x = x1
	end if	
		x = x1
		x1 = x1 - h
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
