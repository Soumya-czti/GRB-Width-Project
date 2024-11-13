       SUBROUTINE width_calc(param,width)
       IMPLICIT NONE

       INTEGER j
       DOUBLE PRECISION PARAM(4)

       DOUBLE PRECISION temp,Eh1,Eh3,alp,beta,e_p,width,epiv
       DOUBLE PRECISION e0,fac,aa,bb,pflux
       DOUBLE PRECISION Bisection_Solve,DiffFlux 
       EXTERNAL Bisection_Solve,DiffFlux
       
       alp = param(1) 
       e_p = param(2)  !in keV
       beta = param(3)  
       epiv = param(4) !in keV

       fac = 2.0+alp
       e0 = e_p/fac
       pflux = (e_p/epiv)**(alp+2.d0)*DEXP(-(2.d0+alp))
       
       aa = 1.0
       bb = e_p
       Eh1 = Bisection_Solve(fac,e_p,e0,aa,bb,DiffFlux)

       Eh3 = e_p**fac/DEXP(2.d0+beta)/2.d0
       Eh3 = Eh3/((alp-beta)*e0)**(alp-beta)
       Eh3 = Eh3**(1.d0/(beta+2.d0))
       width = DLOG10(Eh3/Eh1)
       RETURN
       END

!=======================================================

       DOUBLE PRECISION FUNCTION Bisection_Solve(fac,e_p,e0,aa,bb,func)
       IMPLICIT NONE
       DOUBLE PRECISION emin,emax,val,fac,e_p,e0  
       DOUBLE PRECISION func,fp1,fp2,fmid
       REAL *8 aa,bb,temp1,temp2,temp
       INTEGER i
       EXTERNAL func

       temp1 = func(aa,fac,e_p,e0)
       temp2 = func(bb,fac,e_p,e0)
       DO i = 1,100
          fmid = (aa+bb)/2.0
          temp = func(fmid,fac,e_p,e0)
          IF(temp*temp1.LT.0.d0)THEN
              bb = fmid
          ELSE
              aa = fmid
          ENDIF
          IF(DABS(aa-bb).LT.1.0d-2)EXIT
!         print *, "***", aa,bb
       ENDDO
       Bisection_Solve = fmid

       RETURN
       END

!=======================================================
       DOUBLE PRECISION FUNCTION DiffFlux(ene,fac,e_p,e0)
       IMPLICIT NONE
       DOUBLE PRECISION ene,fac,e_p,e0,k

       k = (DEXP(-fac)*(e_p**fac))/2.0
       DiffFlux = ((ene**fac)*DEXP(-ene/e0))-k
       
       RETURN
       END
!=======================================================

