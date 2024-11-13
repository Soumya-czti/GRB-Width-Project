       SUBROUTINE beta_calc(param,beta)
       IMPLICIT NONE

       DOUBLE PRECISION PARAM(5)

       DOUBLE PRECISION temp,Eh1,Eh2,Eh3,alp,beta,e_p,wd,e_b,epiv
       DOUBLE PRECISION aa,bb,const,pflux,fac,e0
       DOUBLE PRECISION Bisection_Solve,Betafunc,DiffFlux 
       EXTERNAL Bisection_Solve,Betafunc,DiffFlux
       COMMON/betapar/Eh1,wd,epiv
       
       alp = param(1) 
       e_p = param(2)  !in keV
       wd = param(3)  !in log (wd=log(E2/E1))
       epiv = param(4) !in keV
       const = param(5) !beta constant

       fac = 2.0+alp
       e0 = e_p/fac
       pflux = (e_p/epiv)**(alp+2.d0)*DEXP(-(2.d0+alp))
       
       aa = 1.0
       bb = e_p
       Eh1 = Bisection_Solve(fac,e_p,e0,aa,bb,DiffFlux)

       aa = const-1.0
       bb = const+1.0
       beta = Bisection_Solve(fac,pflux,e0,aa,bb,Betafunc)

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

       DOUBLE PRECISION FUNCTION Betafunc(beta,fac,pflux,e0)
       IMPLICIT NONE
       DOUBLE PRECISION pflux,e0,Eh3,beta,alp,fac,Eh1,wd,temp,epiv
       COMMON/betapar/Eh1,wd,epiv

       alp = fac - 2.d0
       Eh3 = wd+DLOG10(Eh1)
       Eh3 = 10.d0**Eh3

       temp = ((alp-beta)*(e0/epiv))**(alp-beta)*DEXP(beta-alp)
       temp = temp*(Eh3/epiv)**(beta+2.d0)
       Betafunc = temp-pflux/2.d0
       
       RETURN
       END
!=======================================================
