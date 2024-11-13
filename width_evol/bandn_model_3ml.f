       SUBROUTINE bandn(EAR,alp,wd,e_p,const,k,PHOTAR)
       IMPLICIT NONE

       DOUBLE PRECISION EAR,PHOTAR

Cf2py  intent(in) EAR
Cf2py  intent(in) alp     
Cf2py  intent(in) wd
Cf2py  intent(in) e_p
Cf2py  intent(in) const
Cf2py  intent(in) k
Cf2py  intent(out) PHOTAR   

       DOUBLE PRECISION Eh1,alp,beta,e_p,wd,e_b,epiv,width
       DOUBLE PRECISION nprime,fac,alb,e0,aa,bb,const,k
       DOUBLE PRECISION  pflux
       DOUBLE PRECISION Bisection_Solve,Betafunc,DiffFlux 
       EXTERNAL Bisection_Solve,Betafunc,DiffFlux
       COMMON/betapar/Eh1,width,epiv
       
       epiv = 100.0  !in keV
       width = wd

       fac = 2.0+alp
       e0 = e_p/fac
       pflux = (e_p/epiv)**(alp+2.d0)*DEXP(-(2.d0+alp))
       
       aa = 1.0
       bb = e_p
       Eh1 = Bisection_Solve(fac,e_p,e0,aa,bb,DiffFlux)

       aa = const-1.0
       bb = const+1.0
       beta = Bisection_Solve(fac,pflux,e0,aa,bb,Betafunc)
       
       alb = alp-beta
       e_b = e0*alb

       nprime = ((e_b/epiv)**alb)*DEXP(-alb)
      
       IF(EAR.LE.e_b) THEN
          PHOTAR = k*((EAR/epiv)**alp)*DEXP(-EAR/e0)
       ELSE
          PHOTAR = k*nprime*(EAR/epiv)**beta
       ENDIF

       RETURN
       END

!=======================================================

       DOUBLE PRECISION FUNCTION Bisection_Solve(fac,e_p,e0,aa,bb,func)
       IMPLICIT NONE
       DOUBLE PRECISION fac,e_p,e0  
       DOUBLE PRECISION func,fmid
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
       DOUBLE PRECISION pflux,e0,Eh3,beta,alp,fac,Eh1,width,temp,epiv
       COMMON/betapar/Eh1,width,epiv

       alp = fac - 2.d0
       Eh3 = width+DLOG10(Eh1)
       Eh3 = 10.d0**Eh3

       temp = ((alp-beta)*(e0/epiv))**(alp-beta)*DEXP(beta-alp)
       temp = temp*(Eh3/epiv)**(beta+2.d0)
       !print *,"------>",beta,temp
       Betafunc = temp-pflux/2.d0
       
       RETURN
       END
!=======================================================
