       SUBROUTINE band(EAR,NE,param,ifl,PHOTAR,width)
       IMPLICIT NONE

       INTEGER NE,IFL,ie,j
       DOUBLE PRECISION EAR(0:NE),PHOTAR(NE),width
       DOUBLE PRECISION PARAM(4),ene

       DOUBLE PRECISION temp,Eh1,Eh2,Eh3,alp,beta,e_p,wd,e_b,epiv
       DOUBLE PRECISION nprime,fac,alb,temi,e0,aa,bb,const
       DOUBLE PRECISION  betap1,elo,vlo,ehi,vhi,foo,pflux,b1
       DOUBLE PRECISION Bisection_Solve,DiffFlux 
       EXTERNAL Bisection_Solve,DiffFlux
       COMMON/betapar/Eh1,wd,epiv
       
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
       

      if (e0 .lt. 1.d0) e0 = 1.d0
      temi = 1.0d0/e0
       
       alb = alp-beta
       if (alb .lt. temi) alb = temi
       e_b = e0*alb

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

      ! DO j = 1,NE
      !    ene = EAR(j)
      !    IF(ene.LE.e_b) THEN
      !       PHOTAR(j) = ((ene/epiv)**alp)*DEXP(-ene/e0)
      !    ELSE
      !       PHOTAR(j) = nprime*(ene/epiv)**beta
      !    ENDIF
      ! ENDDO

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

