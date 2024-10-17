       SUBROUTINE band(EAR,alp,beta,ebr,k,PHOTAR)
       IMPLICIT NONE
      
       DOUBLE PRECISION EAR,PHOTAR

Cf2py  intent(in) EAR
Cf2py  intent(in) alp     
Cf2py  intent(in) beta
Cf2py  intent(in) ebr
Cf2py  intent(in) k
Cf2py  intent(out) PHOTAR   

       DOUBLE PRECISION beta,alp,ebr,fac,k,alp
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
       
       
