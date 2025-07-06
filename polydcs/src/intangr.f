      SUBROUTINE ANGINT(l1,l2,l,mm1,mm2,m,ip1,ip2,ip,result)
      IMPLICIT REAL*4 (A-H, O-Z)
 
      rt = 0.35355339E0

      ip12 = ip1*ip2
      IF ( ip.NE.ip12 ) THEN
        result = 0.0
      ELSE
        xl = l + l + 1
        lmaxa = l1 + l2 + l
        lz = (-1)**lmaxa
        IF ( lz.LE.0 ) THEN
          result = 0.0
        ELSE IF ( l-l1.LE.l2 ) THEN
          lmin = l1 - l2
          IF ( lmin.LT.0 ) lmin = -lmin
          IF ( l.LT.lmin ) THEN
            result = 0.0
          ELSE
c
c     To perform angular integration
c
            mmin = mm1 - mm2
            IF ( mmin.LT.0 ) mmin = -mmin
            IF ( m.NE.mmin ) THEN
              mmax = mm1 + mm2
              IF ( m.NE.mmax ) THEN
                result = 0.0
                GO TO 100
              ELSE
                zc = (-1)**(l1+l2+m)*WIG3J(l1,l2,l,mm1,mm2,-m)
                IF ( mm1.LE.0 ) THEN
                  alpha = 1.0
                ELSE IF ( mm2.LE.0 ) THEN
                  alpha = 1.0
                ELSE
                  alpha = (1+ip1+ip2-ip12)*rt
                END IF
              END IF
            ELSE
              mm = mm1 - mm2
              zc = (-1)**(l1+l2+mm)*WIG3J(l1,l2,l,mm1,-mm2,-mm)
              IF ( mm1.EQ.0 ) THEN
                alpha = 1.0
              ELSE IF ( mm2.EQ.0 ) THEN
                alpha = 1.0
              ELSE IF ( mm.LT.0 ) THEN
                alpha = (1+ip1+ip2-ip12)*ip1*((-1)**mm1)*rt
              ELSE IF ( mm.EQ.0 ) THEN
                alpha = (-1)**mm1
              ELSE
                alpha = (1+ip1+ip2-ip12)*ip2*((-1)**mm2)*rt
              END IF
            END IF
            fa = (l1+l1+1)*(l2+l2+1)/((l+l+1)*12.566371E0)
            fas = SQRT(fa)*xl
            result = alpha*fas*zc*WIG3J0(l1,l2,l)*(-1)**(l1+l2)
          END IF
        ELSE
          result = 0.0
        END IF
      END IF
 
 100  RETURN
      END
 
      FUNCTION RINT(l,lp,n,xk)
      IMPLICIT REAL*4 (A-H, O-Z)
 
      INCLUDE 'par.h'
 
      COMMON /GAMDAT/ gamma(NGAM) , gamaf(NGAM) , gamahf(NGAM)
c
c          if n=0   r-integration for dipol term
c              =1   r-integration for quadrupole term
c              =2   r-integration for octopole term
c
      IF ( n.EQ.2 ) THEN
        is = (l+lp+2)/2 + 1
        RINT = 3.141592654E0*gamaf(is-2)
     &         /(8.0*gamaf(is)*gamahf(is-l)*gamahf(is-lp))
      ELSE IF ( n.EQ.3 ) THEN
        is = (l+lp+3)/2 + 1
        rin = gamaf(is-3)/(gamaf(is)*gamahf(is-l)*gamahf(is-lp))
        RINT = rin*3.141592654E0*xk/8.0
      ELSE
        x1 = l*(l+1) - lp*(lp+1)
        snx = SIN((l-lp)*1.570796327E0)
        RINT = snx/(x1*xk)
      END IF
 
      RETURN
      END

