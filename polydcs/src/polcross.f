      SUBROUTINE POLCROSS(en4)
      IMPLICIT REAL*4 (A-H, O-Z)
 
      INCLUDE 'par.h'

      LOGICAL done
      CHARACTER*6  stat
      CHARACTER*80 blmfn , kmatfn

      DIMENSION ama(NSYM)

      COMMON /BLMCOF/ br(NSYM,NPW,NPW) , bi(NSYM,NPW,NPW)
      COMMON /INPCHR/ stat(NSYM) , blmfn , kmatfn
      COMMON /INPMAT/ akm(NSYM,NPW,NPW) ,
     &                mult(NSYM) , ksiz(NSYM) , ksizb(NSYM) ,
     &                lk(NSYM,NPW) , num(NSYM,NPW) , ml(NSYM,NPW,NPW)
      COMMON /INPPAR/ ipol , iprint , nt , nscat , mol ,
     &                lmax , lmaxb , lbig ,
     &                alpha0 , alpha2 , dpol , dpolex , quad , quadex ,
     &                dpolt , dpole , q1 , q2 , oct1 , oct2
      COMMON /TMAT  / tmr(NSYM,NPW,NPW) , tmi(NSYM,NPW,NPW)

      DATA uunit , pi4/3.5712026E0 , 12.566371E0/
c
c     Calculate the cross sections for polyatomic molecules
c     using T-matrix elements
c
 
      fict = pi4/(en4)
      fict = fict/uunit
c
c     do loop over small l
c
      lma1 = lmax
      DO l = 1 , lma1
        l11 = l - 1
c
c     do loop over l'
c
        DO lp = 1 , lma1
          l22 = lp - 1
c
c     do loop over the molecular IR's
c
          DO k = 1 , nscat
            ip = 0
            nh = 0
            ip1 = 0
            nh1 = 0
            done = .FALSE.
            ks = ksiz(k)
            DO k1 = 1 , ks
              IF ( .NOT.done ) THEN
                IF ( lk(k,k1).EQ.l11 ) THEN
                  ip = k1
                  DO k2 = ip , ks
                    IF ( lk(k,k2).EQ.l11 ) nh = nh + 1
                  END DO
                  done = .TRUE.
                END IF
              END IF
            END DO
            done = .FALSE.
            DO k1 = 1 , ks
              IF ( .NOT.done ) THEN
                IF ( lk(k,k1).EQ.l22 ) THEN
                  ip1 = k1
                  DO k2 = ip1 , ks
                    IF ( lk(k,k2).EQ.l22 ) nh1 = nh1 + 1
                  END DO
                  done = .TRUE.
                END IF
              END IF
            END DO
            IF ( nh.NE.0 .AND. nh1.NE.0 ) THEN
c
c     do loop over h and h'
c
              DO ih = ip , ip + nh - 1
                DO ih1 = ip1 , ip1 + nh1 - 1
                  m1 = num(k,ih)
                  m2 = num(k,ih1)
                  t1r = tmr(k,ih,ih1)
                  t1i = tmi(k,ih,ih1)
c
c     do loop over m
c 
                  DO m = 1 , m1
                    b1r = br(k,ih,m)
                    b1i = bi(k,ih,m)
                    mb1 = ml(k,ih,m)
c 
c     do loop over m'
c 
                    DO mp = 1 , m2
                      b2r = br(k,ih1,mp)
                      b2i = bi(k,ih1,mp)
                      mb2 = ml(k,ih1,mp)
                      c1 = IPARITY(mb1+mb2)
                      c1 = fict*c1
                      f1 = b1r*t1r - b1i*t1i
                      f2 = b1i*t1r + b1r*t1i
 
                      amrpa = (b2r*f1-b2i*f2)
                      amipa = (b2i*f1+b2r*f2)
                      ampa = amrpa**2 + amipa**2
                      ama(k) = ama(k) + ampa*c1
                    END DO
                  END DO
                END DO
              END DO
            END IF
          END DO
        END DO
      END DO
 
      WRITE (6,99000)

      amat = 0
      DO k = 1 , nscat
        amat = amat + ama(k)
          WRITE (6,99001) k , stat(k), ama(k)
      END DO

      WRITE (6,99002) amat 

99000 FORMAT (/2x,'Integral Cross Sections (CS) from T-matrix',
     &        ' in units of 10**-16 cm**2',/2x,68('-'),/)
99001 FORMAT (2x,'Contribution of the IR #',i2,' (',a6,')',
     &        ' to CS = ',f12.8)
99002 FORMAT (//2x,'Integral Cross Section (CS) for all IRs   = ',
     &        f12.8,//)

      RETURN
      END

