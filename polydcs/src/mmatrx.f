      SUBROUTINE MMATRX
      IMPLICIT REAL*4 (A-H, O-Z)
 
      INCLUDE 'par.h'
 
      LOGICAL done
      COMMON /BLMCOF/ br(NSYM,NPW,NPW) , bi(NSYM,NPW,NPW)
      COMMON /INPMAT/ akm(NSYM,NPW,NPW) ,
     &                mult(NSYM) , ksiz(NSYM) , ksizb(NSYM) ,
     &                lk(NSYM,NPW) , num(NSYM,NPW) , ml(NSYM,NPW,NPW)
      COMMON /INPPAR/ ipol , iprint , nt , nscat , mol ,
     &                lmax , lmaxb , lbig ,
     &                alpha0 , alpha2 , dpol , dpolex , quad , quadex ,
     &                dpolt , dpole , q1 , q2 , oct1 , oct2
      COMMON /JKTDAT/ jmin, jmax , mnlamda , mxlamda ,
     &                jint , jfnl , itai , itaf , numti , numtf
      COMMON /MMAT  / amr(0:MAXN) , ami(0:MAXN) ,
     &                amrb(0:MAXN) , amib(0:MAXN)
      COMMON /MMATL / ar(11,NPS,NPS) , aii(11,NPS,NPS) ,
     &                arb(11,NPS,NPS) , aib(11,NPS,NPS)
      COMMON /TMAT  / tmr(NSYM,NPW,NPW) , tmi(NSYM,NPW,NPW)
      COMMON /TMATB / borntr(NSYM,NPW,NPW) , bornti(NSYM,NPW,NPW)
      COMMON /TMATBL/ borntrl(0:NLAM,NPS,NPS) , borntil(0:NLAM,NPS,NPS)

C----------------------------------------------------------------------
C                                          j,mj
C        Calculate the M-matrix elements  M      for poliatomic
C                                          l,l'
C                                j
C              molecules  and   M     for linear molecules.
C                                l,l'
C----------------------------------------------------------------------

      DO i = 0 , MAXN
        amr(i) = 0.0E0
        ami(i) = 0.0E0
      END DO
      IF ( mol.EQ.4 .AND. ipol.NE.0 ) THEN
        DO l = lmax + 1 , 3*lmaxb - 2
          DO lp = 1 , 3*lmaxb - 2
            DO j = 1 , 3
              ar(j,l,lp) = -borntrl(j-1,l,lp)
              aii(j,l,lp) = -borntil(j-1,l,lp)
              arb(j,l,lp) = -borntrl(j-1,l,lp)
              aib(j,l,lp) = -borntil(j-1,l,lp)
              IF ( lp.LT.(lmax+1) ) THEN
                ar(j,lp,l) = ar(j,l,lp)
                aii(j,lp,l) = aii(j,l,lp)
                arb(j,lp,l) = arb(j,l,lp)
                aib(j,lp,l) = aib(j,l,lp)
              END IF
            END DO
          END DO
        END DO
      END IF
 
      IF ( (mol.EQ.2 .OR. mol.EQ.1) .AND. ipol.EQ.1 ) THEN
        lma1 = lmaxb
      ELSE
        lma1 = lmax
      END IF
 
      jmmin = jmin
      jmmax = jmax
      IF ( mol.GT.3 ) THEN
        jmmin = 1
        jmmax = 11
      END IF
c
c     do loop over small l
c
      DO l = 1 , lma1
        l11 = l - 1
        lx = l11*lma1*jmax
c
c     do loop over l'
c
        DO lp = 1 , lma1
          l22 = lp - 1
          l1 = IABS(l11-l22)
          l2 = l11 + l22
          lx1 = l22*jmax
c
c     do loop over j
c 
          DO j = jmmin , jmmax
            jj = j - 1
            IF ( jj.GE.l1 .AND. jj.LE.l2 ) THEN
              lhi = (l11+l22+jj)/2
              lfp = (l11+l22+jj) - 2*lhi
              lldj = l11 + l22 + jj
              jx = 2*jj + 1
c 
c***** Begin code for LINEAR molecule (MOL=4,5) *****
c 
              IF ( mol.GT.3 ) THEN
 
                lmmx = MIN(l,lp,mxlamda+1)
 
                DO lm = mnlamda + 1 , lmmx
                  mm = lm - 1
                  c1 = 0.0
                  IF ( mm.NE.0 .OR. lfp.EQ.0 ) THEN
                    IF ( mm.EQ.0 ) THEN
                      c1 = WIG3J0(l11,l22,jj)
                    ELSE IF ( ABS(mm).LE.l11 .AND. ABS(mm).LE.l2 ) THEN
                      c1 = WIG3J(l11,l22,jj,mm,-mm,0)
                      c1b = WIG3J(l11,l22,jj,-mm,mm,0)
                      c1 = c1 + c1b
                    END IF
                    c11 = FLOAT(jx)
                    ar(j,l,lp) = ar(j,l,lp) + c1*(-1)**mm*tmr(lm,l,lp)
     &                           *SQRT(c11)
                    aii(j,l,lp) = aii(j,l,lp) + c1*(-1)**mm*tmi(lm,l,lp)
     &                            *SQRT(c11)
                  END IF
                END DO
                IF ( ipol.NE.0 ) THEN
                  arb(j,l,lp) = -borntrl(jj,l,lp)
                  aib(j,l,lp) = -borntil(jj,l,lp)
                END IF
              END IF
c
c***** End code for LINEAR molecule (MOL=4,5) *****
c
              IF ( mol.LE.2 ) THEN
c
c     do loop over mj
c 
                DO mjj = 1 , jx
                  mj = j - mjj
                  mjjj = (-1)**mj
                  CALL GAC(jj,mj,g)
                  n = (lx+lx1+jj)*jx + mjj
                  IF ( n.GT.MAXN ) THEN
                    WRITE (6,99001) n , MAXN
                    STOP
                  END IF
                  amr(n) = 0.0E0
                  ami(n) = 0.0E0
                  IF ( ipol.NE.0 ) THEN
                    amrb(n) = 0.0E0
                    amib(n) = 0.0E0
                  END IF
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
                    IF ( ipol.NE.0 ) ks = ksizb(k)
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
                          IF ( ipol.NE.0 ) THEN
                            trb = borntr(k,ih,ih1)
                            tib = bornti(k,ih,ih1)
                          END IF
c
c               do loop over m
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
                              mbj = mb1 + mb2 + mj
                              IF ( mbj.EQ.0 ) THEN
                                IF ( mb1.NE.0 .OR. mb2.NE.0 ) THEN
                                  c1 = WIG3J(l11,l22,jj,mb1,mb2,mj)
                                ELSE
                                  IF ( (lldj/2)*2.NE.lldj ) GO TO 2
                                  c1 = WIG3J0(l11,l22,jj)
                                END IF
                                f1 = b1r*t1r - b1i*t1i
                                f2 = b1i*t1r + b1r*t1i
                                IF ( ipol.NE.0 ) THEN
                                  f1b = b1r*trb - b1i*tib
                                  f2b = b1i*trb + b1r*tib
                                END IF
                                amr(n) = amr(n) + (b2r*f1-b2i*f2)*c1
                                ami(n) = ami(n) + (b2i*f1+b2r*f2)*c1
                                IF ( ipol.NE.0 ) THEN
                                  amrb(n) = amrb(n) + (b2r*f1b-b2i*f2b)
     &                              *c1
                                  amib(n) = amib(n) + (b2i*f1b+b2r*f2b)
     &                              *c1
                                END IF
                              END IF
 2                          END DO
                          END DO
                        END DO
                      END DO
                    END IF
                  END DO
                  IF ( mol.EQ.1 ) THEN
                    amr(n) = amr(n)*g
                    ami(n) = ami(n)*g
                    IF ( ipol.NE.0 ) THEN
                      amrb(n) = amrb(n)*g
                      amib(n) = amib(n)*g
                    END IF
                  ELSE
                    amr(n) = amr(n)*g*mjjj
                    ami(n) = ami(n)*g*mjjj
                    IF ( ipol.NE.0 ) THEN
                      amrb(n) = amrb(n)*g*mjjj
                      amib(n) = amib(n)*g*mjjj
                    END IF
                  END IF
                END DO
              END IF
            END IF
          END DO
        END DO
      END DO 

      RETURN

99001 FORMAT (//,' *** POLYDCS ERROR in MMATRX - - PROGRAM STOP ***',/,
     &        ' N = ',i5,' greater than MAXN = ',i5)

      END
