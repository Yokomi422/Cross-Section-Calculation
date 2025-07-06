      SUBROUTINE BORNTM(xk)
      IMPLICIT REAL*4 (A-H, O-Z)
 
      INCLUDE 'par.h'
 
c
c     This subroutine calculates Born T matrix for all molecules
c     using the subroutine bornks that calculate the Born K-matrix
c     elements for polyatomic molecules and subroutine borntel
c     that calculate the Born T-matrix elements for linear molecules.
c
 
      LOGICAL done
      DIMENSION bornkm(NSYM,NPW,NPW)
      DIMENSION iq(4)
 
      COMMON /INPMAT/ akm(NSYM,NPW,NPW) ,
     &                mult(NSYM) , ksiz(NSYM) , ksizb(NSYM) ,
     &                lk(NSYM,NPW) , num(NSYM,NPW) , ml(NSYM,NPW,NPW)
      COMMON /INPPAR/ ipol , iprint , nt , nscat , mol ,
     &                lmax , lmaxb , lbig ,
     &                alpha0 , alpha2 , dpol , dpolex , quad , quadex ,
     &                dpolt , dpole , q1 , q2 , oct1 , oct2
      COMMON /TMAT  / tmr(NSYM,NPW,NPW) , tmi(NSYM,NPW,NPW)
      COMMON /TMATB / borntr(NSYM,NPW,NPW) , bornti(NSYM,NPW,NPW)
      COMMON /TMATBL/ borntrl(0:NLAM,NPS,NPS) , borntil(0:NLAM,NPS,NPS)
 
      DATA iq/1 , -1 , -1 , 1/
 
c
c     This is for linear molecules
c
      IF ( mol.GT.3 ) THEN
        llmax = 3*lbig - 2
        DO l = 1 , llmax
          l1 = l - 1
          DO lp = 1 , llmax
            l2 = lp - 1
            jmn = IABS(l1-l2)
            jmx = l1 + l2
            DO jt = jmn , jmx
              CALL BORNTEL(l1,l2,jt,xk,bortr,borti)
              borntrl(jt,l,lp) = bortr
              borntil(jt,l,lp) = borti
            END DO
          END DO
        END DO
        GO TO 100
      END IF
c
c*************end for linear molecules************************
c
      DO l = 1 , lmaxb
        DO lp = 1 , lmaxb
          DO i = 1 , nscat
            iqd = iq(i)
            ks = ksizb(i)
            ip = 0
            nh = 0
            ip1 = 0
            nh1 = 0
            l11 = l - 1
            done = .FALSE.
            DO k1 = 1 , ks
              IF ( .NOT.done ) THEN
                ll = lk(i,k1)
                IF ( l11.EQ.ll ) THEN
                  ip = k1
                  DO k2 = ip , ks
                    IF ( lk(i,k2).EQ.l11 ) THEN
                      nh = nh + 1
c                     m1 = ml(i,k1,nh)
                    END IF
                  END DO
                  done = .TRUE.
                END IF
              END IF
            END DO
            l22 = lp - 1
            done = .FALSE.
            DO k1 = 1 , ks
              IF ( .NOT.done ) THEN
                llp = lk(i,k1)
                IF ( l22.EQ.llp ) THEN
                  ip1 = k1
                  DO k2 = ip1 , ks
                    IF ( lk(i,k2).EQ.l22 ) THEN
                      nh1 = nh1 + 1
c                     m2 = ml(i,k2,nh1)
                    END IF
                  END DO
                  done = .TRUE.
                END IF
              END IF
            END DO
            IF ( nh.NE.0 .AND. nh1.NE.0 ) THEN
              DO ih = ip , ip + nh - 1
                DO ih1 = ip1 , ip1 + nh1 - 1
c                 ji = ih - ip + 1
c                 jl = ih1 - ip1 + 1
                  mm1 = ml(i,ih,1)
                  mm2 = ml(i,ih1,1)
                  IF ( ih.LE.ih1 ) THEN
                    CALL BORNKS(ll,llp,mm1,mm2,iqd,xk,bornk)
                    bornkm(i,ih,ih1) = bornk
                    IF ( ll.GE.lmax .OR. llp.GE.lmax ) THEN
                      akm(i,ih,ih1) = bornk
                      akm(i,ih1,ih) = bornk
                    END IF
                  END IF
                END DO
              END DO
              DO j = 1 , ks
                DO k = j , ks
                  IF ( j.NE.k ) bornkm(i,k,j) = bornkm(i,j,k)
                END DO
              END DO
            END IF
          END DO
        END DO
      END DO
c
c     Call tmatrx to generate born t-matrix
c     (Unitarised Born Approximation)
c
      CALL TMATRX(bornkm,borntr,bornti,ksizb)
      IF ( (mol.EQ.2 .OR. mol.EQ.1) .AND. (ipol.EQ.1) )
     &     CALL TMATRX(akm,tmr,tmi,ksizb)
 
 100  RETURN
 
c
c     Format statements
c
99001 FORMAT (//5x,'born(unitarised) k-matrix for state=',i1,2x,a8,//)
99002 FORMAT (//5x,'real born t-matrix for above born k-matrix',//)
99003 FORMAT (//5x,'imag born t-matrix for above born k-matrix',//)
      END
