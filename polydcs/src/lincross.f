      SUBROUTINE LINCROSS(en)
      IMPLICIT REAL*4 (A-H, O-Z)

c
c     Evaluate CS from T-matrices for a linear polar molecule
c
c     Note that in this case the lk(,) matrix is used for 2 rows only:
c     - lk(1,NPW=NSYM) = li() initial L value for each IR
c     - lk(2,NPW=NSYM) = lf() final L value for each IR
c

      INCLUDE 'par.h'

      COMMON /INPMAT/ akm(NSYM,NPW,NPW) ,
     &                mult(NSYM) , ksiz(NSYM) , ksizb(NSYM) ,
     &                lk(NSYM,NPW) , num(NSYM,NPW) , ml(NSYM,NPW,NPW)
      COMMON /INPPAR/ ipol , iprint , nt , nscat , mol ,
     &                lmax , lmaxb , lbig ,
     &                alpha0 , alpha2 , dpol , dpolex , quad , quadex ,
     &                dpolt , dpole , q1 , q2 , oct1 , oct2
      COMMON /TMAT  / tmr(NSYM,NPW,NPW) , tmi(NSYM,NPW,NPW)


      WRITE (6,99000)

      totcsa = 0.0
      totcsb = 0.0
      DO ns = 1 , nscat
 
        parcs = 0.0
        DO l = lk(1,ns) , lk(2,ns)
          DO lp = lk(1,ns) , lk(2,ns)
            cstm = tmr(ns,l,lp)*tmr(ns,l,lp) + tmi(ns,l,lp)
     &              *tmi(ns,l,lp)
            IF ( ns.NE.1 ) cstm = 2*cstm
            parcs = parcs + cstm
          END DO
        END DO

        c34 = 1.0
        IF ( ns.GT.1 ) c34 = 0.5
        csb = c34*parcs*3.1415927/en
        csa = csb/3.5712026
        totcsb = totcsb + csb
        totcsa = totcsa + csa
 
        WRITE (6,99001) ns , csb , csa

      END DO

      WRITE (6,99002) totcsb , totcsa

99000 FORMAT (/2x,'Integral Cross Sections (CS) from T-matrix',
     &        ' ',/2x,42('-'),/)
99001 FORMAT (2x,'Contribution of the IR #',i2,' to CS = ',
     &        f12.8,' Bohr**2 ',f12.8,' Ang.**2')
99002 FORMAT (//2x,'Integral Cross Section (CS) for all IRs   = ',
     &        f12.8,' Bohr**2 ',f12.8,' Ang.**2',//)

      END
