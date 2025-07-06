      SUBROUTINE BORNKS(l,lp,m,mp,iq,xk,bornk)
      IMPLICIT REAL*4 (A-H, O-Z)
 
c
c      This subroutines calculate born k-matrix elements using
c      dipole quadrupole and octopole moments.
c      Note that, if any of the moments is zero,the calculations
c      are not performed to save computer time
c
      COMMON /INPPAR/ ipol , iprint , nt , nscat , mol ,
     &                lmax , lmaxb , lbig ,
     &                alpha0 , alpha2 , dpol , dpolex , quad , quadex ,
     &                dpolt , dpole , q1 , q2 , oct1 , oct2
 
      bornk = 0.0
      pi4 = 12.56637061E0
      a1 = SQRT(pi4/3.0)*dpolt
      a2 = SQRT(pi4/5.0)*q1
      a3 = SQRT(pi4/15.0)*q2
      a4 = SQRT(pi4/7.0)*oct1
      a5 = SQRT(pi4*12.0/35.0)*oct2
      r1 = 0.0E0
      r2 = 0.0E0
      r3 = 0.0E0
      rs1 = 0.0
      rs2 = 0.0
      rs3 = 0.0E0
      rs4 = 0.0E0
      rs5 = 0.0E0
      lmlp = IABS(l-lp)
      IF ( mol.NE.2 ) THEN
c
c     For c3v molecular point group
c
        IF ( lmlp.NE.1 .OR. ABS(dpole).EQ.0.0 ) THEN
          rs1 = 0.0
          r1 = 0.0
        ELSE
          CALL ANGINT(l,lp,1,m,mp,0,iq,iq,1,rs1)
          IF ( rs1.NE.0 ) r1 = RINT(l,lp,1,xk)
        END IF
        IF ( lmlp.NE.2 .OR. q1.EQ.0.0 ) THEN
          rs2 = 0.0
          r2 = 0.0
        ELSE
          CALL ANGINT(l,lp,2,m,mp,0,iq,iq,1,rs2)
          IF ( rs2.NE.0. ) r2 = RINT(l,lp,2,xk)
        END IF
        IF ( lmlp.NE.3 .OR. oct1.EQ.0.0 ) THEN
          rs3 = 0.0
          r3 = 0.0
        ELSE
          CALL ANGINT(l,lp,3,m,mp,0,iq,iq,1,rs3)
          IF ( rs3.NE.0 ) r3 = RINT(l,lp,3,xk)
        END IF
        IF ( lmlp.NE.3 .OR. oct2.EQ.0.0 ) THEN
          rs4 = 0.0
          r3 = 0.0
        ELSE
          CALL ANGINT(l,lp,3,m,mp,3,iq,iq,-1,rs4)
          IF ( rs3.EQ.0 .AND. rs4.NE.0 ) r3 = RINT(l,lp,3,xk)
        END IF
        x1 = a1*rs1*r1
        x2 = a2*rs2*r2
        x3 = (a4*rs3+a5*rs4)*r3
      ELSE
c
c     For C2v molecular point group
c
        IF ( (lmlp.NE.1) .OR. ABS(dpole).EQ.0.0 ) THEN
          rs1 = 0.0
          r1 = 0.0
        ELSE
          CALL ANGINT(l,lp,1,m,mp,0,iq,iq,1,rs1)
          IF ( rs1.NE.0. ) r1 = RINT(l,lp,1,xk)
        END IF
        IF ( lmlp.NE.2 .OR. q1.EQ.0.0 ) THEN
          rs2 = 0.0
          r2 = 0.0
        ELSE
          CALL ANGINT(l,lp,2,m,mp,0,iq,iq,1,rs2)
          IF ( rs2.NE.0. ) r2 = RINT(l,lp,2,xk)
        END IF
        IF ( lmlp.NE.2 .OR. q2.EQ.0.0 ) THEN
          rs3 = 0.0
          r2 = 0.0
        ELSE
          CALL ANGINT(l,lp,2,m,mp,2,iq,iq,1,rs3)
          IF ( rs2.EQ.0 .AND. rs3.NE.0 ) r2 = RINT(l,lp,2,xk)
        END IF
        IF ( lmlp.NE.3 .OR. oct1.EQ.0.0 ) THEN
          rs4 = 0.0
          r3 = 0.0
        ELSE
          CALL ANGINT(l,lp,3,m,mp,0,iq,iq,1,rs4)
          IF ( rs4.NE.0 ) r3 = RINT(l,lp,3,xk)
        END IF
        IF ( lmlp.NE.3 .OR. oct2.EQ.0.0 ) THEN
          rs5 = 0.0
          r3 = 0.0
        ELSE
          CALL ANGINT(l,lp,3,m,mp,2,iq,iq,1,rs5)
          IF ( rs4.EQ.0 .AND. rs5.NE.0 ) r3 = RINT(l,lp,3,xk)
        END IF
        x1 = a1*rs1*r1
        x2 = (a2*rs2+a3*rs3)*r2
        x3 = (a4*rs4+a5*rs5)*r3
      END IF
      bornk = 2.0*xk*(x1+x2+x3)

      RETURN
      END

