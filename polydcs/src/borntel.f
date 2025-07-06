      SUBROUTINE BORNTEL(l,lp,jt,xk,bortr,borti)
      IMPLICIT REAL*4 (A-H, O-Z)
 
      PARAMETER (PI=3.141592,PI2=6.2831853)
c
c     This subroutine calculate born t-matrix elements using
c     dipole, quadrupole and octupole moments for linear molecules.
c     borti represents the T matrix element with index l, l',jt
c
 
      COMMON /INPPAR/ ipol , iprint , nt , nscat , mol ,
     &                lmax , lmaxb , lbig ,
     &                alpha0 , alpha2 , dpol , dpolex , quad , quadex ,
     &                dpolt , dpole , q1 , q2 , oct1 , oct2
 
      bortr = 0.0E0
      borti = 0.0E0
      tm1 = 0.0
      tm2 = 0.0
      tm3 = 0.0
      tm4 = 0.0
 
      l1 = l + 1
      lm1 = l - 1
      l2 = l + 2
      lm2 = l - 2
      lmaxi = MAX(l,lp)
c
c     Spherical polarizability contribution
c
      IF ( alpha0.NE.0.0 .AND. jt.EQ.0 ) THEN
        IF ( l.EQ.lp .AND. l.NE.0 ) THEN
          xn1 = FLOAT((2*l+3)*(2*l-1))
          xn2 = SQRT(FLOAT((2*l+1)))
          tm1 = alpha0*xk**2*PI2*(-1)**(l1)/(xn1*xn2)
        END IF
      END IF
c
c     Dipole contribution
c
      IF ( dpol.NE.0.0 .AND. jt.EQ.1 ) THEN
        IF ( lp.EQ.l1 .OR. lp.EQ.lm1 ) THEN
          xn3 = SQRT(FLOAT(3*lmaxi))
          tm2 = 2*dpol*(-1)**lmaxi/xn3
        END IF
      END IF
c
c     Quadrupole contribution
c
      IF ( quad.NE.0.0 .AND. jt.EQ.2 ) THEN
        xn5 = 0.0E0
        IF ( lp.EQ.l2 .OR. lp.EQ.lm2 .OR. lp.EQ.l ) THEN
          IF ( l.NE.lp ) THEN
            n4a = 30*(2*lmaxi-1)
            n4aa = lmaxi*(lmaxi-1)
            xn4a = SQRT(FLOAT(n4a))
            xn4aa = SQRT(FLOAT(n4aa))
            xn4 = xn4a*xn4aa
            xn5 = 1/xn4
          ELSE
            xn6 = SQRT(FLOAT(2*l+1))
            xn7a = SQRT(FLOAT(5*(2*l+3)*(2*l-1)))
            xn7aa = SQRT(FLOAT((l+1)*l))
            xn7 = xn7a*xn7aa
            xn5 = -xn6/xn7
          END IF
          IF ( xn5.NE.0.0 ) tm3 = 2*xk*xn5*quad*(-1)**l
        END IF
      END IF
c
c     Non spherical polarizability contribution
c
      IF ( alpha2.NE.0.0 .AND. jt.EQ.2 ) THEN
        IF ( lp.EQ.l2 .OR. lp.EQ.lm2 .OR. lp.EQ.l ) THEN
          IF ( lp.NE.l ) THEN
            xn8 = SQRT(FLOAT(6*lmaxi*(lmaxi-1)))
            xn9 = SQRT(FLOAT((2*lmaxi-1)**3))
            xn10 = FLOAT(2*(2*lmaxi+1)*(2*lmaxi-3))
            xn11 = xn8/(xn9*xn10)
          ELSE
            xn12 = 2*SQRT(FLOAT(l*(l+1)))
            xn13aaa = SQRT(FLOAT((2*l+3)*(2*l-1)))
            xn13aa = xn13aaa**3
            xn13a = SQRT(FLOAT(2*l+1))
            xn13 = xn13a*xn13aa
            xn11 = -xn12/xn13
          END IF
          IF ( xn11.NE.0.0 ) tm4 = PI*xk**2*alpha2*xn11*(-1)**(l1)
     &                             /SQRT(5.0E0)
        END IF
      END IF
 
      borti = (tm1+tm2+tm3+tm4)

      RETURN
      END
 
