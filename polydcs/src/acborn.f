      SUBROUTINE ACBORN(m,l1,l2,en,accborn)
      IMPLICIT REAL*4 (A-H, O-Z)
 
c
c                                                             m
c    This subroutine calculated the unitarized Born T matrix T
c                                                              l.l'
c
      PARAMETER (PI=3.1415927)

      COMMON /INPPAR/ ipol , iprint , nt , nscat , mol ,
     &                lmax , lmaxb , lbig ,
     &                alpha0 , alpha2 , dpol , dpolex , quad , quadex ,
     &                dpolt , dpole , q1 , q2 , oct1 , oct2
 
      acb1 = 0.0
      acb2 = 0.0
      acb3 = 0.0
c
c                             m
c    This is for the element T
c                             l,l+1
c
 
      IF ( ABS(l1-l2).EQ.1 ) THEN
        c1 = dpol*IPARITY(l1)
        xc2 = FLOAT(3*(MAX(l1,l2)))
        c2 = 1.0/SQRT(xc2)
        IF ( m.EQ.0 ) THEN
          c3 = -1.0*SQRT(3.0)*WIG3J0(l1,l2,1)
        ELSE
          c3 = -1.0*SQRT(3.0)*IPARITY(m)*WIG3J(l1,l2,1,-m,m,0)
        END IF
        c4 = 2.0
        accborn = c1*c2*c3*c4
      END IF
c
c                             m
c    This is for the element T
c                             l,l
c
      IF ( ABS(l1-l2).EQ.0 ) THEN
c
c    Quadrupole contribution
c
        c1 = quad*IPARITY(l1)
        xc2 = SQRT(FLOAT(5*(2*l1+3)*(2*l1-1)*(l1+1)*l1))
        xc22 = SQRT(FLOAT(2*l1+1))
        c2 = -1.0*xc22/xc2
        IF ( m.EQ.0 ) THEN
          c3 = IPARITY(m)*SQRT(5.0)*WIG3J0(l1,l1,2)
        ELSE
          c3 = IPARITY(m)*SQRT(5.0)*WIG3J(l1,l1,2,m,-m,0)
        END IF
        c4 = 2.0*SQRT(en)
        acb1 = c1*c2*c3*c4
c
c   Spherical polarizability contribution
c
        c1 = alpha0*IPARITY(l1)
        xc22 = FLOAT((2*l1+3)*(2*l1-1))
        xc2 = SQRT(FLOAT(2*l1+1))
        c2 = 1.0/(xc2*xc22)
        IF ( m.EQ.0 ) THEN
          c3 = WIG3J0(l1,l1,0)
        ELSE
          c3 = IPARITY(m)*WIG3J(l1,l1,0,m,-m,0)
        END IF
        c4 = -2.0*PI*en
        acb2 = c1*c2*c3*c4
c
c   Non spherical polarizability contribution
c
        c1 = alpha2*IPARITY(l1)
        xc2 = FLOAT((2*l1+3)*(2*l1-1))
        xc2a = (SQRT(xc2))**3
        xc22 = FLOAT(4*l1*(l1+1))/FLOAT(2*l1+1)
        xc22a = SQRT(xc22)
        c2 = xc22a/xc2a
        IF ( m.EQ.0 ) THEN
          c3 = IPARITY(m)*SQRT(5.0)*WIG3J0(l1,l1,2)
        ELSE
          c3 = IPARITY(m)*SQRT(5.0)*WIG3J(l1,l1,2,m,-m,0)
        END IF
        c4 = PI*en/SQRT(5.0)
        acb3 = c1*c2*c3*c4
 
        accborn = acb1 + acb2 + acb3
      END IF
c
c                             m
c    This is for the element T
c                             l,l+2
c
 
      IF ( ABS(l1-l2).EQ.2 ) THEN
        ll = MAX(l1,l2)
c
c    Quadrupole contribution
c
        c1 = quad*IPARITY(l1)
        xc2 = FLOAT(30*(2*ll-1)*ll*(ll-1))
        c2 = SQRT(1.0/xc2)
        IF ( m.EQ.0 ) THEN
          c3 = SQRT(5.0)*WIG3J0(l1,l2,2)
        ELSE
          c3 = IPARITY(m)*SQRT(5.0)*WIG3J(l1,l2,2,m,-m,0)
        END IF
        c4 = 2.0*SQRT(en)
        acb1 = c1*c2*c3*c4
c
c   Non spherical polarizability contribution
c
        c1 = alpha2*IPARITY(l1)
        xc2a = FLOAT((2*ll-1)**3)
        xc2aa = FLOAT(6*ll*(ll-1))
        xc2 = SQRT(xc2aa/xc2a)
        xc22 = FLOAT(2*(2*l1+1)*(2*ll-3))
        c2 = xc2/xc22
        IF ( m.EQ.0 ) THEN
          c3 = SQRT(5.0)*WIG3J0(l1,l1,2)
        ELSE
          c3 = IPARITY(m)*SQRT(5.0)*WIG3J(l1,l1,2,m,-m,0)
        END IF
        c4 = -PI*en/SQRT(5.0)
        acb2 = c1*c2*c3*c4
 
        accborn = acb1 + acb2
      END IF
 
      RETURN
      END
 
 
