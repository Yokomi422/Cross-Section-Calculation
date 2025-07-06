      SUBROUTINE GAC(jj,mj,g)
      IMPLICIT REAL*4 (A-H, O-Z)
 
      DIMENSION y(2) , x(2)
 
      COMMON /ASYDAT/ ai(11,6,2) , af(11,6,2) , egi(11) , egf(11)
      COMMON /INPPAR/ ipol , iprint , nt , nscat , mol ,
     &                lmax , lmaxb , lbig ,
     &                alpha0 , alpha2 , dpol , dpolex , quad , quadex ,
     &                dpolt , dpole , q1 , q2 , oct1 , oct2
      COMMON /JKTDAT/ jmin, jmax , mnlamda , mxlamda ,
     &                jint , jfnl , itai , itaf , numti , numtf

 
      IF ( mol.EQ.1 ) THEN
c
c     Calculation of g-coefficients for symmetric top molecules
c
        jftaum = -itaf
        ks = itai + jftaum + mj
        IF ( ks.EQ.0 ) THEN
          IF ( mj.EQ.0 .AND. itai.EQ.0 ) THEN
            IF ( IPARITY(jint+jfnl+jj).NE.1 ) GO TO 50
            g1 = WIG3J0(jint,jfnl,jj)
            g = SQRT(FLOAT(2*jfnl+1))*g1
            RETURN
          END IF
          g1 = WIG3J(jint,jfnl,jj,itai,jftaum,mj)
          g = SQRT(FLOAT(2*jfnl+1))*g1
          RETURN
        END IF
 50     g = 0.0
 
        RETURN
      ELSE IF ( mol.EQ.0 ) THEN
c
c     Calculation of g-coefficients for spherical top molecules
c
        xj21 = FLOAT(jint+jint+1)
        g = SQRT(1.0/xj21)
 
        RETURN
      ELSE IF ( mol.NE.2 ) THEN
c
c     Calculation of the g-coefficients for linear molecules
c
        g = WIG3J0(jint,jfnl,jj)
 
        RETURN
      END IF
 
c
c     Calculation of g-coefficients for asymmetric top molecules
c
      rst = SQRT(2.0)*2.0
      j1 = jint + 1
      j2 = jfnl + 1
      jtot = jint + jfnl + jj
      jtoth = (jtot/2)*2
      g = 0.0
c
c     Loop over k
c
      DO k = 1 , j1
        y(1) = ai(numti,k,1)
        y(2) = ai(numti,k,2)
c
c     If both a-coeffs corresponding to nu=0 and 1 are zero
c     then the g-coeffs are zero, too.
c
        IF ( y(1).NE.0 .OR. y(2).NE.0 ) THEN
 
          kk = k - 1
          km = -kk
c
c     do loop over k (prime)
c
          DO kb = 1 , j2
            x(1) = af(numtf,kb,1)
            x(2) = af(numtf,kb,2)
c
c     If both the a-coeffs corresponding to nu=0 and 1 are
c     zero,then the g-coeffs are zero, too.
c
            IF ( x(1).NE.0 .OR. x(2).NE.0. ) THEN
              kkb = kb - 1
              c1 = (-1)**kkb
              kbm = -kkb
              fact = 2.0
              IF ( kk.EQ.0 .AND. kkb.EQ.0 ) fact = 4.0
              IF ( kk.EQ.0 .AND. kkb.GT.0 ) fact = rst
              IF ( kk.GT.0 .AND. kkb.EQ.0 ) fact = rst
              k1 = km + kkb + mj
              IF ( k1.NE.0 ) THEN
                w1 = 0.0
              ELSE IF ( km.NE.0 .OR. kkb.NE.0 ) THEN
                w1 = WIG3J(jint,jfnl,jj,km,kkb,mj)
              ELSE IF ( jtoth.EQ.jtot ) THEN
                w1 = WIG3J0(jint,jfnl,jj)
              ELSE
                w1 = 0.0
              END IF
              k2 = kk + kkb + mj
              IF ( k2.NE.0 ) THEN
                w2 = 0.0
              ELSE IF ( kk.NE.0 .OR. kkb.NE.0 ) THEN
                w2 = WIG3J(jint,jfnl,jj,kk,kkb,mj)
              ELSE IF ( jtoth.EQ.jtot ) THEN
                w2 = WIG3J0(jint,jfnl,jj)
              ELSE
                w2 = 0.0
              END IF
              k3 = km + kbm + mj
              IF ( k3.NE.0 ) THEN
                w3 = 0.0
              ELSE IF ( kk.NE.0 .OR. kkb.NE.0 ) THEN
                w3 = WIG3J(jint,jfnl,jj,km,kbm,mj)
              ELSE IF ( jtoth.EQ.jtot ) THEN
                w3 = WIG3J0(jint,jfnl,jj)
              ELSE
                w3 = 0.0
              END IF
              k4 = kk + kbm + mj
              IF ( k4.NE.0 ) THEN
                w4 = 0.0
              ELSE IF ( km.NE.0 .OR. kkb.NE.0 ) THEN
                w4 = WIG3J(jint,jfnl,jj,kk,kbm,mj)
              ELSE IF ( jtoth.EQ.jtot ) THEN
                w4 = WIG3J0(jint,jfnl,jj)
              ELSE
                w4 = 0.0
              END IF
              gn = 0.0
c
c     do loop over nu
c
              DO nu1 = 1 , 2
                nuu1 = nu1 - 1
                c2 = (-1)**nuu1
                cw1 = c2*w2 + w1
                cw2 = w3 + c2*w4
                a1 = ai(numti,k,nu1)
c
c     do loop over nu (prime)
c
                DO nu2 = 1 , 2
                  nuu2 = nu2 - 1
                  c3 = (-1)**nuu2
                  a2 = x(nu2)
                  gn = gn + a1*c1*a2*(cw1+c3*cw2)
                END DO
              END DO
              gn = gn/fact
              g = g + gn
            END IF
          END DO
        END IF
      END DO
 
      RETURN
      END
 
