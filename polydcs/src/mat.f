      SUBROUTINE PLEG
      IMPLICIT REAL*4 (A-H, O-Z)
 
      INCLUDE 'par.h'

      COMMON /ANGDAT/ th(NANG) , pl(NLAM,NANG)
      COMMON /INPPAR/ ipol , iprint , nt , nscat , mol ,
     &                lmax , lmaxb , lbig ,
     &                alpha0 , alpha2 , dpol , dpolex , quad , quadex ,
     &                dpolt , dpole , q1 , q2 , oct1 , oct2
 
      l1 = lbig - 1
      DO i = 1 , NANG
        t = th(i)*0.017453292E0
        x = COS(t)
        a = 1.0
        b = x
        pl(1,i) = x
        DO l = 2 , l1
          l2 = l - 1
          c = ((l+l2)*x*b-l2*a)/l
          a = b
          b = c
          pl(l,i) = c
        END DO
      END DO
      RETURN
      END

 
      SUBROUTINE MA01A(a,m)
      IMPLICIT REAL*4 (A-H, O-Z)
 
      INCLUDE 'par.h'
 
c********************************************************************
c
c        matrix inversion
c
c **********************************************************************
c
      DIMENSION a(NPW,NPW) , ind(NPW) , c(NPW)  
c
c        a         the m*m matrix being inverted. overwritten on exit by
c                  the inverse matrix
c        m         the order of the a-matrix. this is equal to the
c                  dimension of the matrix (1+k**2). if m is greater
c                  than the dimension of the arrays a, ind and c
c                  recompiling is needed with larger dimension for all
c                  these arrays
c        ind,c     private arrays used by the subroutine. they must
c                  have at least dimension m
c
c
c        if the matrix is 1*1 its inverse is the inverse of its element
c        (instruction 316)
c
      IF ( m.EQ.1 ) THEN
        a(1,1) = 1.E0/a(1,1)
      ELSE
        amax = 0.E0
c
c        find the first pivotal element and store the corresponding row
c        number in i4. ind defines the order of the rows of the original
c        a-matrix before row interchange
c
        DO i = 1 , m
          ind(i) = i
          IF ( ABS(a(i,1)).GT.amax ) THEN
            amax = ABS(a(i,1))
            i4 = i
          END IF
        END DO
        mm = m - 1
c
c        each time through the following loop the a-matrix is
c        reduced by one
c
        DO j = 1 , mm
c
c        interchange the i4-th and the j-th rows of the a-matrix and
c        store order in ind if (i4.ne.j)
c
          IF ( i4.GT.j ) THEN
            isto = ind(j)
            ind(j) = ind(i4)
            ind(i4) = isto
            DO k = 1 , m
              sto = a(i4,k)
              a(i4,k) = a(j,k)
              a(j,k) = sto
            END DO
          END IF
c
c        the j-th row now contains the pivotal element in the j-th
c        position eliminate the j-th element from each of the remaining
c        rows of the a-matrix and store the multipliers in the lower
c        triangle.
c
          amax = 0.E0
          j1 = j + 1
          DO i = j1 , m
            a(i,j) = a(i,j)/a(j,j)
            DO k = j1 , m
              a(i,k) = a(i,k) - a(i,j)*a(j,k)
              IF ( k.LE.j1 ) THEN
c
c       find the next pivotal element and store the corresponding row
c        number in i4
c
                IF ( ABS(a(i,k)).GT.amax ) THEN
                  amax = ABS(a(i,k))
                  i4 = i
                END IF
              END IF
            END DO
          END DO
        END DO
c
c        the elimination is now complete and the a-matrix has been
c        reduced to the product of an upper and lower triangle matrix
c
c
c        replace the a-matrix by its inverse
c
c        first invert the lower triangle matrix and store on itself
c
        DO i1 = 1 , mm
          i = m + 1 - i1
          i2 = i - 1
          DO j1 = 1 , i2
            j = i2 + 1 - j1
            j2 = j + 1
            w1 = -a(i,j)
            IF ( i2.GE.j2 ) THEN
              DO k = j2 , i2
                w1 = w1 - a(k,j)*c(k)
              END DO
            END IF
            c(j) = w1
          END DO
          DO k = 1 , i2
            a(i,k) = c(k)
          END DO
        END DO
c
c        now invert the upper triangle matrix and form the original
c        a-matrix apart from column interchange.this overwrites the
c        original a-matrix.
c
        DO i1 = 1 , m
          i = m + 1 - i1
          i2 = i + 1
          w = 1.0/a(i,i)
          DO j = 1 , m
            IF ( i.LT.j ) THEN
              w1 = 0.0
            ELSE IF ( i.EQ.j ) THEN
              w1 = 1.0
            ELSE
              w1 = a(i,j)
            END IF
            IF ( i1.GT.1 ) THEN
              DO k = i2 , m
                w1 = w1 - a(i,k)*a(k,j)
              END DO
            END IF
            c(j) = w1
          END DO
          DO j = 1 , m
            a(i,j) = c(j)*w
          END DO
c
c        re-order the columns of the inverse a-matrix to coincide with
c        the order of the rows of the a-matrix on input
c
        END DO
        DO i = 1 , m
 20       IF ( ind(i).NE.i ) THEN
            j = ind(i)
            DO k = 1 , m
              sto = a(k,i)
              a(k,i) = a(k,j)
              a(k,j) = sto
            END DO
            isto = ind(j)
            ind(j) = j
            ind(i) = isto
            GO TO 20
          END IF
        END DO
      END IF
      RETURN
      END


      FUNCTION WIG6J(ja,jb,jc,jd,je,jf)
      IMPLICIT REAL*4 (A-H, O-Z)
 
      INCLUDE 'par.h'

c
c        evaluation of wigner 6-j symbol
c
c                           ( ja jb je )
c                           ( jd jc jf )
c
      COMMON /GAMDAT/ gamma(NGAM) , gamaf(NGAM) , gamahf(NGAM)
c
c        gamma contains logarithms of factorials. same dimension as
c        in the main.
c
      lc = ja + jb + jc + jd
      xlc = (-1)**lc
      IF ( ja.EQ.0 ) THEN
        a = (jb+jb+1)*(jc+jc+1)
        WIG6J = (-1)**(jb+jc+jd)/SQRT(a)
        WIG6J = WIG6J*xlc
      ELSE IF ( jf.EQ.0 ) THEN
        a = (ja+ja+1)*(jb+jb+1)
        WIG6J = (-1)**(ja+jb+je)/SQRT(a)
        WIG6J = WIG6J*xlc
      ELSE IF ( je.EQ.0 ) THEN
        a = (ja+ja+1)*(jd+jd+1)
        WIG6J = (-1)**(ja+jd+jf)/SQRT(a)
        WIG6J = WIG6J*xlc
      ELSE
        i1 = ja + jb + je + 2
        i2 = ja + jb - je + 1
        i3 = ja + je - jb + 1
        i4 = je + jb - ja + 1
        d1 = gamma(i2) + gamma(i3) + gamma(i4) - gamma(i1)
        i1 = ja + jc + jf + 2
        i2 = ja + jc - jf + 1
        i3 = ja + jf - jc + 1
        i4 = jf + jc - ja + 1
        d2 = gamma(i2) + gamma(i3) + gamma(i4) - gamma(i1)
        i1 = jb + jd + jf + 2
        i2 = jb + jd - jf + 1
        i3 = jb + jf - jd + 1
        i4 = jf + jd - jb + 1
        d3 = gamma(i2) + gamma(i3) + gamma(i4) - gamma(i1)
        i1 = jc + jd + je + 2
        i2 = jc + jd - je + 1
        i3 = jc + je - jd + 1
        i4 = je + jd - jc + 1
        d4 = gamma(i2) + gamma(i3) + gamma(i4) - gamma(i1)
        delt = (d1+d2+d3+d4)*0.5E0
        k1 = ja + jb + je
        k2 = ja + jc + jf
        k3 = jb + jd + jf
        k4 = jc + jd + je
        k5 = ja + jb + jc + jd
        k6 = jb + jc + je + jf
        k7 = ja + jd + je + jf
        lmini = MAX0(0,k1,k2,k3,k4)
        lmaxi = MIN0(k5,k6,k7)
        l = lmini + 1
        lm = lmini - 1
        l1 = l + 1
        l2 = l - k1
        l3 = l - k2
        l4 = l - k3
        l5 = l - k4
        l6 = k5 - lm
        l7 = k6 - lm
        l8 = k7 - lm
        delt = delt + gamma(l1) - gamma(l2) - gamma(l3) - gamma(l4)
     &         - gamma(l5) - gamma(l6) - gamma(l7) - gamma(l8)
        a = 1.E0
        IF ( lmaxi.NE.lmini ) THEN
          jj = lmaxi
          DO j = l , lmaxi
            jk = jj - 1
            c = (jj-k1)*(jj-k2)
            a = 1.E0 - a*(jj+1)*(k5-jk)*(k6-jk)*(k7-jk)
     &          /(c*(jj-k3)*(jj-k4))
            jj = jk
          END DO
        END IF
        WIG6J = a*EXP(delt)*(-1)**lmini
        WIG6J = WIG6J*xlc
      END IF
      RETURN
      END


      FUNCTION WIG3J(j1,j2,j3,m1,m2,mf)
      IMPLICIT REAL*4 (A-H, O-Z)
 
      INCLUDE 'par.h'

c
c        evaluation of wigner 3-j symbol
c
c                            ( j1 j2 j3 )
c                            ( m1 m2 mf )
c
      COMMON /GAMDAT/ gamma(NGAM) , gamaf(NGAM) , gamahf(NGAM)
c
c        gamma contains logarithms of factorials. same dimension as
c        in the main.
c
      IF ( j3.EQ.0 ) THEN
        a = j1 + j1 + 1
        WIG3J = (-1)**(j1-m1)/SQRT(a)
      ELSE
        m3 = -mf
        i1 = j1 + j2 - j3 + 1
        i2 = j1 - m1 + 1
        i3 = j2 - m2 + 1
        i4 = j3 - m3 + 1
        i5 = j3 + m3 + 1
        i6 = j1 + j2 + j3 + 2
        i7 = j1 - j2 + j3 + 1
        i8 = j2 - j1 + j3 + 1
        i9 = j1 + m1 + 1
        i10 = j2 + m2 + 1
        x = (gamma(i1)+gamma(i2)+gamma(i3)+gamma(i4)+gamma(i5)-gamma(i6)
     &      -gamma(i7)-gamma(i8)-gamma(i9)-gamma(i10))*0.50
        k1 = j1 + m1
        k2 = j2 + j3 - m1
        k3 = j1 - m1
        k4 = j3 - m3
        k5 = j3 - m1 - j2
        lmini = MAX0(0,k5)
        lmaxi = MIN0(k2,k3,k4)
        l = lmini + 1
        lm = lmini - 1
        l1 = l + k1
        l2 = k2 - lm
        l3 = k3 - lm
        l4 = k4 - lm
        l5 = l - k5
        x = x + gamma(l1) + gamma(l2) - gamma(l) - gamma(l3) - gamma(l4)
     &      - gamma(l5)
        a = 1.E0
        IF ( lmaxi.NE.lmini ) THEN
          jj = lmaxi
          DO j = l , lmaxi
            jk = jj - 1
            c = jj*(jj-k5)
            a = 1.E0 - a*(jj+k1)*(k3-jk)*(k4-jk)/(c*(k2-jk))
            jj = jk
          END DO
        END IF
        WIG3J = a*EXP(x)*(-1)**(i10+lm)
      END IF
      RETURN
      END


      FUNCTION WIG3J0(j1,j2,j3)
      IMPLICIT REAL*4 (A-H, O-Z)
 
      INCLUDE 'par.h'

c
c        evaluation of wigner 3-j symbols of the form
c
c                            ( j1 j2 j3 )
c                            (  0  0  0 )
c
c **********************************************************************
c
      COMMON /GAMDAT/ gamma(NGAM) , gamaf(NGAM) , gamahf(NGAM)
c
c        gamma contains logarithms of factorials. same dimension as
c        in the main.
c
      IF ( j1.EQ.0 ) THEN
        a = j3 + j3 + 1
        WIG3J0 = (-1)**j3/SQRT(a)
      ELSE IF ( j3.EQ.0 ) THEN
        a = j1 + j1 + 1
        WIG3J0 = (-1)**j1/SQRT(a)
      ELSE
        i1 = (j1+j2+j3)/2
        fatt = (-1)**i1
        i1 = i1 + 1
        i2 = j1 + j2 - j3 + 1
        i3 = j1 - j2 + j3 + 1
        i4 = -j1 + j2 + j3 + 1
        i5 = j1 + j2 + j3 + 2
        xnum = (gamma(i2)+gamma(i3)+gamma(i4)-gamma(i5))*0.5E0
        i2 = i1 - j3
        i3 = i1 - j2
        i4 = i1 - j1
        xnam = gamma(i1) - gamma(i2) - gamma(i3) - gamma(i4)
        WIG3J0 = fatt*EXP(xnum+xnam)
      END IF
      RETURN
      END


      SUBROUTINE FACTOR
      IMPLICIT REAL*4 (A-H, O-Z)
 
      INCLUDE 'par.h'

c
c          this subroutine calculates the following factorials
c
c              gamma(n)        log of gamma(n)
c              gamaf(n)        gamma(n)
c              gamahf(n)       gamma(n+0.5)
c
      DIMENSION xid(400)
      COMMON /GAMDAT/ gamma(NGAM) , gamaf(NGAM) , gamahf(NGAM)

      gamma(1) = 0.0
      gamma(2) = 0.0
      DO i = 3 , 400
        ix = i - 1
        xi = ix
        gamma(i) = gamma(ix) + ALOG(xi)
      END DO
      gamaf(1) = 1.0
      gamaf(2) = 1.0
      DO i = 2 , 34
        gamaf(i+1) = gamaf(i)*i
      END DO
      gamahf(1) = 1.7724538510
      gamahf(2) = 0.8862269255
      xid(1) = 1.0
      DO i = 2 , 21
        xi = i
        xid(i) = (2.0*xi-1.0)*xid(i-1)/2.0**i
        gamahf(i+1) = xid(i)*gamahf(1)
      END DO
      RETURN
      END


      FUNCTION IPARITY(l)
      IMPLICIT REAL*4 (A-H, O-Z)
 
      k = l/2
      kk = 2*k - l
      IF ( kk.EQ.0 ) THEN
        IPARITY = 1
      ELSE
        IPARITY = -1
      END IF
      RETURN
      END

