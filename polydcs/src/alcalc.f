      SUBROUTINE ALC
      IMPLICIT REAL*4 (A-H, O-Z)
 
c
c     To calculate the A (J-->J') elements for spherical rotors
c                       L
c
 
      INCLUDE 'par.h'
 
      COMMON /ALCOEF/ albcc(2,NLAM) , altcc(2,NLAM,11) , al(NLAM)
      COMMON /INPPAR/ ipol , iprint , nt , nscat , mol ,
     &                lmax , lmaxb , lbig ,
     &                alpha0 , alpha2 , dpol , dpolex , quad , quadex ,
     &                dpolt , dpole , q1 , q2 , oct1 , oct2
      COMMON /JKTDAT/ jmin, jmax , mnlamda , mxlamda ,
     &                jint , jfnl , itai , itaf , numti , numtf
      COMMON /MMAT  / amr(0:MAXN) , ami(0:MAXN) ,
     &                amrb(0:MAXN) , amib(0:MAXN)

      DO ll = 1 , lbig
        al(ll) = 0.0
      END DO
 
      DO ll = 1 , lbig
        ll1 = ll - 1
        DO l = 1 , lmax
          l1 = l - 1
          lx = l1*lmax*jmax
          xl1 = FLOAT(l+l-1)
          DO lb = 1 , lmax
            l2 = lb - 1
            lx1 = l2*lmax*jmax
            xl2 = FLOAT(lb+lb-1)
            lmin = IABS(l1-l2)
            lmx = l1 + l2
            IF ( ll1.GE.lmin .AND. ll1.LE.lmx ) THEN
              lz = l1 + l2 + ll1
              lmz = lz/2
              IF ( lmz*2.EQ.lz ) THEN
                llz = lmz + l2
                w1 = WIG3J0(l2,l1,ll1)
                x1 = (-1)**llz*w1*SQRT(xl1*xl2)
                DO lp = 1 , lmax
                  l3 = lp - 1
                  xl3 = FLOAT(lp+lp-1)
                  lxp = l3*jmax
                  lmlp = IABS(l1-l3)
                  lplp = l1 + l3
                  DO lpb = 1 , lmax
                    l4 = lpb - 1
                    lx1p = l4*jmax
                    lmin1 = IABS(l3-l4)
                    lmx1 = l3 + l4
                    IF ( ll1.GE.lmin1 .AND. ll1.LE.lmx1 ) THEN
                      lz1 = l3 + l4 + ll1
                      lmz1 = lz1/2
                      IF ( lmz1*2.EQ.lz1 ) THEN
                        lbmlpb = IABS(l2-l4)
                        lbplpb = l2 + l4
                        xl4 = FLOAT(l4+l4+1)
                        lzz = lmz1 + l3
                        w2 = WIG3J0(l3,l4,ll1)
                        x2 = (-1)**lzz*w2*SQRT(xl3*xl4)*x1
                        DO j = jmin , jmax
                          jj = j - 1
                          IF ( jj.GE.lmlp .AND. jj.LE.lplp ) THEN
                            IF ( jj.GE.lbmlpb .AND. jj.LE.lbplpb ) THEN
                              x3 = (-1)**jj*x2*WIG6J(l1,l3,l2,l4,jj,ll1)
                              jjx = jj + jj + 1
                              DO m = 1 , jjx
c                               mj = j - m
                                n1 = (lx+lxp+jj)*jjx + m
                                n2 = (lx1+lx1p+jj)*jjx + m
                                IF ( l.EQ.lb .AND. lp.EQ.lpb ) THEN
                                  al(ll) = al(ll)
     &                              + (amr(n1)*amr(n1)+ami(n1)*ami(n1))
     &                              *x3
                                ELSE
                                  al(ll) = al(ll)
     &                              + (amr(n1)*amr(n2)+ami(n1)*ami(n2))
     &                              *x3
                                END IF
                              END DO
                            END IF
                          END IF
                        END DO
                      END IF
                    END IF
                  END DO
                END DO
              END IF
            END IF
          END DO
        END DO
      END DO
      RETURN
      END
 
 
      SUBROUTINE ALI
      IMPLICIT REAL*4 (A-H, O-Z)
 
c
c     To calculate the A (J,K (Tau)-->J',K' (Tau')) elements for 
c                       L
c     symmetric and asymmetric rotors
c
 
      INCLUDE 'par.h'
 
      DIMENSION alb(NLAM)

      COMMON /ALCOEF/ albcc(2,NLAM) , altcc(2,NLAM,11) , al(NLAM)
      COMMON /INPPAR/ ipol , iprint , nt , nscat , mol ,
     &                lmax , lmaxb , lbig ,
     &                alpha0 , alpha2 , dpol , dpolex , quad , quadex ,
     &                dpolt , dpole , q1 , q2 , oct1 , oct2
      COMMON /JKTDAT/ jmin, jmax , mnlamda , mxlamda ,
     &                jint , jfnl , itai , itaf , numti , numtf
      COMMON /MMAT  / amr(0:MAXN) , ami(0:MAXN) ,
     &                amrb(0:MAXN) , amib(0:MAXN)
 
      DO ll = 1 , lbig
        al(ll) = 0.0E0
        alb(ll) = 0.0E0
      END DO

      llm = lmax
      IF ( ipol.NE.0 ) llm = lmaxb
  
c
c     do loop over L
c
      DO ll = 1 , lbig
        ll1 = ll - 1
c
c     do loop over l
c
        DO l = 1 , llm
          l1 = l - 1
          lx = l1*llm*jmax
          xl1 = FLOAT(l+l-1)
c                  _
c     do loop over l
c
          DO lb = 1 , llm
            l2 = lb - 1
            lx1 = l2*llm*jmax
            xl2 = FLOAT(lb+lb-1)
            lmin = IABS(l1-l2)
            lmx = l1 + l2
            IF ( ll1.GE.lmin .AND. ll1.LE.lmx ) THEN
              lz = l1 + l2 + ll1
              lmz = lz/2
              IF ( IPARITY(lz).EQ.1 ) THEN
                llz = lmz + l2
                w1 = WIG3J0(l2,l1,ll1)
                x1 = (-1)**llz*w1*SQRT(xl1*xl2)
c
c     do loop over l'
c
                DO lp = 1 , lmaxb
                  l3 = lp - 1
                  xl3 = FLOAT(lp+lp-1)
                  lxp = l3*jmax
                  lmlp = IABS(l1-l3)
                  lplp = l1 + l3
c                  _
c     do loop over l'
c
                  DO lpb = 1 , lmaxb
                    l4 = lpb - 1
                    lx1p = l4*jmax
                    lmin1 = IABS(l3-l4)
                    lmx1 = l3 + l4
                    IF ( ll1.GE.lmin1 .AND. ll1.LE.lmx1 ) THEN
                      lz1 = l3 + l4 + ll1
                      lmz1 = lz1/2
                      IF ( IPARITY(lz1).EQ.1 ) THEN
                        lbmlpb = IABS(l2-l4)
                        lbplpb = l2 + l4
                        xl4 = FLOAT(l4+l4+1)
                        lzz = lmz1 + l3
                        w2 = WIG3J0(l3,l4,ll1)
                        x2 = (-1)**lzz*w2*SQRT(xl3*xl4)*x1
c
c     do loop over j
c
                        DO j = jmin , jmax
                          jj = j - 1
                          tj = FLOAT(jj+jj+1)
                          IF ( jj.GE.lmlp .AND. jj.LE.lplp ) THEN
                            IF ( jj.GE.lbmlpb .AND. jj.LE.lbplpb ) THEN
                              x3 = (-1)
     &                             **jj*tj*x2*WIG6J(l1,l3,l2,l4,jj,ll1)
                              jjx = jj + jj + 1
c
c     do loop over mj'
c
                              DO m = 1 , jjx
                                mj = j - m
                                n1 = (lx+lxp+jj)*jjx + m
                                IF ( mol.NE.1 .OR. mj.EQ.(itaf-itai) )
     &                               THEN
c                  _
c     do loop over mj'
c
                                  DO md = 1 , jjx
                                    mjd = j - md
                                    IF ( mol.NE.1 .OR. 
     &                                mjd.EQ.(itaf-itai) ) THEN
                                      n2 = (lx1+lx1p+jj)*jjx + md
                                      al(ll) = al(ll)
     &                                  + (amr(n1)*amr(n2)+ami(n1)
     &                                  *ami(n2))*x3
                                      albcc(1,ll) = al(ll)
                                      IF ( ipol.NE.0 ) THEN
                                        alb(ll) = alb(ll)
     &                                    + (amrb(n1)*amrb(n2)+amib(n1)
     &                                    *amib(n2))*x3
                                        albcc(2,ll) = alb(ll)
 
                                      END IF
                                    END IF
                                  END DO
                                END IF
                              END DO
                            END IF
                          END IF
                        END DO
                      END IF
                    END IF
                  END DO
                END DO
              END IF
            END IF
          END DO
        END DO
      END DO
      RETURN
      END


      SUBROUTINE ALL
      IMPLICIT REAL*4 (A-H, O-Z)
 
c
c     To calculate the A (J-->J') elements for linear molecules
c                       L
c

      INCLUDE 'par.h'

      COMMON /ALCOEF/ albcc(2,NLAM) , altcc(2,NLAM,11) , al(NLAM)
      COMMON /INPPAR/ ipol , iprint , nt , nscat , mol ,
     &                lmax , lmaxb , lbig ,
     &                alpha0 , alpha2 , dpol , dpolex , quad , quadex ,
     &                dpolt , dpole , q1 , q2 , oct1 , oct2
      COMMON /JKTDAT/ jmin, jmax , mnlamda , mxlamda ,
     &                jint , jfnl , itai , itaf , numti , numtf
      COMMON /MMATL / ar(11,NPS,NPS) , aii(11,NPS,NPS) ,
     &                arb(11,NPS,NPS) , aib(11,NPS,NPS)

 
      DO ll = 1 , lbig
        al(ll) = 0.0
      END DO
      llmmax = lmax
      IF ( mol.EQ.4 .AND. ipol.NE.0 ) llmmax = 3*lbig - 2
c
c     do loop over L
c
      DO ll = 1 , lbig
        ll1 = ll - 1
        xll1 = FLOAT(ll1+ll1+1)
c
c     do loop over l
c
        DO l = 1 , llmmax
          l1 = l - 1
          xl1 = FLOAT(l+l-1)
c                  _
c     do loop over l
c
          DO lb = 1 , llmmax
            l2 = lb - 1
            xl2 = FLOAT(lb+lb-1)
            lmin = IABS(l1-l2)
            lmx = l1 + l2
            IF ( ll1.GE.lmin .AND. ll1.LE.lmx ) THEN
              lz = l1 + l2 + ll1
              IF ( IPARITY(lz).EQ.1 ) THEN
                w1 = WIG3J0(l1,l2,ll1)
                x1 = w1*SQRT(xl1*xl2)
c
c     do loop over l'
c
                DO lp = 1 , llmmax
                  l3 = lp - 1
c
c     For linear molecules with a center of symmetry only those
c     T-matrix elements with the same parity are different from zero
c
                  IF ( mol.NE.5 .OR. IPARITY(l1).EQ.IPARITY(l3) ) THEN
                    xl3 = FLOAT(lp+lp-1)
                    lmlp = IABS(l1-l3)
                    lplp = l1 + l3
c                  _
c     do loop over l'
c
                    DO lpb = 1 , llmmax
                      l4 = lpb - 1
                      IF ( mol.NE.5 .OR. IPARITY(l2).EQ.IPARITY(l4) )
     &                     THEN
                        lmin1 = IABS(l3-l4)
                        lmx1 = l3 + l4
                        IF ( ll1.GE.lmin1 .AND. ll1.LE.lmx1 ) THEN
                          lz1 = l3 + l4 + ll1
                          IF ( IPARITY(lz1).EQ.1 ) THEN
                            lbmlpb = IABS(l2-l4)
                            lbplpb = l2 + l4
                            xl4 = FLOAT(l4+l4+1)
                            w2 = WIG3J0(l3,l4,ll1)
                            x2 = w2*SQRT(xl3*xl4)
     &                           *x1*IPARITY((l1+l2+l3+l4)/2+ll1+l2+l3)
                            hjm = MAX(l,lb,lp,lpb)
c
c     do loop over j
c
                            DO j = jmin , jmax
                              jj = j - 1
                              hji = 0
                              IF ( hjm.LT.2*lbig .OR. jj.LT.2 ) THEN
 
                                IF ( .NOT.(hjm.GT.2*lbig+ll1 .AND. (jj
     &                               .EQ.1 .OR. jj.EQ.2 .OR. jj.EQ.0)) )
     &                               THEN
                                  IF ( jj.GT.2 ) hji = 1
                                  IF ( hji.NE.1 .OR. 
     &                                 hjm.LE.2*lbig-1 .OR. ipol.EQ.0 )
     &                                 THEN
                                    IF ( jj.GE.lmlp .AND. jj.LE.lplp )
     &                                THEN
                                      IF ( jj.GE.lbmlpb .AND. 
     &                                  jj.LE.lbplpb ) THEN
                                        x3a = 0.0
                                        IF ( ll1.GE.IABS(l1-l2) .AND. 
     &                                    ll1.LE.(l1+l2) .AND. 
     &                                    jj.GE.IABS(l1-l3) .AND. 
     &                                    jj.LE.(l1+l3) .AND. 
     &                                    jj.GE.IABS(l2-l4) .AND. 
     &                                    jj.LE.(l2+l4) .AND. 
     &                                    ll1.GE.IABS(l3-l4) .AND. 
     &                                    ll1.LE.(l3+l4) )
     &                                    x3a = WIG6J(l1,l2,l3,l4,ll1,
     &                                    jj)
 
                                        x3 = (-1)**jj*x2*x3a
c                                       jjx = jj + jj + 1
 
                                        IF 
     &                                    ( (ar(j,l,lp).NE.0.0 .AND. ar(
     &                                    j,lb,lpb).NE.0.0) .OR. 
     &                                    (aii(j,l,lp).NE.0.0 .AND. 
     &                                    aii(j,lb,lpb).NE.0.0) )
     &                                    altcc(1,ll,j) = altcc(1,ll,j)
     &                                    + 
     &                                    (ar(j,l,lp)*ar(j,lb,lpb)+aii(j
     &                                    ,l,lp)*aii(j,lb,lpb))*x3*xll1
                                        IF ( ipol.NE.0 ) THEN
                                         altcc(2,ll,j) = altcc(2,ll,j)
     &                                     + (arb(j,l,lp)*arb(j,lb,lpb)
     &                                     +aib(j,l,lp)*aib(j,lb,lpb))
     &                                     *x3*xll1
                                        END IF
                                      END IF
                                    END IF
                                  END IF
                                END IF
                              END IF
                            END DO
                          END IF
                        END IF
                      END IF
                    END DO
                  END IF
                END DO
              END IF
            END IF
          END DO
        END DO
      END DO

      RETURN
      END
