CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     *****  P O L Y D C S  *****                                      C
C                                                                      C
C     N. Sanna and F. A. Gianturco                                     C
C                                                                      C
C     Comp. Phys. Comm., ...                                           C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM POLYDCS
      IMPLICIT REAL*4 (A-H, O-Z)
 
      INCLUDE 'par.h'
 
      CHARACTER*6  stat
      CHARACTER*80 title , blmfn , kmatfn

      DIMENSION dcsb(NANG) , dcs(NANG)
      DIMENSION sumdcs(2,NANG) , totdcs(2,NANG)
      DIMENSION ji(NDCS) , jf(NDCS)
      DIMENSION sumalbcc(2,NLAM)
      DIMENSION ntaui(11) , ntauf(11) , itaui(11) , itauf(11)

      COMMON /ALCOEF/ albcc(2,NLAM) , altcc(2,NLAM,11) , al(NLAM)
      COMMON /ANGDAT/ th(NANG) , pl(NLAM,NANG)
      COMMON /ASYDAT/ ai(11,6,2) , af(11,6,2) , egi(11) , egf(11)
      COMMON /BLMCOF/ br(NSYM,NPW,NPW) , bi(NSYM,NPW,NPW)
      COMMON /GAMDAT/ gamma(NGAM) , gamaf(NGAM) , gamahf(NGAM)
      COMMON /INPCHR/ stat(NSYM) , blmfn , kmatfn
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
      COMMON /TOTDCS/ cst1 , csm1 , cst2 , csm2 , cstb , csmb

      DATA totdcs/182*0.0E0/
      DATA totcsmb , totcsm , totcsm1 , totcstb , totcst , 
     &     totcst1/6*0.0E0/
      DATA dcsb , roteff/NANG*0.0 , 0.0/

c
c     Read first general input data from stdin
c
      READ (5,99001) title
      READ (5,99001) blmfn
      READ (5,99001) kmatfn
 
      READ (5,*) mol , ipol , ni
 
      IF ( ni.EQ.1 .AND. mol.EQ.4 ) READ (5,*) xom , n0 , n1
 
      IF ( mol.EQ.0 .OR. mol.GT.3 ) READ (5,*) bz
      IF ( mol.EQ.1 ) READ (5,*) az , bz
      IF ( mol.EQ.2 ) READ (5,*) az , bz , cz

      IF ( ipol.NE.0 ) THEN
        IF ( mol.EQ.2 .OR. mol.EQ.1 ) THEN
          READ (5,*) dpolt , dpole , q1 , q2 , oct1 , oct2
        ELSE IF ( mol.EQ.4 ) THEN 
          READ (5,*) alpha0 , alpha2 , dpol , dpolex , quad , quadex
        END IF
      END IF
 
      READ (5,*) energy

      READ (5,*) nscat , lmax , lmaxb , lbig , nt , iprint

      READ (5,*) (ji(i),jf(i),i=1,nt)

      DO i = 1 , nscat
        READ (5,99003) stat(i) , ksiz(i) , ksizb(i) , mult(i)
      END DO
c
c     Write input data on stdout
c
      WRITE (6,99000)
      WRITE (6,99001) title
      WRITE (6,99008) nscat , lmax , lmaxb , lbig , nt , iprint
      IF ( ipol.NE.0 ) WRITE (6,99016) ipol
 
      IF ( iprint.GE.1 ) THEN

        IF ( ni.EQ.1 .AND. mol.EQ.4 ) WRITE (6,99009) n0 , n1 , xom

        WRITE (6,99007) az , bz , cz

        IF ( ipol.NE.0 ) THEN
          IF ( mol.EQ.2 .OR. mol.EQ.1 ) 
     &      WRITE (6,99017) dpolt , dpole , q1 , q2 , oct1 , oct2
          IF ( mol.EQ.4 ) 
     &      WRITE (6,99018) alpha0 , alpha2 , 
     &                      dpol , dpolex , quad , quadex
        END IF
 
        WRITE (6,99004)
        DO i = 1 , nscat
          WRITE (6,99005) i , stat(i) , ksiz(i) , ksizb(i) , mult(i)
        END DO

      END IF
 
      WRITE (6,99010)
      DO i = 1 , nt
        WRITE (6,99011) i , ji(i) , jf(i)
      END DO

c
c     Convert energy (en) in Rydberg (1/2 a.u.) and print it.
c
      en = energy/13.606E0
      en4 = 4.0E0*en
      xki = SQRT(en)
 
      IF ( ni.EQ.1 .AND. mol.EQ.4 ) THEN
        xxom = xom/13.606
        en4 = en4 - 4*xxom*(FLOAT(n0)+0.5)
        xki = SQRT(en-xxom*(FLOAT(n0)+0.5))
      END IF

      WRITE (6,99014) energy , en , xki

      IF ( mol.EQ.0 .OR. mol.GT.3 ) brot = bz
      lmax = lmax + 1
      lmaxb = lmaxb + 1
      lbig = lbig + 1
 
c
c     Check if input data are within the thresholds defined
c     in the parameters list
c
      IF ( nscat.GT.NSYM ) THEN
        WRITE (6,99025) nscat , NSYM
        STOP
      END IF
      IF ( lbig.GT.NLAM ) THEN
        WRITE (6,99026) lbig , NLAM
        STOP
      END IF
      IF ( nt.GT.NDCS ) THEN
        WRITE (6,99027) nt , NDCS
        STOP
      END IF
      DO i = 1 , nscat
        IF ( ksiz(i).GT.NPW ) THEN
          WRITE (6,99028) ksiz(i) , NPW , i
          STOP
        END IF
      END DO
      DO i = 1 , nscat
        IF ( ksizb(i).GT.NPW ) THEN
          WRITE (6,99029) ksizb(i) , NPW , i
          STOP
        END IF
      END DO
c
c    Call blmread to read the blm data for spherical,
c    symmetric and asymmetric top molecules
c
      IF ( mol.LE.2 ) CALL BLMREAD
c
c     Call factor to generate the factorial
c
      CALL FACTOR
c
c     Call kread to read the K-matrices and calculate T-matrices.
c
      CALL KREAD(en)
c
c     For all but linear molecules, check if lmax/lmaxb are less
c     then the maximum values of lmax/lmaxb read in from k-mat/blm files.
c     If this is the case, issue a warning message and set
c     ksiz(), ksizb() accordingly.
c
      IF ( mol.LT.3 ) THEN
        DO i=1, nscat
          newsiz  = 0
          ks  = ksiz(i)
          DO j=1,ks
            IF ( (lk(i,j).GT.lmax-1) ) THEN 
              newsiz = j - 1
              GO TO 10
            END IF
          END DO
 10       IF ( newsiz.NE.0 ) THEN
            ksiz(i) = newsiz
            IF ( i.EQ.1 ) WRITE (6,99039) lmax-1 , lk(i,ks)
            WRITE (6,99040) i , ksiz(i) , ksiz(i)
          END IF 
        END DO

        IF ( ipol.NE.0 ) THEN
          DO i=1, nscat
            newsizb = 0
            ksb = ksizb(i)
            DO j=1,ksb
              IF ( (lk(i,j).GT.lmaxb-1) ) THEN 
                newsizb = j - 1
                GO TO 11
              END IF
            END DO
 11         IF ( newsizb.NE.0 ) THEN
              ksizb(i) = newsizb
              IF ( i.EQ.1 ) WRITE (6,99041) lmaxb-1 , lk(i,ksb)
              WRITE (6,99042) i , ksizb(i) , ksizb(i)
            END IF
          END DO
        END IF 
      END IF
c
c     Print out blms data for spherical, symmetric and 
c     asymmetric top molecules and k-matrices data for all
c     molecular types.
c
      IF ( iprint.GT.1 ) CALL BKPRINT

C--------------------------------------------------------------
C     Here begin the DCS calculation
C--------------------------------------------------------------
 
c
c     To generate the scattering angles
c
      step = 180.0E0/FLOAT(NANG-1)
      th(1) = 0.0E0
      DO i = 2 , NANG
        th(i) = th(i-1) + step
      END DO
c
c     To call subroutine pleg and generate Legendre Polynomials
c     pl(l,theta)
c
      CALL PLEG
c
c     For polar molecules (IPOL=1), calculates the Born K-matrices.
c
      IF ( ipol.NE.0 ) CALL BORNTM(xki)
 
c
c     To calculate cross-section for elastic and inelastic case
c
      DO ntr = 1 , nt
        DO i = 1 , NANG
          DO k = 1 , 2
            sumdcs(k,i) = 0.0E0
          END DO
        END DO
        sumcsmb = 0.0E0
        sumcstb = 0.0E0
        sumcst1 = 0.0E0
        sumcsm1 = 0.0E0
        sumcsm = 0.0E0
        sumcst = 0.0E0
        jint = ji(ntr)
        jint1 = jint + 1
        jfnl = jf(ntr)
        jfnl1 = jfnl + 1
        j21 = jint + jint + 1
        jp21 = jfnl + jfnl + 1
        xjp21 = FLOAT(jp21)
        fact = xjp21/en4
        xjk1 = FLOAT(jint*(jint+1))
        xjk2 = FLOAT(jfnl*(jfnl+1))
        jmfl = jint - jfnl

        jmin = IABS(jmfl) + 1
        jmax = jint + jfnl + 1

c
c***** Begin code for SPHERICAL TOP molecule (MOL=0) *****
c
        IF ( mol.EQ.0 ) THEN
          rote = brot*(xjk2-xjk1)/13.606E0
          IF ( en.LE.rote ) GO TO 100
          xkf = SQRT(en-rote)
c
c     Calculate and print CS from T-matrices
c
          IF ( ntr.EQ.1 ) CALL POLCROSS(en4)
c
c     Call mmatrx to generate the M-matrix
c
          CALL MMATRX
c
c     Call alc (spherical top) to generate the Al coefficients
c     of the DCS's expansion
c
          CALL ALC
 
          DO ll = 1 , lbig
            ll1 = ll - 1
            xl = FLOAT(ll1+ll1+1)
            albcc(1,ll) = al(ll)*fact*xl
          END DO
c
c     Call cross to calculate and print Al coefficients, DCS, CST and CSM
c
          WRITE (6,*)
          WRITE (6,99002)
          WRITE (6,99006) jint , jfnl
          WRITE (6,99002)

          CALL CROSS(albcc,en,xki,xkf,dcsb,dcs,csm,cst)
 
          DO ii = 1 , NANG
            totdcs(1,ii) = totdcs(1,ii) + dcs(ii)
          END DO
 
          WRITE (6,99022) csm , cst
 
          totcsm = totcsm + csm
          totcst = totcst + cst
 
          roteff = roteff + cst*rote/2

        END IF
c 
c***** End code for SPHERICAL TOP molecule (MOL=0) *****
c
c
c***** Begin code for LINEAR molecule (MOL=4,5) *****
c
        IF ( mol.GT.3 ) THEN
          rote = brot*(xjk2-xjk1)/13.606E0
          IF ( en.LE.rote ) GO TO 100
          xkf = SQRT(en-rote)
 
          IF ( ni.EQ.1 .AND. mol.EQ.4 )
     &         xkf = SQRT(en-rote-xxom*(FLOAT(n1)+0.5))
 
          IF ( ntr.EQ.1 ) THEN
c
c     Calculate and print CS from T-matrices
c
            CALL LINCROSS(en4)
c 
c     Call mmatrx to generate the M-matrix
c
            CALL MMATRX
c
c     Call all (linear molecules) only once to generate
c     the Al coefficients of the DCS's expansion
c
            CALL ALL

          END IF
 
          DO ll = 1 , lbig
            ll1 = ll - 1
            DO lt = jmin , jmax
              IF ( lt.EQ.jmin ) albcc(1,ll) = 0.0
              xx1 = WIG3J0(jint,lt-1,jfnl)
              xx1 = xx1**2
              albcc(1,ll) = albcc(1,ll) + altcc(1,ll,lt)*xx1*fact
              IF ( ipol.NE.0 .AND. lt.LE.3 ) THEN
                IF ( lt.EQ.jmin ) albcc(2,ll) = 0.0
                albcc(2,ll) = albcc(2,ll) + altcc(2,ll,lt)*fact*xx1
              END IF
            END DO
            sumalbcc(1,ll) = sumalbcc(1,ll) + albcc(1,ll)
          END DO

          WRITE (6,*) 
          WRITE (6,99002)
          WRITE (6,99006) jint , jfnl
          WRITE (6,99002)
c
c     Call cross for generate and print DCS, CST and CSM
c
          CALL CROSS(albcc,en,xki,xkf,dcsb,dcs,csm,cst)
 
          roteff = roteff + cst*rote/2
 
          DO ii = 1 , NANG
            totdcs(1,ii) = totdcs(1,ii) + dcs(ii)
            IF ( ipol.NE.0 ) totdcs(2,ii) = totdcs(2,ii) + dcsb(ii)
          END DO

          WRITE (6,99022) csm , cst
 
          totcst1 = totcst1 + cst1
          totcsm1 = totcsm1 + csm1
 
          totcstb = totcstb + cstb
          totcsmb = totcsmb + csmb
 
          totcsm = totcsm + csm
          totcst = totcst + cst
 
        END IF
c 
c***** End code for LINEAR molecule (MOL=4,5) *****
c
c
c***** Begin code for SYMMETRIC/ASYMMETRIC TOP molecule (MOL=1,2) *****
c
        IF ( mol.EQ.1 .OR. mol.EQ.2 ) THEN
c
c     Calculate and print CS from T-matrices
c
          IF ( ntr.EQ.1 ) CALL POLCROSS(en4)
c
c     For asymmetric top molecules, call "asymtop" to have
c     the energy levels and eigenfunction expansion coefficients.
c     For other molecular symmetries, the rotational eigenfunction
c     expansion coefficients and energies are trivial.
c
          IF ( mol.EQ.2 ) THEN

            DO ii1 = 1 , 11
              DO ii2 = 1 , 6
                DO ii3 = 1 , 2
                  ai(ii1,ii2,ii3) = 0.0
                  af(ii1,ii2,ii3) = 0.0
                END DO
              END DO
            END DO
 
            CALL ASYMTOP(jint,itaui,ntaui,numji,ai,egi)
            CALL ASYMTOP(jfnl,itauf,ntauf,numjf,af,egf)

          END IF
c
c     For symmetric top -- itai is K
c                          itaf is K'
c     For asymmetric top - itai is Tau
c                          itaf is Tau'
c
          DO iti = 1 , j21
            IF ( mol.EQ.1 ) THEN
              itai = jint1 - iti
              eji = (bz*xjk1+(az-bz)*itai*itai)
            ELSE
              itai = itaui(iti)
              eji = egi(iti)
              numti = ntaui(iti)
            END IF
 
            DO itf = 1 , jp21
              IF ( mol.EQ.1 ) THEN
                itaf = jfnl1 - itf
                ejf = (bz*xjk2+(az-bz)*itaf*itaf)
              ELSE
                itaf = itauf(itf)
                ejf = egf(itf)
                numtf = ntauf(itf)
              END IF
              rote = ejf - eji
c
c     Rote is the energy change of the transition J,K(Tau)-->J^,K^(Tau^)
c     Transitions are allowed only if energy is greater than rote
c 
              IF ( en.GT.rote/13.606 ) THEN
c
c     Selection rule:
c                     Transition allowed only if
c                     delta K (symmetric top) or
c                     delta Tau (asymmetric top) is EVEN.
c
                IF ( IPARITY(itaf-itai).EQ.1 ) THEN
                  xkf = SQRT(en-(ejf-eji)/13.606)
c
c     Call mmatrx to generate the M-matrix
c
                  CALL MMATRX

c
c     Call ali for generate the Al coefficient of the DCS's expansions
c
                  CALL ALI
 
                  DO lll = 1 , lbig
                    lll1 = lll - 1
                    xll = FLOAT(lll1+lll1+1)
                    albcc(1,lll) = albcc(1,lll)*fact*xll
                    IF ( mol.EQ.1 ) albcc(1,lll) = albcc(1,lll)/xjp21
                    albcc(2,lll) = albcc(2,lll)*fact*xll
                    IF ( mol.EQ.1 ) albcc(2,lll) = albcc(2,lll)/xjp21
                  END DO
c
c     Print header and rotational energies summary
c
                  IF ( iprint.EQ.2 ) THEN
                    WRITE (6,99002)
                    WRITE (6,99021) jint , itai , jfnl , itaf
                    WRITE (6,99002)
                    WRITE(6,99043)
                    WRITE(6,99044) eji*1.0E+03 , ejf*1.0E+03 , 
     &                             rote*1.0E+03 , xki, xkf
                  END IF
c
c     Call cross for generate DCS ,CST and CSM
c
                  CALL CROSS(albcc,en,xki,xkf,dcsb,dcs,csm,cst)
 
                  DO i = 1 , NANG
                    sumdcs(1,i) = sumdcs(1,i) + dcs(i)
                    IF ( ipol.NE.0 ) sumdcs(2,i) = sumdcs(2,i) + dcsb(i)
                  END DO
 
                  IF ( ipol.NE.0 ) THEN
                    sumcsmb = sumcsmb + csmb
                    sumcstb = sumcstb + cstb
                  END IF
 
                  roteff = roteff + cst*rote
 
                  sumcsm = sumcsm + csm
                  sumcst = sumcst + cst
                  sumcst1 = sumcst1 + cst1
                  sumcsm1 = sumcsm1 + csm1
                END IF
              END IF
            END DO
          END DO
 
c +++++++++++++++ exit from K e K^ loop+++++++++++++++++++++++++

c
c     Calculate  DCS, CS and CSM for j-->j' summed 
c     over K(Tau) and K^(Tau^)
c
          WRITE (6,*)
          WRITE (6,99002)
          WRITE (6,99006) jint , jfnl
          WRITE (6,99002)
          WRITE (6,99019)
 
          DO i = 1 , NANG
            totdcs(1,i) = totdcs(1,i) + sumdcs(1,i)
            WRITE (6,99020) th(i) , sumdcs(1,i)
            IF ( ipol.NE.0 ) totdcs(2,i) = totdcs(2,i) + sumdcs(2,i)
          END DO
 
          IF ( ipol.NE.0 ) THEN
            totcsmb = totcsmb + sumcsmb
            totcstb = totcstb + sumcstb
          END IF
 
          WRITE (6,99022) sumcsm , sumcst
 
          totcsm = totcsm + sumcsm
          totcst = totcst + sumcst
          totcst1 = totcst1 + sumcst1
          totcsm1 = totcsm1 + sumcsm1
        END IF
c
c***** End code for SYMMETRIC/ASYMMETRIC TOP molecule (MOL=1,2) *****
c
 100  END DO

c+++++++++++++++++ exit from J e J^ loop++++++++++++++++++++++++

c
c     Print DCS,CS and CSM summed on all transitions for all molecules
c
      WRITE (6,99002)
      WRITE (6,99030)
      WRITE (*,99002)

      IF ( iprint.EQ.2 ) THEN
        IF ( ipol.NE.0 ) THEN
          WRITE (6,99031) 
           DO i = 1 , NANG
             WRITE (6,99033) th(i) , totdcs(1,i), totdcs(2,i)
           END DO
           WRITE (6,99034)
           WRITE (6,99037) totcsmb , totcstb
           WRITE (6,99035)
           WRITE (6,99037) totcsm1 , totcst1
        ELSE
          WRITE (6,99032)
           DO i = 1 , NANG
             WRITE (6,99033) th(i) , totdcs(1,i)
           END DO
        END IF
      END IF

      WRITE (6,99036)
      WRITE (6,99037) totcsm , totcst

      IF ( iprint.EQ.2 ) WRITE (6,99038) roteff/totcst
      

99000 FORMAT (20x,'******************************'/,
     &        20x,'*                            *'/,
     &        20x,'*        POLYDCS V1.0        *'/,
     &        20x,'*         June 1998          *'/,
     &        20x,'*                            *'/,
     &        20x,'******************************'//)
99001 FORMAT (a80)
99002 FORMAT (2x,74('*'))
99003 FORMAT (a6,3I6)
99004 FORMAT (/2x,'Irreducible Representations (IR)',/
     &        /8x,'IR #',5x,'IR Type',6x,'KSIZ',6x,'KSIZB',5x,
     &        'IR Read (0=No)')
99005 FORMAT (6x,i5,5x,a6,2(6x,I5),8x,I5)
99006 FORMAT (6x,'Results for ',I2,' --> ',I2,
     &        "  J-J' transition summed on all K and K'")
99007 FORMAT (//2x,'MOLECULAR PROPERTIES:',/,
     &        /2x,'Rotational Constant(s) in eV',
     &        /6x,'(Az) =',E13.6,5x,'(Bz) =',E13.6,5x,'(Cz) =',E13.6)
99008 FORMAT (/2x,'No. of Scattering K-matrices         (NSCAT)  =',i4,
     &        /2x,'Max Value of L (from K-matrices)     (LMAX)   =',i4,
     &        /2x,'Max Value of L (Born, <>0 if IPOL=1) (LMAXB)  =',i4,
     &        /2x,'Max Value of DCS PLeg Expansion      (LBIG)   =',i4,
     &        /2x,'No. of Rotational Transitions        (NT)     =',i4,
     &        /2x,'Print Flag (0:Min, 1:Norm, 2:Full)   (IPRINT) =',i4)
99009 FORMAT (/2x,'Vibrational VCC calculation:',/,2x,'  nu = ',i3,
     &        " nu' = ",i3,' Frequency = ',f10.4,' eV')
99010 FORMAT (//2x,'*** Rotational Transitions to be calculated ***',
     &        //6x,'  #       ji       jf',/7x,'---------------------')
99011 FORMAT (3(5x,i4))
99014 FORMAT (/2x,'Energy (eV)  = ',f12.7,
     &        /2x,'Energy (ryd) = ',f12.7,
     $        '  sqrt(Energy) [xki] = ',f12.7,//)
99016 FORMAT (/2x,'Polar molecule found (IPOL=',i2,')'/)
99017 FORMAT (/2x,'Electric moments in a.u.',
     &        /2x,'Dipole          Theor.(DPOLT) =',f9.5,
     &                         '    Exp. (DPOLE) =',f9.5,
     &        /2x,'Quadrupole               (Q1) =',f9.5,
     &                         '            (Q2) =',f9.5,
     &        /2x,'Octupole               (OCT1) =',f9.5,
     &                         '          (OCT2) =',f9.5)
99018 FORMAT (/2x,'Electric moments and polarizabilities in a.u.',
     &        /2x,'Polarizabilities     (ALPHA0) =',f9.5,
     &                         '        (ALPHA2) =',f9.5,
     &        /2x,'Dipole           Theor.(DPOL) =',f9.5,
     &                         '   Exp. (DPOLEX) =',f9.5,
     &        /2x,'Quadrupole             (QUAD) =',f9.5,
     &                         '        (QUADEX) =',f9.5)
99019 FORMAT (/6x,'THETA(deg)',10x,'DCS(10**-16 cm**2/sr)')
99020 FORMAT (6x,f7.2,16x,E15.8)
99021 FORMAT (6x,'Results for J=',I2,' K=',I2, 
     &        " --> J'=",I2," K'=",I2,' transition')
99022 FORMAT (/2x,'Partial Momentum Transfer Cross Section       =',
     &        f13.8,' 10**-16 cm**2',
     &        /2x,'Partial Integrated Differential Cross Section =',
     &        f13.8,' 10**-16 cm**2',/)
99024 FORMAT (f12.4,5E16.6)
99025 FORMAT (//,' *** POLYDCS ERROR in main - - PROGRAM STOP ***',/,
     &        ' NSCAT = ',i5,' greater than NSYM = ',i5)
99026 FORMAT (//,' *** POLYDCS ERROR in main - - PROGRAM STOP ***',/,
     &        ' LBIG  = ',i5,' greater than NLAM = ',i5)
99027 FORMAT (//,' *** POLYDCS ERROR in main - - PROGRAM STOP ***',/,
     &        ' NT    = ',i5,' greater than NDCS = ',i5)
99028 FORMAT (//,' *** POLYDCS ERROR in main - - PROGRAM STOP ***',/,
     &        ' KSIZ  = ',i5,' greater than NPW  = ',i5,
     &        ' for IR # ',i5)
99029 FORMAT (//,' *** POLYDCS ERROR in main - - PROGRAM STOP ***',/,
     &        ' KSIZB  = ',i5,' greater than NPW  = ',i5,
     &        ' for IR # ',i5)
99030 FORMAT (20x,'DATA SUMMED ON ALL TRANSITIONS') 
99031 FORMAT (/6x,'THETA(deg)',8x,'DCSC(10**-16 cm**2/sr)',
     &        5x,'DCSB(10**-16 cm**2/sr)')
99032 FORMAT (/6x,'THETA(deg)',8x,'DCSC(10**-16 cm**2/sr)')
99033 FORMAT (6x,f7.2,2(12x,E15.8))
99034 FORMAT (//2x,'Quantity Summed Over States (Born):')
99035 FORMAT (//2x,'Quantity Summed Over States (Close-Coupling):')
99036 FORMAT (//2x,'Quantity Summed Over States (Complete Rep.):')
99037 FORMAT (/2x,'  Momentum Transfer Cross Section (CSM)      =',
     &        f15.8,' 10**-16 cm**2',
     &        /2x,'  Integrated Differential Cross Section (CS) =',
     &        f15.8,' 10**-16 cm**2')
99038 FORMAT (//2x,'Rotational Efficiency =',f15.8)
99039 FORMAT (/2x,'*** WARNING ***',
     &        /2x,'==> Lmax = ',i4,'. Lmax from k-matrix file = ',i4)
99040 FORMAT (2x,'    IR # ',i3,'    K-matrix size reduced to ',
     &        i4,' x ',i4)
99041 FORMAT (/2x,'*** WARNING ***',
     &        /2x,'==> Lmaxb = ',i4,'. Lmaxb from Blm file = ',i4)
99042 FORMAT (2x,'    IR # ',i3,'    Blm size reduced to ',
     &        i4,' x ',i4)
99043 FORMAT (//2x,'Rotational energy spacing and wavevectors:')
99044 FORMAT (/2x,"  Initial Energy          [Ejk] =",E13.6,' meV',
     &        /2x,"  Final Energy           [Ejk'] =",E13.6,' meV',
     &        /2x,"  Delta E            [Ejk'-Ejk] =",E13.6,' meV',
     &        /2x,"  Initial wavevector        [k] =",E13.6,' Ryd',
     &        /2x,"  Final wavevector         [k'] =",E13.6,' Ryd',//)

      END

