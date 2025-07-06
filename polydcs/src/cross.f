      SUBROUTINE CROSS(an,en,xki,xkf,dcsb,dcs,csm,cst)
      IMPLICIT REAL*4 (A-H, O-Z)
 
      INCLUDE 'par.h'
 
      DIMENSION an(2,NLAM)
      DIMENSION dcs1(NANG) , dcs2(NANG) , dcs(NANG) , dcsb(NANG)
 
      COMMON /ANGDAT/ th(NANG) , pl(NLAM,NANG)
      COMMON /INPPAR/ ipol , iprint , nt , nscat , mol ,
     &                lmax , lmaxb , lbig ,
     &                alpha0 , alpha2 , dpol , dpolex , quad , quadex ,
     &                dpolt , dpole , q1 , q2 , oct1 , oct2
      COMMON /JKTDAT/ jmin, jmax , mnlamda , mxlamda ,
     &                jint , jfnl , itai , itaf , numti , numtf
      COMMON /TOTDCS/ cst1 , csm1 , cst2 , csm2 , cstb , csmb
 
      DATA pi4 , unit /12.5637061E0 , 3.5712026E0/
 
      csmb = 0.0
      cstb = 0.0
      csm2 = 0.0
      cst2 = 0.0
 
      jdel = jfnl - jint
      jif = ABS(jdel)
c     jj1 = jint + jfnl + 1
c     jj2 = jint + jfnl + 2
 
      ita = 0
      IF ( mol.EQ.1 .OR. mol.EQ.2 ) ita = IABS(itaf+itai)
 
c
c     If jval=1 does the Born correction
c
      jval = 0
      IF ( ipol.NE.0 ) THEN
        IF ( (mol.EQ.1 .OR. mol.EQ.2) .AND. jif.EQ.1 .AND. ita.EQ.0 )
     &       jval = 1
        IF ( mol.EQ.4 .AND. (jif.EQ.1 .OR. jif.EQ.2) ) jval = 1
      END IF
 
      xkr = xkf/xki
 
      IF ( jval.NE.0 ) THEN

        xkpl = xki + xkf
        xkmn = ABS(xki-xkf)
        sumb = pi4*2.0/3.0*ALOG(xkpl/xkmn)/xki**2
        IF ( mol.GT.3 ) THEN
          d2 = dpolex*dpolex
          q2 = quadex*quadex
        ELSE
          d2 = dpole*dpole
        END IF
 
        xk3 = 3.0*en
        a1 = 2.0*d2/xk3
        a2 = a1*pi4
        a3 = 4.0*d2/3.0
        IF ( mol.GT.3 ) THEN
c         a4 = 4.0*pi4*q2/45.0
          a5 = alpha2**2
        END IF
        xk1 = (en+xkf*xkf)
        xk2 = 2.0*xki*xkf

      END IF

      IF ( iprint.EQ.2 ) THEN
        IF ( jval.EQ.0 ) WRITE (6,99003)
        IF ( jval.EQ.1 ) WRITE (6,99002)
      END IF
 
c
c     Print the Al coefficients
c
      IF ( iprint.EQ.2 ) THEN
        DO l = 1 , lbig
          lla = l - 1
          IF ( jval.EQ.1 ) THEN
            WRITE (6,99006) lla , an(1,l) , an(2,l)
          ELSE
            WRITE (6,99006) lla , an(1,l)/unit
          END IF
        END DO
      END IF
 
      DO i = 1 , NANG
        dcsb(i) = 0.0
        dcs2(i) = 0.0
        csth = COS(th(i)*0.0174532E0)
        sum1 = an(1,1)
        IF ( jval.EQ.1 ) THEN
          sum2 = an(2,1)
        END IF
        DO l = 2 , lbig
          lla = l - 1
          pll = pl(lla,i)
          an1 = an(1,l)
          sum1 = sum1 + an1*pll
          IF ( jval.EQ.1 ) THEN
            an2 = an(2,l)
            sum2 = sum2 + an2*pll
          END IF
        END DO
        dcs1(i) = sum1/unit*xkr
        IF ( jval.EQ.1 ) dcs2(i) = sum2/unit*xkr
c
c     This is for polar molecules
c
        IF ( jval.NE.0 ) THEN
 
c     This is for the dipole transition correction
 
          IF ( jif.EQ.1 ) THEN
            IF ( i.NE.1 ) THEN
              dcsb(i) = a3*xkr/((xk1-xk2*csth)*unit)
            ELSE
              xfg = (SQRT(xki)-SQRT(xkf))
              dcsb(i) = a3*xkr/((xfg**2)*unit)
            END IF
            dcsb(i) = dcsb(i)*(2*jfnl+1)*(WIG3J0(jint,jfnl,1))**2
          END IF
c
c     This is for the quadrupole and alpha2 transition correction.
c     At the moment, implemented for linear molecules only
c
          IF ( jif.EQ.2 ) THEN
c
c     aca is the contribute of the alpha2
c
            aca = a5*xkr*(xk1-xk2*csth)*(WIG3J0(jfnl,jint,2))**2
            aca = ((pi4/64)**2)*aca/unit
            cist = (WIG3J0(jfnl,jint,2))**2
            aaa = ((pi4/64)**2)*a5*xkr*xk1*cist/unit
            bbb = ((pi4/64)**2)*a5*xkr*xk2*cist/unit
c
c     ala is the contribute of the quadrupole
c
            ala = 4.0*q2*xkr/(9.0*5.0*unit)
            ala = ala*(2*2+1)*(WIG3J0(jint,jfnl,2))**2
            IF ( i.EQ.1 .AND. ala.LT.0.0 ) ala = ABS(ala)
c
c     aga is the contribute of the quadrupole and alpha2 togheter
c
            aga = 2*SQRT(aca*ala)
            dcsb(i) = ala + aca + aga
          END IF
        END IF
c
c     dcs(i) is the differential cross-section in the closure
c     formula or in the close-coupling (if ipol=0) approx.
c
        dcs(i) = dcsb(i) + dcs1(i) - dcs2(i)

      END DO
c
c     To print differential cross-sections
c
      IF ( iprint.EQ.2 .AND. mol.NE.3 ) THEN
        WRITE (6,99004) 
        IF ( jval.EQ.0 ) THEN
          DO i = 1 , NANG
            WRITE (6,99008) th(i) , dcs(i)
          END DO
        ELSE IF ( jval.EQ.1 ) THEN
          WRITE (6,99005) 
          DO i = 1 , NANG
            WRITE (6,99008) th(i) , dcsb(i), dcs1(i), dcs2(i), dcs(i)
          END DO
        END IF
      END IF
c
c
c         Momentum transfer cross sections and total
c         (integrated) cross section in born and close
c         coupling theory (born only for polar molecules)
c
c
c          csm1  momentum transfer cross-section in close-coupling
c          cst1             total cross-section  "     "
c          csm2  momentum transfer cross-section in unitarised born
c          cst2              total cross-section  "     "
c          csmb  momentum transfer cross-section in analytic born
c          cstb             total cross-section  "     "
c          csm   momentum transfer cross-section in closure formula
c          cst              total cross-section  "     "
c
      csm1 = (an(1,1)-an(1,2)/3.0)/unit*pi4
      cst1 = pi4*an(1,1)/unit
      IF ( jval.NE.0 ) THEN
        csm2 = (an(2,1)-an(2,2)/3.0)/unit*pi4
        cst2 = pi4*an(2,1)/unit
 
        IF ( jif.EQ.1 ) THEN
          csmb = a2/unit
          cstb = sumb*d2/unit
        END IF
 
        IF ( jif.EQ.2 ) THEN
          abm = 1.0 - bbb/aaa
          ab = 1.0 + bbb/aaa
          cyou = SQRT(ab**3) + aaa/(5.0*bbb)*(SQRT(abm**5)-SQRT(ab**5))
          csmb = 2.0*(aaa+ala+bbb/3) + 8.0*aaa*SQRT(ala*aaa)/(3.0*bbb)
     &           *cyou
          csmb = pi4/2.0*csmb
          tyou = SQRT((1-bbb/aaa)**3) - SQRT((1+bbb/aaa)**3)
          cstb = 2.0*(aaa+ala) - 4.0/3.0*SQRT(aaa*ala)*(aaa/bbb)*tyou
          cstb = cstb*pi4/2.0
        END IF

      END IF
 
 
      csm = csmb + csm1 - csm2
      cst = cstb + cst1 - cst2
c
c     Print momentum transfer cross sections
c
      IF ( iprint.EQ.2 ) WRITE (6,99007)
      IF ( jval.EQ.0 ) THEN
        IF ( iprint.EQ.2 ) WRITE (6,99012) csm
      ELSE
        IF ( iprint.EQ.2 ) WRITE (6,99011) csm1 , csm2 , csmb , csm
      END IF
c
c     Print integrated transfer cross sections
c
      IF ( iprint.EQ.2 ) WRITE (6,99010)
      IF ( jval.EQ.0 ) THEN
        IF ( iprint.EQ.2 ) WRITE (6,99012) cst
      ELSE
        IF ( iprint.EQ.2 ) WRITE (6,99011) cst1 , cst2 , cstb , cst
      END IF
 
      RETURN
c
c     Format statements
c
99002 FORMAT (/2x,'Close-Coupling (CC) and Unitarised Born (BORN)', 
     &        ' AL(L=0-Lbig) coefficients',/2x,72('-'),/,28x,'CC',
     &        22x,'BORN')
99003 FORMAT (/6x,'Close-Coupling AL(L=0-Lbig) coefficients',/6x,
     &        40('-'))
99004 FORMAT (/6x,'State-to-State Differential Cross Sections',
     &        /6x,42('-'),/,8x,'THETA(deg)',3x,'DCS(10**-16 cm**2/sr)')
99005 FORMAT (21x,'Born Analitic',2x,'Close-Coupling',2x,
     &        'Born Unit.',2x,'Complete Rep.')
99006 FORMAT (10x,'A(L=',i2,')',e20.8,5x,e20.8)
99007 FORMAT (/2x,'Momentum Transfer Cross Section (in units of',
     &        ' 10**-16 cm**2)')
99008 FORMAT (8x,f7.2,5x,4(2x,E12.4))
99010 FORMAT (/2x,'Total Cross Section (in units of 10**-16 cm**2)')
99011 FORMAT (/5x,'In Close-Coupling        =',e13.6,
     &        /5x,'In Unitarised Born       =',e13.6,
     &        /5x,'In Analytic Born         =',e13.6,
     &        /5x,'In Closure Approximation =',e13.6)
99012 FORMAT (/5x,'In Close-Coupling        =',e13.6)
      END
 
