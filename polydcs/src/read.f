      SUBROUTINE KREAD(en)
      IMPLICIT REAL*4 (A-H, O-Z)
 
c
c     Read the K-matrices and generate the T-matrices
c
      INCLUDE 'par.h'

      CHARACTER*6  stat
      CHARACTER*80 blmfn , kmatfn

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
      COMMON /TMAT  / tmr(NSYM,NPW,NPW) , tmi(NSYM,NPW,NPW)

c
c     Read the K-matrices from file
c
      OPEN (UNIT=7,FILE=kmatfn,STATUS='old')
      imax = 0
      DO i = 1 , nscat
        IF ( mult(i).EQ.0 ) THEN
          DO j = 1 , ks
            DO k = 1 , ks
              akm(i,j,k) = akm(i-1,j,k)
            END DO
          END DO
        ELSE
          READ (7,99001) nter , energy
          IF ( ksiz(i) .NE. nter ) THEN
            WRITE (6,99003) nter , ksiz(i) , nscat
            STOP
          END IF
          IF ( ABS(energy-en*13.606).GT.10E-04) THEN
            WRITE (6,99008) energy , en*13.606 , nscat
            STOP
          END IF

          ks = ksiz(i)
          DO j = 1 , ks
            DO k = 1 , ks , 4
              imax = k + 3
              IF ( imax.GT.ks ) imax = ks
              READ (7,99002) (akm(i,j,l),l=k,imax)
            END DO
          END DO
        END IF
      END DO
c
c     Case of non-polar non-linear molecules:
c     generate T-matrices and exit kread.
c
      IF ( (IPOL.EQ.0) .AND. (MOL.LT.3) ) THEN
        CALL TMATRX(akm,tmr,tmi,ksiz)
        RETURN
      END IF
c
c     Case of linear non-polar molecules:
c     - Set the min/max value of L for Blm (Bl in this case ...).
c     - Generate T-matrices and exit kread.
c
      IF ( MOL.EQ.5 ) THEN
        mnlamda = 0
        mxlamda = 0
        DO i = 1 , nscat
          mnlamda = MIN(mnlamda,ksiz(i))
          mxlamda = MAX(mxlamda,ksiz(i))
        END DO 
        DO i = 1 , nscat
          lk(1,i) = i
          lk(2,i) = mxlamda
        END DO
        CALL TMATRX(akm,tmr,tmi,ksiz)
        RETURN
      END IF
c
c     Case of linear polar molecules:
c     - Increase the no. of IRs up to Lmax and set the K-matrix elements
c       by calling the ACBORN routine.
c     - Increase the dimension of each IR up to Lmaxb*2
c     - Define the new indices for each IR (no. of them and dimensions)
c     - Generate T-matrices and exit kread.
c
c     Note that in this case the lk(,) matrix is used for 2 rows only:
c     - lk(1,NPW=NSYM) = li() initial L value for each IR 
c     - lk(2,NPW=NSYM) = lf() final L value for each IR 
c     where the lk(,) has to be used by the LINCROSS subroutine only.
c
      IF ( (IPOL.EQ.1) .AND. (MOL.EQ.4) ) THEN
c
c************* To increase the number of molecular IR's **************
c
        IF ( lmax.GT.NPW ) THEN
          WRITE (6,99004) lmax , NPW
          STOP
        END IF 
        IF ( lmax.GT.NSYM ) THEN
          WRITE (6,99005) lmax , NSYM
          STOP
        END IF 
        IF ( lmax.GT.nscat ) THEN
          DO m = nscat + 1 , lmax
            DO i = m - 1 , lmax - 1
              CALL ACBORN(m-1,i,i,en,accborn)
              akm(m,i-m+2,i-m+2) = -0.5*accborn
              IF ( i.EQ.lmax-1 ) GO TO 10
              CALL ACBORN(m-1,i,i+1,en,accborn)
              akm(m,i-m+2,i-m+3) = -0.5*accborn
              IF ( i.NE.lmax-2 ) THEN
                CALL ACBORN(m-1,i,i+2,en,accborn)
                akm(m,i-m+2,i-m+4) = -0.5*accborn
              END IF
            END DO
 10       END DO
        END IF
c
c************** To increase the dimension of each IR's ***************
c
        IF ( (2*lmaxb).GT.NPW ) THEN
          WRITE (6,99006) (2*lmaxb) , NPW
          STOP
        END IF
        IF ( (2*lmaxb).GT.NSYM ) THEN
          WRITE (6,99007) (2*lmaxb) , NSYM
          STOP
        END IF
        IF ( (2*lmaxb-2).GT.(lmax-2) ) THEN
          llty = 2*lmaxb - 2
          DO m = 1 , lmax
            DO i = lmax - 2 , llty
              IF ( i.NE.lmax-2 ) THEN
                IF ( i.NE.lmax-1 ) THEN
                  CALL ACBORN(m-1,i,i,en,accborn)
                  akm(m,i-m+2,i-m+2) = -0.5*accborn
                END IF
                IF ( i.EQ.(2*lmaxb-2) ) GO TO 20
                CALL ACBORN(m-1,i,i+1,en,accborn)
                akm(m,i-m+2,i-m+3) = -0.5*accborn
              END IF
              IF ( i.NE.(2*lmaxb-3) ) THEN
                CALL ACBORN(m-1,i,i+2,en,accborn)
                akm(m,i-m+2,i-m+4) = -0.5*accborn
              END IF
            END DO
 20       END DO
        END IF
c
c************* Reassign indeces and call tmatrx *****
c
        nscat = lmax
        mnlamda = 0
        mxlamda = nscat - 1
        DO i = 1 , nscat
          ksiz(i) = 2*lmaxb - i
          lk(1,i) = i
          lk(2,i) = 2*lmaxb - 2
        END DO
c
c     This is for the lower part of the Kmatrix
c
        DO i = 1 , nscat
          DO j = 2 , ksiz(i)
            DO l = 1 , j - 1
              akm(i,j,l) = akm(i,l,j)
            END DO
          END DO
        END DO

        CALL TMATRX(akm,tmr,tmi,ksiz)
   
c
c     Reassign the new indeces to T-matrices and exit
c
 
        llmax = 2*lmaxb - 2
        DO i = 2 , nscat
          DO j = llmax - i + 1 , 1 , -1
            DO k = llmax - i + 1 , 1 , -1
              tmr(i,j+i-1,k+i-1) = tmr(i,j,k)
              tmi(i,j+i-1,k+i-1) = tmi(i,j,k)
            END DO
          END DO
          DO j = 1 , i - 1
            DO k = 1 , i - 1
              tmr(i,j,k) = 0.0
              tmi(i,j,k) = 0.0
            END DO
          END DO
        END DO
   
        RETURN

      END IF


      RETURN
 
99001 FORMAT (15x,i5,e16.8)
99002 FORMAT (1x,4(1x,e18.10))
99003 FORMAT (//,' *** POLYDCS ERROR in KREAD - - PROGRAM STOP ***',/,
     &        ' Read KSIZ is',i5,' not equal to KSIZ = ',i5,
     &        ' defined in input for IR # ',i5)
99004 FORMAT (//,' *** POLYDCS ERROR in KREAD - - PROGRAM STOP ***',/,
     &        ' Lmax = ',i4,' greater than NPW = ',i4,/,
     &        ' No storage available to increase Born K-matrices')
99005 FORMAT (//,' *** POLYDCS ERROR in KREAD - - PROGRAM STOP ***',/,
     &        ' Lmax = ',i4,' greater than NSYM = ',i4,/,
     &        ' No storage available to increase Born K-matrices')
99006 FORMAT (//,' *** POLYDCS ERROR in KREAD - - PROGRAM STOP ***',/,
     &        ' 2*Lmaxb = ',i4,' greater than NPW = ',i4,/,
     &        ' No storage available to increase Born K-matrices')
99007 FORMAT (//,' *** POLYDCS ERROR in KREAD - - PROGRAM STOP ***',/,
     &        ' 2*Lmaxb = ',i4,' greater than NSYM = ',i4,/,
     &        ' No storage available to increase Born K-matrices')
99008 FORMAT (//,' *** POLYDCS ERROR in KREAD - - PROGRAM STOP ***',/,
     &        ' Read ENERGY is',f7.3,' eV not equal to EN = ',f10.5,
     &        ' eV defined in input for IR # ',i5)

      END


      SUBROUTINE BLMREAD
      IMPLICIT REAL*4 (A-H, O-Z)
 
c
c     Read blm data for spherical, symmetric and 
c     asymmetric top molecules
c

      INCLUDE 'par.h'

      CHARACTER*6  stat
      CHARACTER*80 blmfn , kmatfn
      DIMENSION mdummy(NPW) , brdummy(NPW) , bidummy(NPW)

      COMMON /BLMCOF/ br(NSYM,NPW,NPW) , bi(NSYM,NPW,NPW)
      COMMON /INPCHR/ stat(NSYM) , blmfn , kmatfn
      COMMON /INPMAT/ akm(NSYM,NPW,NPW) ,
     &                mult(NSYM) , ksiz(NSYM) , ksizb(NSYM) ,
     &                lk(NSYM,NPW) , num(NSYM,NPW) , ml(NSYM,NPW,NPW)
      COMMON /INPPAR/ ipol , iprint , nt , nscat , mol ,
     &                lmax , lmaxb , lbig ,
     &                alpha0 , alpha2 , dpol , dpolex , quad , quadex ,
     &                dpolt , dpole , q1 , q2 , oct1 , oct2


      DATA sqr2/1.414213562E0/

      OPEN (UNIT=8,FILE=blmfn,STATUS='old')
c
c     Read in the blm's and convert them from real to complex.
c
      imax = 0
      DO i = 1 , nscat
        READ (8,99001) iblmsize

        IF ( ksizb(i) .NE. iblmsize ) THEN
          WRITE (6,99002) iblmsize , ksizb(i) , nscat
          STOP
        END IF

        ksizb(i) = iblmsize
        ks = ksizb(i) 
 
        DO j = 1 , ks
          READ (8,*) lk(i,j) , num(i,j)
          ms = num(i,j)
          DO jj = 1 , ms , 4
            imax = jj + 3
            IF ( imax.GT.ms ) imax = ms
            READ (8,*) (ml(i,j,k),br(i,j,k),k=jj,imax)
          END DO
          ip = 0
          DO k = 1 , ms
            IF ( ml(i,j,k).EQ.0 ) THEN
              ip = ip + 1
              mdummy(ip) = ml(i,j,k)
              brdummy(ip) = br(i,j,k)
              bidummy(ip) = 0.0E0
            ELSE
              ip = ip + 2
              IF ( ml(i,j,k).GT.0 ) THEN
                mdummy(ip) = -ml(i,j,k)
                mdummy(ip-1) = ml(i,j,k)
                coeff = FLOAT((-1)**ABS(mdummy(ip)))
                brdummy(ip) = br(i,j,k)/sqr2
                brdummy(ip-1) = coeff*br(i,j,k)/sqr2
                bidummy(ip) = 0.0E0
                bidummy(ip-1) = 0.0E0
              ELSE
                mdummy(ip) = ml(i,j,k)
                mdummy(ip-1) = -ml(i,j,k)
                coeff = FLOAT((-1)**ABS(mdummy(ip)+1))
                brdummy(ip) = 0.0E0
                brdummy(ip-1) = 0.0E0
                bidummy(ip) = br(i,j,k)/sqr2
                bidummy(ip-1) = coeff*br(i,j,k)/sqr2
              END IF
            END IF
          END DO
          num(i,j) = ip
          ms = num(i,j)
          DO k = 1 , ms
            ml(i,j,k) = mdummy(k)
            br(i,j,k) = brdummy(k)
            bi(i,j,k) = bidummy(k)
          END DO
        END DO
      END DO

99001 FORMAT (22x,i6)
99002 FORMAT (//,' *** POLYDCS ERROR in BLMREAD - - PROGRAM STOP ***',/,
     &        ' Read KSIZB is',i5,' not equal to KSIZB = ',i5,
     &        ' defined in input for IR # ',i5)

      END

