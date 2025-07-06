      SUBROUTINE BKPRINT
      IMPLICIT REAL*4 (A-H, O-Z)
 
c
c     Write out Blms and K-matrices data.
c
      INCLUDE 'par.h'
 
      CHARACTER*6  stat
      CHARACTER*80 blmfn , kmatfn
      
      COMMON /BLMCOF/ br(NSYM,NPW,NPW) , bi(NSYM,NPW,NPW)
      COMMON /INPCHR/ stat(NSYM) , blmfn , kmatfn
      COMMON /INPMAT/ akm(NSYM,NPW,NPW) ,
     &                mult(NSYM) , ksiz(NSYM) , ksizb(NSYM) ,
     &                lk(NSYM,NPW) , num(NSYM,NPW) , ml(NSYM,NPW,NPW)
      COMMON /INPPAR/ ipol , iprint , nt , nscat , mol ,
     &                lmax , lmaxb , lbig ,
     &                alpha0 , alpha2 , dpol , dpolex , quad , quadex ,
     &                dpolt , dpole , q1 , q2 , oct1 , oct2


c
c     For linear molecules print only k-matrices
c
      IF ( mol.GT.3 ) THEN
        WRITE (6,99002)
        DO i = 1 , nscat
          CALL PRTKML(i)
        END DO
        RETURN
      END IF
c
c     Here is the printout of blms and k-matrices
c     for all the other molecular types.
c

      WRITE (6,99001)
      DO i = 1 , nscat
        WRITE (6,*)
        WRITE (6,99003)
        WRITE (6,99004)
        WRITE (6,99003)
        WRITE (6,99005) stat(i) , ksiz(i) , ksiz(i)
        ks = ksiz(i)
        DO j = 1 , ks
          lakj = lk(i,j)
          WRITE (6,99006) lakj
          nma = num(i,j)
          DO k = 1 , nma
            WRITE (6,99007) ml(i,j,k),br(i,j,k),bi(i,j,k),i,j,k 
          END DO
        END DO
      END DO
c
c     Below istep is fixed to 6. It can be changed as needed together
c     with the format statement at 99009.
c
      istep = 6
      DO k = 1 , nscat
        IF ( mult(k).NE.0 ) THEN
          WRITE (6,99008) stat(k)
          ks = ksiz(k)
          DO J = 1 , ks , istep
            imin = j
            imax = j + istep - 1
            IF ( imax.GT.ks ) imax = ks
            WRITE (6,99009) (lk(k,l),l=imin,imax)
            DO I = 1 , ks
              WRITE (6,99010) lk(k,i) , (akm(k,i,l),l=imin,imax)
            END DO
            WRITE (6,*)
          END DO
        END IF
      END DO
 
99001 FORMAT (///2x,'*** Blms and K-matrices elements printout ***',//)
99002 FORMAT (///2x,'*** K-matrices elements printout ***',//)
99003 FORMAT (2x,72('-'))
99004 FORMAT (2x,'   IR   ',2x,'l',4x,'m',6x,'b(real)',8x,'b(imag)')
99005 FORMAT (2x,a6,' dimension = ',i2,' x ',i2)
99006 FORMAT (11x,i2)
99007 FORMAT (15x,i3,2F15.8,3I3)
99008 FORMAT (/2x,'K-matrix for state =',a6,/)
99009 FORMAT (8x,6('    L = ',i3))
99010 FORMAT (' L = ',i3,20(2x,f9.5))

      END

      SUBROUTINE PRTKML(K)
      IMPLICIT REAL*4 (A-H, O-Z)

      INCLUDE 'par.h'

      CHARACTER*6  stat
      CHARACTER*80 blmfn , kmatfn

      COMMON /INPCHR/ stat(NSYM) , blmfn , kmatfn
      COMMON /INPMAT/ akm(NSYM,NPW,NPW) ,
     &                mult(NSYM) , ksiz(NSYM) , ksizb(NSYM) ,
     &                lk(NSYM,NPW) , num(NSYM,NPW) , ml(NSYM,NPW,NPW)

      istep = 6
      WRITE (6,99008) stat(k)
      ks = ksiz(k)
      DO j = 1 , ks , istep
        imin = j
        imax = j + istep - 1
        IF ( imax.GT.ks ) imax = ks
        WRITE (6,99009) (l,l=imin,imax)
        DO i = 1 , ks
          WRITE (6,99010) i , (akm(k,i,l),l=imin,imax)
        END DO
        WRITE (6,*)
      END DO

99008 FORMAT (/2x,'K-matrix for state =',a6,/)
99009 FORMAT (8x,6('    # = ',i3))
99010 FORMAT (' # = ',i3,20(2x,f9.5))

      END
