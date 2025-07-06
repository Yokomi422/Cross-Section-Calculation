      SUBROUTINE TMATRX(akmat,atr,ati,kks)
      IMPLICIT REAL*4 (A-H, O-Z)
 
      INCLUDE 'par.h'

c
c     This subroutine generate T-matrix from K-matrix
c

      DIMENSION kks(NSYM)
      DIMENSION ak(NPW,NPW) , ak2(NPW,NPW) , ak21(NPW,NPW)
      DIMENSION atr(NSYM,NPW,NPW) , ati(NSYM,NPW,NPW)
      DIMENSION akmat(NSYM,NPW,NPW)

      COMMON /INPPAR/ ipol , iprint , nt , nscat , mol ,
     &                lmax , lmaxb , lbig ,
     &                alpha0 , alpha2 , dpol , dpolex , quad , quadex ,
     &                dpolt , dpole , q1 , q2 , oct1 , oct2


      DO ns = 1 , nscat
        ks = kks(ns)
        n = ks - 1
        DO j = 1 , n
          jj = j + 1
          DO k = jj , ks
            sum = 0.0E0
            DO kk = 1 , ks
              ak(j,kk) = akmat(ns,j,kk)
              ak(kk,k) = akmat(ns,kk,k)
              sum = sum + ak(j,kk)*ak(kk,k)
            END DO
            ak2(k,j) = sum
            ak2(j,k) = sum
            ak21(j,k) = sum
            ak21(k,j) = sum
          END DO
        END DO
        DO j = 1 , ks
          sum = 0.0E0
          DO k = 1 , ks
            ak(j,k) = akmat(ns,j,k)
            a = ak(j,k)
            sum = sum + a*a
          END DO
          ak2(j,j) = sum
          ak21(j,j) = sum + 1.0E0
        END DO
c
c     Call ma01j for inverse of the matrix
c
        CALL MA01A(ak21,ks)

        DO j = 1 , ks
          DO k = 1 , ks
            smr = 0.0E0
            smi = 0.0E0
            DO kk = 1 , ks
              aki = ak21(j,kk)
              smr = smr - aki*ak2(kk,k)
              smi = smi + aki*ak(kk,k)
            END DO
            atr(ns,j,k) = smr + smr
            ati(ns,j,k) = smi + smi
          END DO
        END DO
      END DO
 
      RETURN
      END

