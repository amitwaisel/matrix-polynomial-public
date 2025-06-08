!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      SUBROUTINE RECSY_SOLVE_33(C, D, INFO, MACHINE)
!     .. Scalar Arguments ..
      IMPLICIT NONE
      INTEGER           INFO, UPLO, LDA, LDB, LDC
      DOUBLE PRECISION  C(0:2), MACHINE(0:*)
!     .. Array Arguments ..
      DOUBLE PRECISION  D(0:2,0:2)
!     .. Local Scalars ..
      INTEGER           I, PI, J, K
      DOUBLE PRECISION  P, AP, AP2, TMP, AD, AC, SMIN
!     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY
!     .. Intrinsic Functions ..
      INTRINSIC         MAX
!     .. Executable Statements ..
!
      SMIN = MACHINE(0)
      INFO = 0
      DO I = 0,2
         AP = ABS(D(I,I))
         PI = I
         DO J = I+1,2
            AP2 = ABS(D(J,I))
            IF (AP2.GT.AP) THEN
               AP = AP2
               PI = J
            END IF
         END DO
         IF (I.EQ.0) SMIN = MAX(SMIN, AP*MACHINE(2))
         IF (PI.NE.I) THEN
            DO J = I,2
               TMP = D(I,J)
               D(I,J) = D(PI,J)
               D(PI,J) = TMP
            END DO
            TMP = C(I)
            C(I) = C(PI)
            C(PI) = TMP
         END IF
         AD = AP
         IF (AD.LT.SMIN) THEN
            INFO = 1
            RETURN
         END IF
         DO J = I+1,2
            P = D(J,I) / D(I,I)
            DO K = I+1,2
               D(J,K) = D(J,K) - P * D(I,K)
            END DO
            C(J) = C(J) - P * C(I)
         END DO
      END DO
      DO I = 2,0,-1
         AD = ABS(D(I,I))
         AC = ABS(C(I))
         IF (AD.LT.1D0.AND.AC.GT.1D0) THEN
            IF (AC.GT.MACHINE(1)*AD) THEN
               INFO = 1
               RETURN
            END IF
         END IF
         C(I) = C(I) / D(I,I)
         DO J=I-1,0,-1
            C(J) = C(J) - C(I) * D(J,I)
         END DO
      END DO
      RETURN
      END
