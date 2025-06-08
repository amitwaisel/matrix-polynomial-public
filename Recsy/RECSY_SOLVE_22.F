!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      SUBROUTINE RECSY_SOLVE_22(C0, C1, D, INFO, MACHINE)
!     .. Scalar Arguments ..
      IMPLICIT NONE
      INTEGER           INFO, UPLO, LDA, LDB, LDC, M, N
      DOUBLE PRECISION  C0, C1
!     .. Array Arguments ..
      DOUBLE PRECISION  D(0:1,0:1), MACHINE(0:*)
!     .. Local Scalars ..
      INTEGER           I, PI, J, K
      DOUBLE PRECISION  P, AP, AP2, TMP, AD, AC, AD0, AD1
!     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY
!     .. Intrinsic Functions ..
      INTRINSIC         MAX
!     .. Executable Statements ..
!
      INFO = 0
      AD0 = ABS(D(0,0))
      AD1 = ABS(D(1,0))
      IF (AD0.GE.AD1) THEN
         IF (AD0.LT.MACHINE(0)) THEN
            INFO = 1
            RETURN
         END IF
         P = D(1,0) / D(0,0)
         D(1,0) = D(1,0) - D(0,0) * P
         D(1,1) = D(1,1) - D(0,1) * P
         C1 = C1 - C0 * P
         AD = ABS(D(1,1))
         AC = ABS(C1)
         IF (AD.LT.MACHINE(0).OR.AD.LT.MACHINE(2)*AD0) THEN
            INFO = 1
            RETURN
         END IF
         IF (AD.LT.1D0.AND.AC.GT.1D0) THEN
            IF (AC.GT.MACHINE(1)*AD) THEN
               INFO = 1
               RETURN
            END IF
         END IF
         C1 = C1 / D(1,1)
         C0 = C0 - D(0,1) * C1 
         AC = ABS(C0)
         IF (AD0.LT.1D0.AND.AC.GT.1D0) THEN
            IF (AC.GT.MACHINE(1)*AD0) THEN
               INFO = 1
               RETURN
            END IF
         END IF
         C0 = C0 / D(0,0)
      ELSE
         IF (AD1.LT.MACHINE(0)) THEN
            INFO = 1
            RETURN
         END IF
         P = D(0,0) / D(1,0)
         D(0,0) = D(0,0) - D(1,0) * P
         D(0,1) = D(0,1) - D(1,1) * P
         C0 = C0 - C1 * P
         AD = ABS(D(0,1))
         AC = ABS(C0)
         IF (AD.LT.MACHINE(0).OR.AD.LT.MACHINE(2)*AD1) THEN
            INFO = 1
            RETURN
         END IF
         IF (AD.LT.1D0.AND.AC.GT.1D0) THEN
            IF (AC.GT.MACHINE(1)*AD) THEN
               INFO = 1
               RETURN
            END IF
         END IF
         C0 = C0 / D(0,1)
         C1 = C1 - D(1,1) * C0 
         AC = ABS(C1)
         IF (AD1.LT.1D0.AND.AC.GT.1D0) THEN
            IF (AC.GT.MACHINE(1)*AD1) THEN
               INFO = 1
               RETURN
            END IF
         END IF
         C1 = C1 / D(1,0)
         AC = C1
         C1 = C0
         C0 = AC
      END IF
      RETURN
     END
