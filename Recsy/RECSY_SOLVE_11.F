!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      SUBROUTINE RECSY_SOLVE_11(C, D, INFO, MACHINE)
!     .. Scalar Arguments ..
      IMPLICIT NONE
      INTEGER           INFO
      DOUBLE PRECISION  C, D, MACHINE(0:*)
!     .. Array Arguments ..
!     .. Local Scalars ..
      DOUBLE PRECISION  AD, AC
!     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY
!     .. Intrinsic Functions ..
      INTRINSIC         MAX, ABS
!     .. Executable Statements ..
!
      INFO = 0
      AD = ABS(D)
      AC = ABS(C)
      IF (AD.LT.MACHINE(0)) THEN
         INFO = 1
         RETURN
      END IF
      IF (AD.LT.1D0.AND.AC.GT.1D0) THEN
         IF (AC.GT.MACHINE(1)*AD) THEN
            INFO = 1
            RETURN
         END IF
      END IF

      C = C / D
      RETURN
     END
