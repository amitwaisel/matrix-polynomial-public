!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      SUBROUTINE RECSY_SOLVE_CP(M, LD, C, D, SCALE, INFO, MACHINE)
!     .. Scalar Arguments ..
      IMPLICIT NONE
      INTEGER           INFO, UPLO, LDA, LDB, LDC, M, N, LD
      DOUBLE PRECISION  C(0:LD-1), MACHINE(0:*)
!     .. Array Arguments ..
      DOUBLE PRECISION  D(0:LD-1,0:LD-1)
!     .. Local Scalars ..
      INTEGER           I, II, I2, J, K, K1, K2, M1
      INTEGER           IPIV(0:M-1), JPIV(0:M-1)
      DOUBLE PRECISION  SCALE
!     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY
!     .. Intrinsic Functions ..
      INTRINSIC         MAX
!     .. Executable Statements ..
!
      CALL DGETC2(M, D, LD, IPIV, JPIV, INFO)
      CALL DGESC2(M, D, LD, C, IPIV, JPIV, SCALE)
      RETURN
      END
