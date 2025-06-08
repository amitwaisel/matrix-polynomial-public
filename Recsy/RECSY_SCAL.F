!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      SUBROUTINE RECSY_SCAL(M, N, SCALE, A, LDA)
!     .. Scalar Arguments ..
      IMPLICIT NONE
      INTEGER           LDA, M, N
      DOUBLE PRECISION  SCALE
!     .. Array Arguments ..
      DOUBLE PRECISION  A(0:LDA-1,0:*)
!     .. Local Scalars ..
      INTEGER           I
!     .. External Subroutines ..
      EXTERNAL          DSCAL
!     .. Intrinsic Functions ..
!
      CALL DGEMM('N', 'N', M, N, 0, 0D0, 0D0, M, 0D0, 1, SCALE, A, LDA)
!      DO I=0,N-1
!         CALL DSCAL(M, SCALE, A(0,I), 1)
!      END DO
      RETURN
      END
