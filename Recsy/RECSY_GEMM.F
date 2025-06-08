!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      RECURSIVE SUBROUTINE RECSY_GEMM(T, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, MACHINE )
      IMPLICIT NONE
!     .. Scalar Arguments ..
      INTEGER          :: TRANS, LDA, LDB, LDC, M, N, K
!     .. Array Arguments ..
      DOUBLE PRECISION :: A(0:LDA-1,0:*), B(0:LDB-1,0:*), C(0:LDC-1,0:*)
      DOUBLE PRECISION :: ALPHA, BETA, MACHINE(0:*)
!     .. Local Scalars ..
      INTEGER          :: AL, I, J, L, T
      DOUBLE PRECISION :: T00, T01, T02, T03
      DOUBLE PRECISION :: T10, T11, T12, T13
      DOUBLE PRECISION :: T20, T21, T22, T23
      DOUBLE PRECISION :: T30, T31, T32, T33
!     .. External Subroutines ..
      EXTERNAL          DGEMM
!     .. Intrinsic Functions ..
      INTRINSIC         MAX
      CHARACTER         TA, TB
!     .. Executable Statements ..
!
      trans = t
      IF (TRANS.EQ.0.OR.TRANS.EQ.2) THEN
         TA = 'N'
      ELSE
         TA = 'T'
      END IF
      IF (TRANS.EQ.0.OR.TRANS.EQ.1) THEN
         TB = 'N'
      ELSE
         TB = 'T'
      END IF

      IF (M.GE.MACHINE(3).OR.N.GE.MACHINE(3).OR.K.GE.MACHINE(3)) THEN
         CALL DGEMM(TA, TB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
         RETURN
      END IF

      trans = t
      if (trans.eq.1) then
        trans=2
      else if (trans.eq.2) then
        trans=1
      end if

      IF (ALPHA.EQ.1D0) THEN
         AL=1
      ELSEIF (ALPHA.EQ.-1D0) THEN
         AL=-1
      ELSE
         print*,'ALPHA .NE. +/- 1 not supported'
         stop 99
      END IF
      IF (BETA.EQ.0D0) THEN
         SELECT CASE (TRANS)
         CASE (0)        ! 'N','N'
            SELECT CASE (AL)
            CASE (1)
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     T33 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T03 = T03 + A(I+0,L) * B(L,J+3)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T13 = T13 + A(I+1,L) * B(L,J+3)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T22 = T22 + A(I+2,L) * B(L,J+2)
                        T23 = T23 + A(I+2,L) * B(L,J+3)
                        T30 = T30 + A(I+3,L) * B(L,J+0)
                        T31 = T31 + A(I+3,L) * B(L,J+1)
                        T32 = T32 + A(I+3,L) * B(L,J+2)
                        T33 = T33 + A(I+3,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                     C(I+3,J+3) = T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T22 = T22 + A(I+2,L) * B(L,J+2)
                        T30 = T30 + A(I+3,L) * B(L,J+0)
                        T31 = T31 + A(I+3,L) * B(L,J+1)
                        T32 = T32 + A(I+3,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T30 = T30 + A(I+3,L) * B(L,J+0)
                        T31 = T31 + A(I+3,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     T30 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T30 = T30 + A(I+3,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                     C(I+3,J+0) = T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T03 = T03 + A(I+0,L) * B(L,J+3)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T13 = T13 + A(I+1,L) * B(L,J+3)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T22 = T22 + A(I+2,L) * B(L,J+2)
                        T23 = T23 + A(I+2,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T22 = T22 + A(I+2,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T03 = T03 + A(I+0,L) * B(L,J+3)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T13 = T13 + A(I+1,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T03 = T03 + A(I+0,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                  CASE (1)
                     T00 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                  END SELECT
               END SELECT
            CASE (-1)
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     T33 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T02 = T02 - A(I+0,L) * B(L,J+2)
                        T03 = T03 - A(I+0,L) * B(L,J+3)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                        T12 = T12 - A(I+1,L) * B(L,J+2)
                        T13 = T13 - A(I+1,L) * B(L,J+3)
                        T20 = T20 - A(I+2,L) * B(L,J+0)
                        T21 = T21 - A(I+2,L) * B(L,J+1)
                        T22 = T22 - A(I+2,L) * B(L,J+2)
                        T23 = T23 - A(I+2,L) * B(L,J+3)
                        T30 = T30 - A(I+3,L) * B(L,J+0)
                        T31 = T31 - A(I+3,L) * B(L,J+1)
                        T32 = T32 - A(I+3,L) * B(L,J+2)
                        T33 = T33 - A(I+3,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                     C(I+3,J+3) = T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T02 = T02 - A(I+0,L) * B(L,J+2)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                        T12 = T12 - A(I+1,L) * B(L,J+2)
                        T20 = T20 - A(I+2,L) * B(L,J+0)
                        T21 = T21 - A(I+2,L) * B(L,J+1)
                        T22 = T22 - A(I+2,L) * B(L,J+2)
                        T30 = T30 - A(I+3,L) * B(L,J+0)
                        T31 = T31 - A(I+3,L) * B(L,J+1)
                        T32 = T32 - A(I+3,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                        T20 = T20 - A(I+2,L) * B(L,J+0)
                        T21 = T21 - A(I+2,L) * B(L,J+1)
                        T30 = T30 - A(I+3,L) * B(L,J+0)
                        T31 = T31 - A(I+3,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     T30 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T20 = T20 - A(I+2,L) * B(L,J+0)
                        T30 = T30 - A(I+3,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                     C(I+3,J+0) = T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T02 = T02 - A(I+0,L) * B(L,J+2)
                        T03 = T03 - A(I+0,L) * B(L,J+3)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                        T12 = T12 - A(I+1,L) * B(L,J+2)
                        T13 = T13 - A(I+1,L) * B(L,J+3)
                        T20 = T20 - A(I+2,L) * B(L,J+0)
                        T21 = T21 - A(I+2,L) * B(L,J+1)
                        T22 = T22 - A(I+2,L) * B(L,J+2)
                        T23 = T23 - A(I+2,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T02 = T02 - A(I+0,L) * B(L,J+2)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                        T12 = T12 - A(I+1,L) * B(L,J+2)
                        T20 = T20 - A(I+2,L) * B(L,J+0)
                        T21 = T21 - A(I+2,L) * B(L,J+1)
                        T22 = T22 - A(I+2,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                        T20 = T20 - A(I+2,L) * B(L,J+0)
                        T21 = T21 - A(I+2,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T20 = T20 - A(I+2,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T02 = T02 - A(I+0,L) * B(L,J+2)
                        T03 = T03 - A(I+0,L) * B(L,J+3)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                        T12 = T12 - A(I+1,L) * B(L,J+2)
                        T13 = T13 - A(I+1,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T02 = T02 - A(I+0,L) * B(L,J+2)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                        T12 = T12 - A(I+1,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T02 = T02 - A(I+0,L) * B(L,J+2)
                        T03 = T03 - A(I+0,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T02 = T02 - A(I+0,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                  CASE (1)
                     T00 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                  END SELECT
               END SELECT
            CASE DEFAULT
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0 
                     T13 = 0 
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     T33 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T03 = T03 + A(I+0,L) * B(L,J+3)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T13 = T13 + A(I+1,L) * B(L,J+3)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T22 = T22 + A(I+2,L) * B(L,J+2)
                        T23 = T23 + A(I+2,L) * B(L,J+3)
                        T30 = T30 + A(I+3,L) * B(L,J+0)
                        T31 = T31 + A(I+3,L) * B(L,J+1)
                        T32 = T32 + A(I+3,L) * B(L,J+2)
                        T33 = T33 + A(I+3,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+0,J+3) = 0 + ALPHA*T03
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+1,J+3) = 0 + ALPHA*T13
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+2,J+2) = 0 + ALPHA*T22
                     C(I+2,J+3) = 0 + ALPHA*T23
                     C(I+3,J+0) = 0 + ALPHA*T30
                     C(I+3,J+1) = 0 + ALPHA*T31
                     C(I+3,J+2) = 0 + ALPHA*T32
                     C(I+3,J+3) = 0 + ALPHA*T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T22 = T22 + A(I+2,L) * B(L,J+2)
                        T30 = T30 + A(I+3,L) * B(L,J+0)
                        T31 = T31 + A(I+3,L) * B(L,J+1)
                        T32 = T32 + A(I+3,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+2,J+2) = 0 + ALPHA*T22
                     C(I+3,J+0) = 0 + ALPHA*T30
                     C(I+3,J+1) = 0 + ALPHA*T31
                     C(I+3,J+2) = 0 + ALPHA*T32
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T30 = T30 + A(I+3,L) * B(L,J+0)
                        T31 = T31 + A(I+3,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+3,J+0) = 0 + ALPHA*T30
                     C(I+3,J+1) = 0 + ALPHA*T31
                  CASE (1)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T30 = T30 + A(I+3,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+3,J+0) = 0 + ALPHA*T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T03 = T03 + A(I+0,L) * B(L,J+3)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T13 = T13 + A(I+1,L) * B(L,J+3)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T22 = T22 + A(I+2,L) * B(L,J+2)
                        T23 = T23 + A(I+2,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+0,J+3) = 0 + ALPHA*T03
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+1,J+3) = 0 + ALPHA*T13
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+2,J+2) = 0 + ALPHA*T22
                     C(I+2,J+3) = 0 + ALPHA*T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T22 = T22 + A(I+2,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+2,J+2) = 0 + ALPHA*T22
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+2,J+0) = 0 + ALPHA*T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T03 = T03 + A(I+0,L) * B(L,J+3)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T13 = T13 + A(I+1,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+0,J+3) = 0 + ALPHA*T03
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+1,J+3) = 0 + ALPHA*T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+1,J+0) = 0 + ALPHA*T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T03 = T03 + A(I+0,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+0,J+3) = 0 + ALPHA*T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                  CASE (1)
                     T00 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                  END SELECT
               END SELECT
            END SELECT
         CASE (1)        ! 'N','T'
            SELECT CASE (AL)
            CASE (1)
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     T33 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T03 = T03 + A(I+0,L) * B(J+3,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T13 = T13 + A(I+1,L) * B(J+3,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T22 = T22 + A(I+2,L) * B(J+2,L)
                        T23 = T23 + A(I+2,L) * B(J+3,L)
                        T30 = T30 + A(I+3,L) * B(J+0,L)
                        T31 = T31 + A(I+3,L) * B(J+1,L)
                        T32 = T32 + A(I+3,L) * B(J+2,L)
                        T33 = T33 + A(I+3,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                     C(I+3,J+3) = T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T22 = T22 + A(I+2,L) * B(J+2,L)
                        T30 = T30 + A(I+3,L) * B(J+0,L)
                        T31 = T31 + A(I+3,L) * B(J+1,L)
                        T32 = T32 + A(I+3,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T30 = T30 + A(I+3,L) * B(J+0,L)
                        T31 = T31 + A(I+3,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     T30 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T30 = T30 + A(I+3,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                     C(I+3,J+0) = T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T03 = T03 + A(I+0,L) * B(J+3,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T13 = T13 + A(I+1,L) * B(J+3,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T22 = T22 + A(I+2,L) * B(J+2,L)
                        T23 = T23 + A(I+2,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T22 = T22 + A(I+2,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T03 = T03 + A(I+0,L) * B(J+3,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T13 = T13 + A(I+1,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T03 = T03 + A(I+0,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                  CASE (1)
                     T00 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                  END SELECT
               END SELECT
            CASE (-1)
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     T33 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T02 = T02 - A(I+0,L) * B(J+2,L)
                        T03 = T03 - A(I+0,L) * B(J+3,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                        T12 = T12 - A(I+1,L) * B(J+2,L)
                        T13 = T13 - A(I+1,L) * B(J+3,L)
                        T20 = T20 - A(I+2,L) * B(J+0,L)
                        T21 = T21 - A(I+2,L) * B(J+1,L)
                        T22 = T22 - A(I+2,L) * B(J+2,L)
                        T23 = T23 - A(I+2,L) * B(J+3,L)
                        T30 = T30 - A(I+3,L) * B(J+0,L)
                        T31 = T31 - A(I+3,L) * B(J+1,L)
                        T32 = T32 - A(I+3,L) * B(J+2,L)
                        T33 = T33 - A(I+3,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                     C(I+3,J+3) = T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T02 = T02 - A(I+0,L) * B(J+2,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                        T12 = T12 - A(I+1,L) * B(J+2,L)
                        T20 = T20 - A(I+2,L) * B(J+0,L)
                        T21 = T21 - A(I+2,L) * B(J+1,L)
                        T22 = T22 - A(I+2,L) * B(J+2,L)
                        T30 = T30 - A(I+3,L) * B(J+0,L)
                        T31 = T31 - A(I+3,L) * B(J+1,L)
                        T32 = T32 - A(I+3,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                        T20 = T20 - A(I+2,L) * B(J+0,L)
                        T21 = T21 - A(I+2,L) * B(J+1,L)
                        T30 = T30 - A(I+3,L) * B(J+0,L)
                        T31 = T31 - A(I+3,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     T30 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T20 = T20 - A(I+2,L) * B(J+0,L)
                        T30 = T30 - A(I+3,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                     C(I+3,J+0) = T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T02 = T02 - A(I+0,L) * B(J+2,L)
                        T03 = T03 - A(I+0,L) * B(J+3,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                        T12 = T12 - A(I+1,L) * B(J+2,L)
                        T13 = T13 - A(I+1,L) * B(J+3,L)
                        T20 = T20 - A(I+2,L) * B(J+0,L)
                        T21 = T21 - A(I+2,L) * B(J+1,L)
                        T22 = T22 - A(I+2,L) * B(J+2,L)
                        T23 = T23 - A(I+2,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T02 = T02 - A(I+0,L) * B(J+2,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                        T12 = T12 - A(I+1,L) * B(J+2,L)
                        T20 = T20 - A(I+2,L) * B(J+0,L)
                        T21 = T21 - A(I+2,L) * B(J+1,L)
                        T22 = T22 - A(I+2,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                        T20 = T20 - A(I+2,L) * B(J+0,L)
                        T21 = T21 - A(I+2,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T20 = T20 - A(I+2,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T02 = T02 - A(I+0,L) * B(J+2,L)
                        T03 = T03 - A(I+0,L) * B(J+3,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                        T12 = T12 - A(I+1,L) * B(J+2,L)
                        T13 = T13 - A(I+1,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T02 = T02 - A(I+0,L) * B(J+2,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                        T12 = T12 - A(I+1,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T02 = T02 - A(I+0,L) * B(J+2,L)
                        T03 = T03 - A(I+0,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T02 = T02 - A(I+0,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                  CASE (1)
                     T00 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                  END SELECT
               END SELECT
            CASE DEFAULT
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0 
                     T13 = 0 
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     T33 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T03 = T03 + A(I+0,L) * B(J+3,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T13 = T13 + A(I+1,L) * B(J+3,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T22 = T22 + A(I+2,L) * B(J+2,L)
                        T23 = T23 + A(I+2,L) * B(J+3,L)
                        T30 = T30 + A(I+3,L) * B(J+0,L)
                        T31 = T31 + A(I+3,L) * B(J+1,L)
                        T32 = T32 + A(I+3,L) * B(J+2,L)
                        T33 = T33 + A(I+3,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+0,J+3) = 0 + ALPHA*T03
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+1,J+3) = 0 + ALPHA*T13
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+2,J+2) = 0 + ALPHA*T22
                     C(I+2,J+3) = 0 + ALPHA*T23
                     C(I+3,J+0) = 0 + ALPHA*T30
                     C(I+3,J+1) = 0 + ALPHA*T31
                     C(I+3,J+2) = 0 + ALPHA*T32
                     C(I+3,J+3) = 0 + ALPHA*T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T22 = T22 + A(I+2,L) * B(J+2,L)
                        T30 = T30 + A(I+3,L) * B(J+0,L)
                        T31 = T31 + A(I+3,L) * B(J+1,L)
                        T32 = T32 + A(I+3,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+2,J+2) = 0 + ALPHA*T22
                     C(I+3,J+0) = 0 + ALPHA*T30
                     C(I+3,J+1) = 0 + ALPHA*T31
                     C(I+3,J+2) = 0 + ALPHA*T32
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T30 = T30 + A(I+3,L) * B(J+0,L)
                        T31 = T31 + A(I+3,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+3,J+0) = 0 + ALPHA*T30
                     C(I+3,J+1) = 0 + ALPHA*T31
                  CASE (1)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T30 = T30 + A(I+3,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+3,J+0) = 0 + ALPHA*T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T03 = T03 + A(I+0,L) * B(J+3,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T13 = T13 + A(I+1,L) * B(J+3,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T22 = T22 + A(I+2,L) * B(J+2,L)
                        T23 = T23 + A(I+2,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+0,J+3) = 0 + ALPHA*T03
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+1,J+3) = 0 + ALPHA*T13
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+2,J+2) = 0 + ALPHA*T22
                     C(I+2,J+3) = 0 + ALPHA*T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T22 = T22 + A(I+2,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+2,J+2) = 0 + ALPHA*T22
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+2,J+0) = 0 + ALPHA*T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T03 = T03 + A(I+0,L) * B(J+3,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T13 = T13 + A(I+1,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+0,J+3) = 0 + ALPHA*T03
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+1,J+3) = 0 + ALPHA*T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+1,J+0) = 0 + ALPHA*T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T03 = T03 + A(I+0,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+0,J+3) = 0 + ALPHA*T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                  CASE (1)
                     T00 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                  END SELECT
               END SELECT
            END SELECT
         CASE (2)        ! 'T','N'
            SELECT CASE (AL)
            CASE (1)
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     T33 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T03 = T03 + A(L,I+0) * B(L,J+3)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T13 = T13 + A(L,I+1) * B(L,J+3)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T22 = T22 + A(L,I+2) * B(L,J+2)
                        T23 = T23 + A(L,I+2) * B(L,J+3)
                        T30 = T30 + A(L,I+3) * B(L,J+0)
                        T31 = T31 + A(L,I+3) * B(L,J+1)
                        T32 = T32 + A(L,I+3) * B(L,J+2)
                        T33 = T33 + A(L,I+3) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                     C(I+3,J+3) = T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T22 = T22 + A(L,I+2) * B(L,J+2)
                        T30 = T30 + A(L,I+3) * B(L,J+0)
                        T31 = T31 + A(L,I+3) * B(L,J+1)
                        T32 = T32 + A(L,I+3) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T30 = T30 + A(L,I+3) * B(L,J+0)
                        T31 = T31 + A(L,I+3) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     T30 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T30 = T30 + A(L,I+3) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                     C(I+3,J+0) = T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T03 = T03 + A(L,I+0) * B(L,J+3)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T13 = T13 + A(L,I+1) * B(L,J+3)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T22 = T22 + A(L,I+2) * B(L,J+2)
                        T23 = T23 + A(L,I+2) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T22 = T22 + A(L,I+2) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T03 = T03 + A(L,I+0) * B(L,J+3)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T13 = T13 + A(L,I+1) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T03 = T03 + A(L,I+0) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                  CASE (1)
                     T00 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                  END SELECT
               END SELECT
            CASE (-1)
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     T33 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T02 = T02 - A(L,I+0) * B(L,J+2)
                        T03 = T03 - A(L,I+0) * B(L,J+3)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                        T12 = T12 - A(L,I+1) * B(L,J+2)
                        T13 = T13 - A(L,I+1) * B(L,J+3)
                        T20 = T20 - A(L,I+2) * B(L,J+0)
                        T21 = T21 - A(L,I+2) * B(L,J+1)
                        T22 = T22 - A(L,I+2) * B(L,J+2)
                        T23 = T23 - A(L,I+2) * B(L,J+3)
                        T30 = T30 - A(L,I+3) * B(L,J+0)
                        T31 = T31 - A(L,I+3) * B(L,J+1)
                        T32 = T32 - A(L,I+3) * B(L,J+2)
                        T33 = T33 - A(L,I+3) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                     C(I+3,J+3) = T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T02 = T02 - A(L,I+0) * B(L,J+2)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                        T12 = T12 - A(L,I+1) * B(L,J+2)
                        T20 = T20 - A(L,I+2) * B(L,J+0)
                        T21 = T21 - A(L,I+2) * B(L,J+1)
                        T22 = T22 - A(L,I+2) * B(L,J+2)
                        T30 = T30 - A(L,I+3) * B(L,J+0)
                        T31 = T31 - A(L,I+3) * B(L,J+1)
                        T32 = T32 - A(L,I+3) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                        T20 = T20 - A(L,I+2) * B(L,J+0)
                        T21 = T21 - A(L,I+2) * B(L,J+1)
                        T30 = T30 - A(L,I+3) * B(L,J+0)
                        T31 = T31 - A(L,I+3) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     T30 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T20 = T20 - A(L,I+2) * B(L,J+0)
                        T30 = T30 - A(L,I+3) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                     C(I+3,J+0) = T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T02 = T02 - A(L,I+0) * B(L,J+2)
                        T03 = T03 - A(L,I+0) * B(L,J+3)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                        T12 = T12 - A(L,I+1) * B(L,J+2)
                        T13 = T13 - A(L,I+1) * B(L,J+3)
                        T20 = T20 - A(L,I+2) * B(L,J+0)
                        T21 = T21 - A(L,I+2) * B(L,J+1)
                        T22 = T22 - A(L,I+2) * B(L,J+2)
                        T23 = T23 - A(L,I+2) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T02 = T02 - A(L,I+0) * B(L,J+2)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                        T12 = T12 - A(L,I+1) * B(L,J+2)
                        T20 = T20 - A(L,I+2) * B(L,J+0)
                        T21 = T21 - A(L,I+2) * B(L,J+1)
                        T22 = T22 - A(L,I+2) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                        T20 = T20 - A(L,I+2) * B(L,J+0)
                        T21 = T21 - A(L,I+2) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T20 = T20 - A(L,I+2) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T02 = T02 - A(L,I+0) * B(L,J+2)
                        T03 = T03 - A(L,I+0) * B(L,J+3)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                        T12 = T12 - A(L,I+1) * B(L,J+2)
                        T13 = T13 - A(L,I+1) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T02 = T02 - A(L,I+0) * B(L,J+2)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                        T12 = T12 - A(L,I+1) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T02 = T02 - A(L,I+0) * B(L,J+2)
                        T03 = T03 - A(L,I+0) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T02 = T02 - A(L,I+0) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                  CASE (1)
                     T00 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                  END SELECT
               END SELECT
            CASE DEFAULT
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0 
                     T13 = 0 
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     T33 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T03 = T03 + A(L,I+0) * B(L,J+3)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T13 = T13 + A(L,I+1) * B(L,J+3)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T22 = T22 + A(L,I+2) * B(L,J+2)
                        T23 = T23 + A(L,I+2) * B(L,J+3)
                        T30 = T30 + A(L,I+3) * B(L,J+0)
                        T31 = T31 + A(L,I+3) * B(L,J+1)
                        T32 = T32 + A(L,I+3) * B(L,J+2)
                        T33 = T33 + A(L,I+3) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+0,J+3) = 0 + ALPHA*T03
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+1,J+3) = 0 + ALPHA*T13
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+2,J+2) = 0 + ALPHA*T22
                     C(I+2,J+3) = 0 + ALPHA*T23
                     C(I+3,J+0) = 0 + ALPHA*T30
                     C(I+3,J+1) = 0 + ALPHA*T31
                     C(I+3,J+2) = 0 + ALPHA*T32
                     C(I+3,J+3) = 0 + ALPHA*T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T22 = T22 + A(L,I+2) * B(L,J+2)
                        T30 = T30 + A(L,I+3) * B(L,J+0)
                        T31 = T31 + A(L,I+3) * B(L,J+1)
                        T32 = T32 + A(L,I+3) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+2,J+2) = 0 + ALPHA*T22
                     C(I+3,J+0) = 0 + ALPHA*T30
                     C(I+3,J+1) = 0 + ALPHA*T31
                     C(I+3,J+2) = 0 + ALPHA*T32
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T30 = T30 + A(L,I+3) * B(L,J+0)
                        T31 = T31 + A(L,I+3) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+3,J+0) = 0 + ALPHA*T30
                     C(I+3,J+1) = 0 + ALPHA*T31
                  CASE (1)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T30 = T30 + A(L,I+3) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+3,J+0) = 0 + ALPHA*T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T03 = T03 + A(L,I+0) * B(L,J+3)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T13 = T13 + A(L,I+1) * B(L,J+3)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T22 = T22 + A(L,I+2) * B(L,J+2)
                        T23 = T23 + A(L,I+2) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+0,J+3) = 0 + ALPHA*T03
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+1,J+3) = 0 + ALPHA*T13
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+2,J+2) = 0 + ALPHA*T22
                     C(I+2,J+3) = 0 + ALPHA*T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T22 = T22 + A(L,I+2) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+2,J+2) = 0 + ALPHA*T22
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+2,J+0) = 0 + ALPHA*T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T03 = T03 + A(L,I+0) * B(L,J+3)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T13 = T13 + A(L,I+1) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+0,J+3) = 0 + ALPHA*T03
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+1,J+3) = 0 + ALPHA*T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+1,J+0) = 0 + ALPHA*T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T03 = T03 + A(L,I+0) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+0,J+3) = 0 + ALPHA*T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                  CASE (1)
                     T00 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                  END SELECT
               END SELECT
            END SELECT
         CASE (3)        ! 'T','T'
            SELECT CASE (AL)
            CASE (1)
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     T33 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T03 = T03 + A(L,I+0) * B(J+3,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T13 = T13 + A(L,I+1) * B(J+3,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T22 = T22 + A(L,I+2) * B(J+2,L)
                        T23 = T23 + A(L,I+2) * B(J+3,L)
                        T30 = T30 + A(L,I+3) * B(J+0,L)
                        T31 = T31 + A(L,I+3) * B(J+1,L)
                        T32 = T32 + A(L,I+3) * B(J+2,L)
                        T33 = T33 + A(L,I+3) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                     C(I+3,J+3) = T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T22 = T22 + A(L,I+2) * B(J+2,L)
                        T30 = T30 + A(L,I+3) * B(J+0,L)
                        T31 = T31 + A(L,I+3) * B(J+1,L)
                        T32 = T32 + A(L,I+3) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T30 = T30 + A(L,I+3) * B(J+0,L)
                        T31 = T31 + A(L,I+3) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     T30 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T30 = T30 + A(L,I+3) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                     C(I+3,J+0) = T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T03 = T03 + A(L,I+0) * B(J+3,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T13 = T13 + A(L,I+1) * B(J+3,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T22 = T22 + A(L,I+2) * B(J+2,L)
                        T23 = T23 + A(L,I+2) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T22 = T22 + A(L,I+2) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T03 = T03 + A(L,I+0) * B(J+3,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T13 = T13 + A(L,I+1) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T03 = T03 + A(L,I+0) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                  CASE (1)
                     T00 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                  END SELECT
               END SELECT
            CASE (-1)
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     T33 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T02 = T02 - A(L,I+0) * B(J+2,L)
                        T03 = T03 - A(L,I+0) * B(J+3,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                        T12 = T12 - A(L,I+1) * B(J+2,L)
                        T13 = T13 - A(L,I+1) * B(J+3,L)
                        T20 = T20 - A(L,I+2) * B(J+0,L)
                        T21 = T21 - A(L,I+2) * B(J+1,L)
                        T22 = T22 - A(L,I+2) * B(J+2,L)
                        T23 = T23 - A(L,I+2) * B(J+3,L)
                        T30 = T30 - A(L,I+3) * B(J+0,L)
                        T31 = T31 - A(L,I+3) * B(J+1,L)
                        T32 = T32 - A(L,I+3) * B(J+2,L)
                        T33 = T33 - A(L,I+3) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                     C(I+3,J+3) = T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T02 = T02 - A(L,I+0) * B(J+2,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                        T12 = T12 - A(L,I+1) * B(J+2,L)
                        T20 = T20 - A(L,I+2) * B(J+0,L)
                        T21 = T21 - A(L,I+2) * B(J+1,L)
                        T22 = T22 - A(L,I+2) * B(J+2,L)
                        T30 = T30 - A(L,I+3) * B(J+0,L)
                        T31 = T31 - A(L,I+3) * B(J+1,L)
                        T32 = T32 - A(L,I+3) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                        T20 = T20 - A(L,I+2) * B(J+0,L)
                        T21 = T21 - A(L,I+2) * B(J+1,L)
                        T30 = T30 - A(L,I+3) * B(J+0,L)
                        T31 = T31 - A(L,I+3) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     T30 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T20 = T20 - A(L,I+2) * B(J+0,L)
                        T30 = T30 - A(L,I+3) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                     C(I+3,J+0) = T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T02 = T02 - A(L,I+0) * B(J+2,L)
                        T03 = T03 - A(L,I+0) * B(J+3,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                        T12 = T12 - A(L,I+1) * B(J+2,L)
                        T13 = T13 - A(L,I+1) * B(J+3,L)
                        T20 = T20 - A(L,I+2) * B(J+0,L)
                        T21 = T21 - A(L,I+2) * B(J+1,L)
                        T22 = T22 - A(L,I+2) * B(J+2,L)
                        T23 = T23 - A(L,I+2) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T02 = T02 - A(L,I+0) * B(J+2,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                        T12 = T12 - A(L,I+1) * B(J+2,L)
                        T20 = T20 - A(L,I+2) * B(J+0,L)
                        T21 = T21 - A(L,I+2) * B(J+1,L)
                        T22 = T22 - A(L,I+2) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                        T20 = T20 - A(L,I+2) * B(J+0,L)
                        T21 = T21 - A(L,I+2) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T20 = T20 - A(L,I+2) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T02 = T02 - A(L,I+0) * B(J+2,L)
                        T03 = T03 - A(L,I+0) * B(J+3,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                        T12 = T12 - A(L,I+1) * B(J+2,L)
                        T13 = T13 - A(L,I+1) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T02 = T02 - A(L,I+0) * B(J+2,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                        T12 = T12 - A(L,I+1) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T02 = T02 - A(L,I+0) * B(J+2,L)
                        T03 = T03 - A(L,I+0) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T02 = T02 - A(L,I+0) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                  CASE (1)
                     T00 = 0
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                  END SELECT
               END SELECT
            CASE DEFAULT
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0 
                     T13 = 0 
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     T33 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T03 = T03 + A(L,I+0) * B(J+3,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T13 = T13 + A(L,I+1) * B(J+3,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T22 = T22 + A(L,I+2) * B(J+2,L)
                        T23 = T23 + A(L,I+2) * B(J+3,L)
                        T30 = T30 + A(L,I+3) * B(J+0,L)
                        T31 = T31 + A(L,I+3) * B(J+1,L)
                        T32 = T32 + A(L,I+3) * B(J+2,L)
                        T33 = T33 + A(L,I+3) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+0,J+3) = 0 + ALPHA*T03
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+1,J+3) = 0 + ALPHA*T13
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+2,J+2) = 0 + ALPHA*T22
                     C(I+2,J+3) = 0 + ALPHA*T23
                     C(I+3,J+0) = 0 + ALPHA*T30
                     C(I+3,J+1) = 0 + ALPHA*T31
                     C(I+3,J+2) = 0 + ALPHA*T32
                     C(I+3,J+3) = 0 + ALPHA*T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T22 = T22 + A(L,I+2) * B(J+2,L)
                        T30 = T30 + A(L,I+3) * B(J+0,L)
                        T31 = T31 + A(L,I+3) * B(J+1,L)
                        T32 = T32 + A(L,I+3) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+2,J+2) = 0 + ALPHA*T22
                     C(I+3,J+0) = 0 + ALPHA*T30
                     C(I+3,J+1) = 0 + ALPHA*T31
                     C(I+3,J+2) = 0 + ALPHA*T32
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T30 = T30 + A(L,I+3) * B(J+0,L)
                        T31 = T31 + A(L,I+3) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+3,J+0) = 0 + ALPHA*T30
                     C(I+3,J+1) = 0 + ALPHA*T31
                  CASE (1)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T30 = T30 + A(L,I+3) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+3,J+0) = 0 + ALPHA*T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T03 = T03 + A(L,I+0) * B(J+3,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T13 = T13 + A(L,I+1) * B(J+3,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T22 = T22 + A(L,I+2) * B(J+2,L)
                        T23 = T23 + A(L,I+2) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+0,J+3) = 0 + ALPHA*T03
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+1,J+3) = 0 + ALPHA*T13
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+2,J+2) = 0 + ALPHA*T22
                     C(I+2,J+3) = 0 + ALPHA*T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T22 = T22 + A(L,I+2) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                     C(I+2,J+2) = 0 + ALPHA*T22
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+2,J+0) = 0 + ALPHA*T20
                     C(I+2,J+1) = 0 + ALPHA*T21
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+2,J+0) = 0 + ALPHA*T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T03 = T03 + A(L,I+0) * B(J+3,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T13 = T13 + A(L,I+1) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+0,J+3) = 0 + ALPHA*T03
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                     C(I+1,J+3) = 0 + ALPHA*T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                     C(I+1,J+2) = 0 + ALPHA*T12
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+1,J+0) = 0 + ALPHA*T10
                     C(I+1,J+1) = 0 + ALPHA*T11
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+1,J+0) = 0 + ALPHA*T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T03 = T03 + A(L,I+0) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                     C(I+0,J+3) = 0 + ALPHA*T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                     C(I+0,J+2) = 0 + ALPHA*T02
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                     C(I+0,J+1) = 0 + ALPHA*T01
                  CASE (1)
                     T00 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = 0 + ALPHA*T00
                  END SELECT
               END SELECT
            END SELECT
         END SELECT
      ELSE
         SELECT CASE (TRANS)
         CASE (0)        ! 'N','N'
            SELECT CASE (AL)
            CASE (1)
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T23 = BETA*C(I+2,J+3)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     T32 = BETA*C(I+3,J+2)
                     T33 = BETA*C(I+3,J+3)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T03 = T03 + A(I+0,L) * B(L,J+3)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T13 = T13 + A(I+1,L) * B(L,J+3)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T22 = T22 + A(I+2,L) * B(L,J+2)
                        T23 = T23 + A(I+2,L) * B(L,J+3)
                        T30 = T30 + A(I+3,L) * B(L,J+0)
                        T31 = T31 + A(I+3,L) * B(L,J+1)
                        T32 = T32 + A(I+3,L) * B(L,J+2)
                        T33 = T33 + A(I+3,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                     C(I+3,J+3) = T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     T32 = BETA*C(I+3,J+2)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T22 = T22 + A(I+2,L) * B(L,J+2)
                        T30 = T30 + A(I+3,L) * B(L,J+0)
                        T31 = T31 + A(I+3,L) * B(L,J+1)
                        T32 = T32 + A(I+3,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T30 = T30 + A(I+3,L) * B(L,J+0)
                        T31 = T31 + A(I+3,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     T20 = BETA*C(I+2,J+0)
                     T30 = BETA*C(I+3,J+0)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T30 = T30 + A(I+3,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                     C(I+3,J+0) = T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T23 = BETA*C(I+2,J+3)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T03 = T03 + A(I+0,L) * B(L,J+3)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T13 = T13 + A(I+1,L) * B(L,J+3)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T22 = T22 + A(I+2,L) * B(L,J+2)
                        T23 = T23 + A(I+2,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T22 = T22 + A(I+2,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     T20 = BETA*C(I+2,J+0)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T03 = T03 + A(I+0,L) * B(L,J+3)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T13 = T13 + A(I+1,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T03 = T03 + A(I+0,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                  END SELECT
               END SELECT
            CASE (-1)
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T23 = BETA*C(I+2,J+3)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     T32 = BETA*C(I+3,J+2)
                     T33 = BETA*C(I+3,J+3)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T02 = T02 - A(I+0,L) * B(L,J+2)
                        T03 = T03 - A(I+0,L) * B(L,J+3)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                        T12 = T12 - A(I+1,L) * B(L,J+2)
                        T13 = T13 - A(I+1,L) * B(L,J+3)
                        T20 = T20 - A(I+2,L) * B(L,J+0)
                        T21 = T21 - A(I+2,L) * B(L,J+1)
                        T22 = T22 - A(I+2,L) * B(L,J+2)
                        T23 = T23 - A(I+2,L) * B(L,J+3)
                        T30 = T30 - A(I+3,L) * B(L,J+0)
                        T31 = T31 - A(I+3,L) * B(L,J+1)
                        T32 = T32 - A(I+3,L) * B(L,J+2)
                        T33 = T33 - A(I+3,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                     C(I+3,J+3) = T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     T32 = BETA*C(I+3,J+2)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T02 = T02 - A(I+0,L) * B(L,J+2)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                        T12 = T12 - A(I+1,L) * B(L,J+2)
                        T20 = T20 - A(I+2,L) * B(L,J+0)
                        T21 = T21 - A(I+2,L) * B(L,J+1)
                        T22 = T22 - A(I+2,L) * B(L,J+2)
                        T30 = T30 - A(I+3,L) * B(L,J+0)
                        T31 = T31 - A(I+3,L) * B(L,J+1)
                        T32 = T32 - A(I+3,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                        T20 = T20 - A(I+2,L) * B(L,J+0)
                        T21 = T21 - A(I+2,L) * B(L,J+1)
                        T30 = T30 - A(I+3,L) * B(L,J+0)
                        T31 = T31 - A(I+3,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     T20 = BETA*C(I+2,J+0)
                     T30 = BETA*C(I+3,J+0)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T20 = T20 - A(I+2,L) * B(L,J+0)
                        T30 = T30 - A(I+3,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                     C(I+3,J+0) = T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T23 = BETA*C(I+2,J+3)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T02 = T02 - A(I+0,L) * B(L,J+2)
                        T03 = T03 - A(I+0,L) * B(L,J+3)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                        T12 = T12 - A(I+1,L) * B(L,J+2)
                        T13 = T13 - A(I+1,L) * B(L,J+3)
                        T20 = T20 - A(I+2,L) * B(L,J+0)
                        T21 = T21 - A(I+2,L) * B(L,J+1)
                        T22 = T22 - A(I+2,L) * B(L,J+2)
                        T23 = T23 - A(I+2,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T02 = T02 - A(I+0,L) * B(L,J+2)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                        T12 = T12 - A(I+1,L) * B(L,J+2)
                        T20 = T20 - A(I+2,L) * B(L,J+0)
                        T21 = T21 - A(I+2,L) * B(L,J+1)
                        T22 = T22 - A(I+2,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                        T20 = T20 - A(I+2,L) * B(L,J+0)
                        T21 = T21 - A(I+2,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     T20 = BETA*C(I+2,J+0)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T20 = T20 - A(I+2,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T02 = T02 - A(I+0,L) * B(L,J+2)
                        T03 = T03 - A(I+0,L) * B(L,J+3)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                        T12 = T12 - A(I+1,L) * B(L,J+2)
                        T13 = T13 - A(I+1,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T02 = T02 - A(I+0,L) * B(L,J+2)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                        T12 = T12 - A(I+1,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                        T11 = T11 - A(I+1,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T10 = T10 - A(I+1,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T02 = T02 - A(I+0,L) * B(L,J+2)
                        T03 = T03 - A(I+0,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                        T02 = T02 - A(I+0,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                        T01 = T01 - A(I+0,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                  END SELECT
               END SELECT
            CASE DEFAULT
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0 
                     T13 = 0 
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     T33 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T03 = T03 + A(I+0,L) * B(L,J+3)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T13 = T13 + A(I+1,L) * B(L,J+3)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T22 = T22 + A(I+2,L) * B(L,J+2)
                        T23 = T23 + A(I+2,L) * B(L,J+3)
                        T30 = T30 + A(I+3,L) * B(L,J+0)
                        T31 = T31 + A(I+3,L) * B(L,J+1)
                        T32 = T32 + A(I+3,L) * B(L,J+2)
                        T33 = T33 + A(I+3,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+0,J+3) = BETA*C(I+0,J+3) + ALPHA*T03
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+1,J+3) = BETA*C(I+1,J+3) + ALPHA*T13
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+2,J+2) = BETA*C(I+2,J+2) + ALPHA*T22
                     C(I+2,J+3) = BETA*C(I+2,J+3) + ALPHA*T23
                     C(I+3,J+0) = BETA*C(I+3,J+0) + ALPHA*T30
                     C(I+3,J+1) = BETA*C(I+3,J+1) + ALPHA*T31
                     C(I+3,J+2) = BETA*C(I+3,J+2) + ALPHA*T32
                     C(I+3,J+3) = BETA*C(I+3,J+3) + ALPHA*T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T22 = T22 + A(I+2,L) * B(L,J+2)
                        T30 = T30 + A(I+3,L) * B(L,J+0)
                        T31 = T31 + A(I+3,L) * B(L,J+1)
                        T32 = T32 + A(I+3,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+2,J+2) = BETA*C(I+2,J+2) + ALPHA*T22
                     C(I+3,J+0) = BETA*C(I+3,J+0) + ALPHA*T30
                     C(I+3,J+1) = BETA*C(I+3,J+1) + ALPHA*T31
                     C(I+3,J+2) = BETA*C(I+3,J+2) + ALPHA*T32
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T30 = T30 + A(I+3,L) * B(L,J+0)
                        T31 = T31 + A(I+3,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+3,J+0) = BETA*C(I+3,J+0) + ALPHA*T30
                     C(I+3,J+1) = BETA*C(I+3,J+1) + ALPHA*T31
                  CASE (1)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T30 = T30 + A(I+3,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+3,J+0) = BETA*C(I+3,J+0) + ALPHA*T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T03 = T03 + A(I+0,L) * B(L,J+3)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T13 = T13 + A(I+1,L) * B(L,J+3)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T22 = T22 + A(I+2,L) * B(L,J+2)
                        T23 = T23 + A(I+2,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+0,J+3) = BETA*C(I+0,J+3) + ALPHA*T03
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+1,J+3) = BETA*C(I+1,J+3) + ALPHA*T13
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+2,J+2) = BETA*C(I+2,J+2) + ALPHA*T22
                     C(I+2,J+3) = BETA*C(I+2,J+3) + ALPHA*T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                        T22 = T22 + A(I+2,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+2,J+2) = BETA*C(I+2,J+2) + ALPHA*T22
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                        T21 = T21 + A(I+2,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T20 = T20 + A(I+2,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T03 = T03 + A(I+0,L) * B(L,J+3)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                        T13 = T13 + A(I+1,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+0,J+3) = BETA*C(I+0,J+3) + ALPHA*T03
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+1,J+3) = BETA*C(I+1,J+3) + ALPHA*T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                        T12 = T12 + A(I+1,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                        T11 = T11 + A(I+1,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T10 = T10 + A(I+1,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                        T03 = T03 + A(I+0,L) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+0,J+3) = BETA*C(I+0,J+3) + ALPHA*T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                        T02 = T02 + A(I+0,L) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                        T01 = T01 + A(I+0,L) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                  CASE (1)
                     T00 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                  END SELECT
               END SELECT
            END SELECT
         CASE (1)        ! 'N','T'
            SELECT CASE (AL)
            CASE (1)
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T23 = BETA*C(I+2,J+3)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     T32 = BETA*C(I+3,J+2)
                     T33 = BETA*C(I+3,J+3)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T03 = T03 + A(I+0,L) * B(J+3,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T13 = T13 + A(I+1,L) * B(J+3,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T22 = T22 + A(I+2,L) * B(J+2,L)
                        T23 = T23 + A(I+2,L) * B(J+3,L)
                        T30 = T30 + A(I+3,L) * B(J+0,L)
                        T31 = T31 + A(I+3,L) * B(J+1,L)
                        T32 = T32 + A(I+3,L) * B(J+2,L)
                        T33 = T33 + A(I+3,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                     C(I+3,J+3) = T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     T32 = BETA*C(I+3,J+2)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T22 = T22 + A(I+2,L) * B(J+2,L)
                        T30 = T30 + A(I+3,L) * B(J+0,L)
                        T31 = T31 + A(I+3,L) * B(J+1,L)
                        T32 = T32 + A(I+3,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T30 = T30 + A(I+3,L) * B(J+0,L)
                        T31 = T31 + A(I+3,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     T20 = BETA*C(I+2,J+0)
                     T30 = BETA*C(I+3,J+0)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T30 = T30 + A(I+3,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                     C(I+3,J+0) = T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T23 = BETA*C(I+2,J+3)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T03 = T03 + A(I+0,L) * B(J+3,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T13 = T13 + A(I+1,L) * B(J+3,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T22 = T22 + A(I+2,L) * B(J+2,L)
                        T23 = T23 + A(I+2,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T22 = T22 + A(I+2,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     T20 = BETA*C(I+2,J+0)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T03 = T03 + A(I+0,L) * B(J+3,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T13 = T13 + A(I+1,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T03 = T03 + A(I+0,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                  END SELECT
               END SELECT
            CASE (-1)
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T23 = BETA*C(I+2,J+3)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     T32 = BETA*C(I+3,J+2)
                     T33 = BETA*C(I+3,J+3)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T02 = T02 - A(I+0,L) * B(J+2,L)
                        T03 = T03 - A(I+0,L) * B(J+3,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                        T12 = T12 - A(I+1,L) * B(J+2,L)
                        T13 = T13 - A(I+1,L) * B(J+3,L)
                        T20 = T20 - A(I+2,L) * B(J+0,L)
                        T21 = T21 - A(I+2,L) * B(J+1,L)
                        T22 = T22 - A(I+2,L) * B(J+2,L)
                        T23 = T23 - A(I+2,L) * B(J+3,L)
                        T30 = T30 - A(I+3,L) * B(J+0,L)
                        T31 = T31 - A(I+3,L) * B(J+1,L)
                        T32 = T32 - A(I+3,L) * B(J+2,L)
                        T33 = T33 - A(I+3,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                     C(I+3,J+3) = T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     T32 = BETA*C(I+3,J+2)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T02 = T02 - A(I+0,L) * B(J+2,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                        T12 = T12 - A(I+1,L) * B(J+2,L)
                        T20 = T20 - A(I+2,L) * B(J+0,L)
                        T21 = T21 - A(I+2,L) * B(J+1,L)
                        T22 = T22 - A(I+2,L) * B(J+2,L)
                        T30 = T30 - A(I+3,L) * B(J+0,L)
                        T31 = T31 - A(I+3,L) * B(J+1,L)
                        T32 = T32 - A(I+3,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                        T20 = T20 - A(I+2,L) * B(J+0,L)
                        T21 = T21 - A(I+2,L) * B(J+1,L)
                        T30 = T30 - A(I+3,L) * B(J+0,L)
                        T31 = T31 - A(I+3,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     T20 = BETA*C(I+2,J+0)
                     T30 = BETA*C(I+3,J+0)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T20 = T20 - A(I+2,L) * B(J+0,L)
                        T30 = T30 - A(I+3,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                     C(I+3,J+0) = T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T23 = BETA*C(I+2,J+3)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T02 = T02 - A(I+0,L) * B(J+2,L)
                        T03 = T03 - A(I+0,L) * B(J+3,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                        T12 = T12 - A(I+1,L) * B(J+2,L)
                        T13 = T13 - A(I+1,L) * B(J+3,L)
                        T20 = T20 - A(I+2,L) * B(J+0,L)
                        T21 = T21 - A(I+2,L) * B(J+1,L)
                        T22 = T22 - A(I+2,L) * B(J+2,L)
                        T23 = T23 - A(I+2,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T02 = T02 - A(I+0,L) * B(J+2,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                        T12 = T12 - A(I+1,L) * B(J+2,L)
                        T20 = T20 - A(I+2,L) * B(J+0,L)
                        T21 = T21 - A(I+2,L) * B(J+1,L)
                        T22 = T22 - A(I+2,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                        T20 = T20 - A(I+2,L) * B(J+0,L)
                        T21 = T21 - A(I+2,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     T20 = BETA*C(I+2,J+0)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T20 = T20 - A(I+2,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T02 = T02 - A(I+0,L) * B(J+2,L)
                        T03 = T03 - A(I+0,L) * B(J+3,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                        T12 = T12 - A(I+1,L) * B(J+2,L)
                        T13 = T13 - A(I+1,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T02 = T02 - A(I+0,L) * B(J+2,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                        T12 = T12 - A(I+1,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                        T11 = T11 - A(I+1,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T10 = T10 - A(I+1,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T02 = T02 - A(I+0,L) * B(J+2,L)
                        T03 = T03 - A(I+0,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                        T02 = T02 - A(I+0,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                        T01 = T01 - A(I+0,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     DO L = 0,K-1
                        T00 = T00 - A(I+0,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                  END SELECT
               END SELECT
            CASE DEFAULT
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0 
                     T13 = 0 
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     T33 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T03 = T03 + A(I+0,L) * B(J+3,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T13 = T13 + A(I+1,L) * B(J+3,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T22 = T22 + A(I+2,L) * B(J+2,L)
                        T23 = T23 + A(I+2,L) * B(J+3,L)
                        T30 = T30 + A(I+3,L) * B(J+0,L)
                        T31 = T31 + A(I+3,L) * B(J+1,L)
                        T32 = T32 + A(I+3,L) * B(J+2,L)
                        T33 = T33 + A(I+3,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+0,J+3) = BETA*C(I+0,J+3) + ALPHA*T03
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+1,J+3) = BETA*C(I+1,J+3) + ALPHA*T13
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+2,J+2) = BETA*C(I+2,J+2) + ALPHA*T22
                     C(I+2,J+3) = BETA*C(I+2,J+3) + ALPHA*T23
                     C(I+3,J+0) = BETA*C(I+3,J+0) + ALPHA*T30
                     C(I+3,J+1) = BETA*C(I+3,J+1) + ALPHA*T31
                     C(I+3,J+2) = BETA*C(I+3,J+2) + ALPHA*T32
                     C(I+3,J+3) = BETA*C(I+3,J+3) + ALPHA*T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T22 = T22 + A(I+2,L) * B(J+2,L)
                        T30 = T30 + A(I+3,L) * B(J+0,L)
                        T31 = T31 + A(I+3,L) * B(J+1,L)
                        T32 = T32 + A(I+3,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+2,J+2) = BETA*C(I+2,J+2) + ALPHA*T22
                     C(I+3,J+0) = BETA*C(I+3,J+0) + ALPHA*T30
                     C(I+3,J+1) = BETA*C(I+3,J+1) + ALPHA*T31
                     C(I+3,J+2) = BETA*C(I+3,J+2) + ALPHA*T32
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T30 = T30 + A(I+3,L) * B(J+0,L)
                        T31 = T31 + A(I+3,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+3,J+0) = BETA*C(I+3,J+0) + ALPHA*T30
                     C(I+3,J+1) = BETA*C(I+3,J+1) + ALPHA*T31
                  CASE (1)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T30 = T30 + A(I+3,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+3,J+0) = BETA*C(I+3,J+0) + ALPHA*T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T03 = T03 + A(I+0,L) * B(J+3,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T13 = T13 + A(I+1,L) * B(J+3,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T22 = T22 + A(I+2,L) * B(J+2,L)
                        T23 = T23 + A(I+2,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+0,J+3) = BETA*C(I+0,J+3) + ALPHA*T03
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+1,J+3) = BETA*C(I+1,J+3) + ALPHA*T13
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+2,J+2) = BETA*C(I+2,J+2) + ALPHA*T22
                     C(I+2,J+3) = BETA*C(I+2,J+3) + ALPHA*T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                        T22 = T22 + A(I+2,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+2,J+2) = BETA*C(I+2,J+2) + ALPHA*T22
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                        T21 = T21 + A(I+2,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T20 = T20 + A(I+2,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T03 = T03 + A(I+0,L) * B(J+3,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                        T13 = T13 + A(I+1,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+0,J+3) = BETA*C(I+0,J+3) + ALPHA*T03
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+1,J+3) = BETA*C(I+1,J+3) + ALPHA*T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                        T12 = T12 + A(I+1,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                        T11 = T11 + A(I+1,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T10 = T10 + A(I+1,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                        T03 = T03 + A(I+0,L) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+0,J+3) = BETA*C(I+0,J+3) + ALPHA*T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                        T02 = T02 + A(I+0,L) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                        T01 = T01 + A(I+0,L) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                  CASE (1)
                     T00 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(I+0,L) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                  END SELECT
               END SELECT
            END SELECT
         CASE (2)        ! 'T','N'
            SELECT CASE (AL)
            CASE (1)
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T23 = BETA*C(I+2,J+3)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     T32 = BETA*C(I+3,J+2)
                     T33 = BETA*C(I+3,J+3)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T03 = T03 + A(L,I+0) * B(L,J+3)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T13 = T13 + A(L,I+1) * B(L,J+3)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T22 = T22 + A(L,I+2) * B(L,J+2)
                        T23 = T23 + A(L,I+2) * B(L,J+3)
                        T30 = T30 + A(L,I+3) * B(L,J+0)
                        T31 = T31 + A(L,I+3) * B(L,J+1)
                        T32 = T32 + A(L,I+3) * B(L,J+2)
                        T33 = T33 + A(L,I+3) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                     C(I+3,J+3) = T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     T32 = BETA*C(I+3,J+2)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T22 = T22 + A(L,I+2) * B(L,J+2)
                        T30 = T30 + A(L,I+3) * B(L,J+0)
                        T31 = T31 + A(L,I+3) * B(L,J+1)
                        T32 = T32 + A(L,I+3) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T30 = T30 + A(L,I+3) * B(L,J+0)
                        T31 = T31 + A(L,I+3) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     T20 = BETA*C(I+2,J+0)
                     T30 = BETA*C(I+3,J+0)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T30 = T30 + A(L,I+3) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                     C(I+3,J+0) = T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T23 = BETA*C(I+2,J+3)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T03 = T03 + A(L,I+0) * B(L,J+3)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T13 = T13 + A(L,I+1) * B(L,J+3)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T22 = T22 + A(L,I+2) * B(L,J+2)
                        T23 = T23 + A(L,I+2) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T22 = T22 + A(L,I+2) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     T20 = BETA*C(I+2,J+0)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T03 = T03 + A(L,I+0) * B(L,J+3)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T13 = T13 + A(L,I+1) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T03 = T03 + A(L,I+0) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                  END SELECT
               END SELECT
            CASE (-1)
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T23 = BETA*C(I+2,J+3)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     T32 = BETA*C(I+3,J+2)
                     T33 = BETA*C(I+3,J+3)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T02 = T02 - A(L,I+0) * B(L,J+2)
                        T03 = T03 - A(L,I+0) * B(L,J+3)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                        T12 = T12 - A(L,I+1) * B(L,J+2)
                        T13 = T13 - A(L,I+1) * B(L,J+3)
                        T20 = T20 - A(L,I+2) * B(L,J+0)
                        T21 = T21 - A(L,I+2) * B(L,J+1)
                        T22 = T22 - A(L,I+2) * B(L,J+2)
                        T23 = T23 - A(L,I+2) * B(L,J+3)
                        T30 = T30 - A(L,I+3) * B(L,J+0)
                        T31 = T31 - A(L,I+3) * B(L,J+1)
                        T32 = T32 - A(L,I+3) * B(L,J+2)
                        T33 = T33 - A(L,I+3) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                     C(I+3,J+3) = T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     T32 = BETA*C(I+3,J+2)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T02 = T02 - A(L,I+0) * B(L,J+2)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                        T12 = T12 - A(L,I+1) * B(L,J+2)
                        T20 = T20 - A(L,I+2) * B(L,J+0)
                        T21 = T21 - A(L,I+2) * B(L,J+1)
                        T22 = T22 - A(L,I+2) * B(L,J+2)
                        T30 = T30 - A(L,I+3) * B(L,J+0)
                        T31 = T31 - A(L,I+3) * B(L,J+1)
                        T32 = T32 - A(L,I+3) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                        T20 = T20 - A(L,I+2) * B(L,J+0)
                        T21 = T21 - A(L,I+2) * B(L,J+1)
                        T30 = T30 - A(L,I+3) * B(L,J+0)
                        T31 = T31 - A(L,I+3) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     T20 = BETA*C(I+2,J+0)
                     T30 = BETA*C(I+3,J+0)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T20 = T20 - A(L,I+2) * B(L,J+0)
                        T30 = T30 - A(L,I+3) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                     C(I+3,J+0) = T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T23 = BETA*C(I+2,J+3)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T02 = T02 - A(L,I+0) * B(L,J+2)
                        T03 = T03 - A(L,I+0) * B(L,J+3)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                        T12 = T12 - A(L,I+1) * B(L,J+2)
                        T13 = T13 - A(L,I+1) * B(L,J+3)
                        T20 = T20 - A(L,I+2) * B(L,J+0)
                        T21 = T21 - A(L,I+2) * B(L,J+1)
                        T22 = T22 - A(L,I+2) * B(L,J+2)
                        T23 = T23 - A(L,I+2) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T02 = T02 - A(L,I+0) * B(L,J+2)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                        T12 = T12 - A(L,I+1) * B(L,J+2)
                        T20 = T20 - A(L,I+2) * B(L,J+0)
                        T21 = T21 - A(L,I+2) * B(L,J+1)
                        T22 = T22 - A(L,I+2) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                        T20 = T20 - A(L,I+2) * B(L,J+0)
                        T21 = T21 - A(L,I+2) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     T20 = BETA*C(I+2,J+0)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T20 = T20 - A(L,I+2) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T02 = T02 - A(L,I+0) * B(L,J+2)
                        T03 = T03 - A(L,I+0) * B(L,J+3)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                        T12 = T12 - A(L,I+1) * B(L,J+2)
                        T13 = T13 - A(L,I+1) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T02 = T02 - A(L,I+0) * B(L,J+2)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                        T12 = T12 - A(L,I+1) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                        T11 = T11 - A(L,I+1) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T10 = T10 - A(L,I+1) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T02 = T02 - A(L,I+0) * B(L,J+2)
                        T03 = T03 - A(L,I+0) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                        T02 = T02 - A(L,I+0) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                        T01 = T01 - A(L,I+0) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = T00
                  END SELECT
               END SELECT
            CASE DEFAULT
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0 
                     T13 = 0 
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     T33 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T03 = T03 + A(L,I+0) * B(L,J+3)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T13 = T13 + A(L,I+1) * B(L,J+3)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T22 = T22 + A(L,I+2) * B(L,J+2)
                        T23 = T23 + A(L,I+2) * B(L,J+3)
                        T30 = T30 + A(L,I+3) * B(L,J+0)
                        T31 = T31 + A(L,I+3) * B(L,J+1)
                        T32 = T32 + A(L,I+3) * B(L,J+2)
                        T33 = T33 + A(L,I+3) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+0,J+3) = BETA*C(I+0,J+3) + ALPHA*T03
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+1,J+3) = BETA*C(I+1,J+3) + ALPHA*T13
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+2,J+2) = BETA*C(I+2,J+2) + ALPHA*T22
                     C(I+2,J+3) = BETA*C(I+2,J+3) + ALPHA*T23
                     C(I+3,J+0) = BETA*C(I+3,J+0) + ALPHA*T30
                     C(I+3,J+1) = BETA*C(I+3,J+1) + ALPHA*T31
                     C(I+3,J+2) = BETA*C(I+3,J+2) + ALPHA*T32
                     C(I+3,J+3) = BETA*C(I+3,J+3) + ALPHA*T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T22 = T22 + A(L,I+2) * B(L,J+2)
                        T30 = T30 + A(L,I+3) * B(L,J+0)
                        T31 = T31 + A(L,I+3) * B(L,J+1)
                        T32 = T32 + A(L,I+3) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+2,J+2) = BETA*C(I+2,J+2) + ALPHA*T22
                     C(I+3,J+0) = BETA*C(I+3,J+0) + ALPHA*T30
                     C(I+3,J+1) = BETA*C(I+3,J+1) + ALPHA*T31
                     C(I+3,J+2) = BETA*C(I+3,J+2) + ALPHA*T32
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T30 = T30 + A(L,I+3) * B(L,J+0)
                        T31 = T31 + A(L,I+3) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+3,J+0) = BETA*C(I+3,J+0) + ALPHA*T30
                     C(I+3,J+1) = BETA*C(I+3,J+1) + ALPHA*T31
                  CASE (1)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T30 = T30 + A(L,I+3) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+3,J+0) = BETA*C(I+3,J+0) + ALPHA*T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T03 = T03 + A(L,I+0) * B(L,J+3)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T13 = T13 + A(L,I+1) * B(L,J+3)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T22 = T22 + A(L,I+2) * B(L,J+2)
                        T23 = T23 + A(L,I+2) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+0,J+3) = BETA*C(I+0,J+3) + ALPHA*T03
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+1,J+3) = BETA*C(I+1,J+3) + ALPHA*T13
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+2,J+2) = BETA*C(I+2,J+2) + ALPHA*T22
                     C(I+2,J+3) = BETA*C(I+2,J+3) + ALPHA*T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                        T22 = T22 + A(L,I+2) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+2,J+2) = BETA*C(I+2,J+2) + ALPHA*T22
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                        T21 = T21 + A(L,I+2) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T20 = T20 + A(L,I+2) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T03 = T03 + A(L,I+0) * B(L,J+3)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                        T13 = T13 + A(L,I+1) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+0,J+3) = BETA*C(I+0,J+3) + ALPHA*T03
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+1,J+3) = BETA*C(I+1,J+3) + ALPHA*T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                        T12 = T12 + A(L,I+1) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                        T11 = T11 + A(L,I+1) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T10 = T10 + A(L,I+1) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                        T03 = T03 + A(L,I+0) * B(L,J+3)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+0,J+3) = BETA*C(I+0,J+3) + ALPHA*T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                        T02 = T02 + A(L,I+0) * B(L,J+2)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                        T01 = T01 + A(L,I+0) * B(L,J+1)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                  CASE (1)
                     T00 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(L,J+0)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                  END SELECT
               END SELECT
            END SELECT
         CASE (3)        ! 'T','T'
            SELECT CASE (AL)
            CASE (1)
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T23 = BETA*C(I+2,J+3)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     T32 = BETA*C(I+3,J+2)
                     T33 = BETA*C(I+3,J+3)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T03 = T03 + A(L,I+0) * B(J+3,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T13 = T13 + A(L,I+1) * B(J+3,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T22 = T22 + A(L,I+2) * B(J+2,L)
                        T23 = T23 + A(L,I+2) * B(J+3,L)
                        T30 = T30 + A(L,I+3) * B(J+0,L)
                        T31 = T31 + A(L,I+3) * B(J+1,L)
                        T32 = T32 + A(L,I+3) * B(J+2,L)
                        T33 = T33 + A(L,I+3) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                     C(I+3,J+3) = T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     T32 = BETA*C(I+3,J+2)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T22 = T22 + A(L,I+2) * B(J+2,L)
                        T30 = T30 + A(L,I+3) * B(J+0,L)
                        T31 = T31 + A(L,I+3) * B(J+1,L)
                        T32 = T32 + A(L,I+3) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T30 = T30 + A(L,I+3) * B(J+0,L)
                        T31 = T31 + A(L,I+3) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     T20 = BETA*C(I+2,J+0)
                     T30 = BETA*C(I+3,J+0)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T30 = T30 + A(L,I+3) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                     C(I+3,J+0) = T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T23 = BETA*C(I+2,J+3)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T03 = T03 + A(L,I+0) * B(J+3,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T13 = T13 + A(L,I+1) * B(J+3,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T22 = T22 + A(L,I+2) * B(J+2,L)
                        T23 = T23 + A(L,I+2) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T22 = T22 + A(L,I+2) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     T20 = BETA*C(I+2,J+0)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T03 = T03 + A(L,I+0) * B(J+3,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T13 = T13 + A(L,I+1) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T03 = T03 + A(L,I+0) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                  END SELECT
               END SELECT
            CASE (-1)
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T23 = BETA*C(I+2,J+3)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     T32 = BETA*C(I+3,J+2)
                     T33 = BETA*C(I+3,J+3)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T02 = T02 - A(L,I+0) * B(J+2,L)
                        T03 = T03 - A(L,I+0) * B(J+3,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                        T12 = T12 - A(L,I+1) * B(J+2,L)
                        T13 = T13 - A(L,I+1) * B(J+3,L)
                        T20 = T20 - A(L,I+2) * B(J+0,L)
                        T21 = T21 - A(L,I+2) * B(J+1,L)
                        T22 = T22 - A(L,I+2) * B(J+2,L)
                        T23 = T23 - A(L,I+2) * B(J+3,L)
                        T30 = T30 - A(L,I+3) * B(J+0,L)
                        T31 = T31 - A(L,I+3) * B(J+1,L)
                        T32 = T32 - A(L,I+3) * B(J+2,L)
                        T33 = T33 - A(L,I+3) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                     C(I+3,J+3) = T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     T32 = BETA*C(I+3,J+2)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T02 = T02 - A(L,I+0) * B(J+2,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                        T12 = T12 - A(L,I+1) * B(J+2,L)
                        T20 = T20 - A(L,I+2) * B(J+0,L)
                        T21 = T21 - A(L,I+2) * B(J+1,L)
                        T22 = T22 - A(L,I+2) * B(J+2,L)
                        T30 = T30 - A(L,I+3) * B(J+0,L)
                        T31 = T31 - A(L,I+3) * B(J+1,L)
                        T32 = T32 - A(L,I+3) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                     C(I+3,J+2) = T32
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T30 = BETA*C(I+3,J+0)
                     T31 = BETA*C(I+3,J+1)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                        T20 = T20 - A(L,I+2) * B(J+0,L)
                        T21 = T21 - A(L,I+2) * B(J+1,L)
                        T30 = T30 - A(L,I+3) * B(J+0,L)
                        T31 = T31 - A(L,I+3) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+3,J+0) = T30
                     C(I+3,J+1) = T31
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     T20 = BETA*C(I+2,J+0)
                     T30 = BETA*C(I+3,J+0)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T20 = T20 - A(L,I+2) * B(J+0,L)
                        T30 = T30 - A(L,I+3) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                     C(I+3,J+0) = T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     T23 = BETA*C(I+2,J+3)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T02 = T02 - A(L,I+0) * B(J+2,L)
                        T03 = T03 - A(L,I+0) * B(J+3,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                        T12 = T12 - A(L,I+1) * B(J+2,L)
                        T13 = T13 - A(L,I+1) * B(J+3,L)
                        T20 = T20 - A(L,I+2) * B(J+0,L)
                        T21 = T21 - A(L,I+2) * B(J+1,L)
                        T22 = T22 - A(L,I+2) * B(J+2,L)
                        T23 = T23 - A(L,I+2) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                     C(I+2,J+3) = T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     T22 = BETA*C(I+2,J+2)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T02 = T02 - A(L,I+0) * B(J+2,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                        T12 = T12 - A(L,I+1) * B(J+2,L)
                        T20 = T20 - A(L,I+2) * B(J+0,L)
                        T21 = T21 - A(L,I+2) * B(J+1,L)
                        T22 = T22 - A(L,I+2) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                     C(I+2,J+2) = T22
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T20 = BETA*C(I+2,J+0)
                     T21 = BETA*C(I+2,J+1)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                        T20 = T20 - A(L,I+2) * B(J+0,L)
                        T21 = T21 - A(L,I+2) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+2,J+0) = T20
                     C(I+2,J+1) = T21
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     T20 = BETA*C(I+2,J+0)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T20 = T20 - A(L,I+2) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                     C(I+2,J+0) = T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     T13 = BETA*C(I+1,J+3)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T02 = T02 - A(L,I+0) * B(J+2,L)
                        T03 = T03 - A(L,I+0) * B(J+3,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                        T12 = T12 - A(L,I+1) * B(J+2,L)
                        T13 = T13 - A(L,I+1) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                     C(I+1,J+3) = T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     T12 = BETA*C(I+1,J+2)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T02 = T02 - A(L,I+0) * B(J+2,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                        T12 = T12 - A(L,I+1) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                     C(I+1,J+2) = T12
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T10 = BETA*C(I+1,J+0)
                     T11 = BETA*C(I+1,J+1)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                        T11 = T11 - A(L,I+1) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+1,J+0) = T10
                     C(I+1,J+1) = T11
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     T10 = BETA*C(I+1,J+0)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T10 = T10 - A(L,I+1) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+1,J+0) = T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     T03 = BETA*C(I+0,J+3)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T02 = T02 - A(L,I+0) * B(J+2,L)
                        T03 = T03 - A(L,I+0) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                     C(I+0,J+3) = T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     T02 = BETA*C(I+0,J+2)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                        T02 = T02 - A(L,I+0) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                     C(I+0,J+2) = T02
                  CASE (2)
                     T00 = BETA*C(I+0,J+0)
                     T01 = BETA*C(I+0,J+1)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                        T01 = T01 - A(L,I+0) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = T00
                     C(I+0,J+1) = T01
                  CASE (1)
                     T00 = BETA*C(I+0,J+0)
                     DO L = 0,K-1
                        T00 = T00 - A(L,I+0) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = T00
                  END SELECT
               END SELECT
            CASE DEFAULT
               DO I = 0,M-4,4
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0 
                     T13 = 0 
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     T33 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T03 = T03 + A(L,I+0) * B(J+3,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T13 = T13 + A(L,I+1) * B(J+3,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T22 = T22 + A(L,I+2) * B(J+2,L)
                        T23 = T23 + A(L,I+2) * B(J+3,L)
                        T30 = T30 + A(L,I+3) * B(J+0,L)
                        T31 = T31 + A(L,I+3) * B(J+1,L)
                        T32 = T32 + A(L,I+3) * B(J+2,L)
                        T33 = T33 + A(L,I+3) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+0,J+3) = BETA*C(I+0,J+3) + ALPHA*T03
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+1,J+3) = BETA*C(I+1,J+3) + ALPHA*T13
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+2,J+2) = BETA*C(I+2,J+2) + ALPHA*T22
                     C(I+2,J+3) = BETA*C(I+2,J+3) + ALPHA*T23
                     C(I+3,J+0) = BETA*C(I+3,J+0) + ALPHA*T30
                     C(I+3,J+1) = BETA*C(I+3,J+1) + ALPHA*T31
                     C(I+3,J+2) = BETA*C(I+3,J+2) + ALPHA*T32
                     C(I+3,J+3) = BETA*C(I+3,J+3) + ALPHA*T33
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T30 = 0
                     T31 = 0
                     T32 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T22 = T22 + A(L,I+2) * B(J+2,L)
                        T30 = T30 + A(L,I+3) * B(J+0,L)
                        T31 = T31 + A(L,I+3) * B(J+1,L)
                        T32 = T32 + A(L,I+3) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+2,J+2) = BETA*C(I+2,J+2) + ALPHA*T22
                     C(I+3,J+0) = BETA*C(I+3,J+0) + ALPHA*T30
                     C(I+3,J+1) = BETA*C(I+3,J+1) + ALPHA*T31
                     C(I+3,J+2) = BETA*C(I+3,J+2) + ALPHA*T32
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T30 = T30 + A(L,I+3) * B(J+0,L)
                        T31 = T31 + A(L,I+3) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+3,J+0) = BETA*C(I+3,J+0) + ALPHA*T30
                     C(I+3,J+1) = BETA*C(I+3,J+1) + ALPHA*T31
                  CASE (1)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     T30 = 0
                     T31 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T30 = T30 + A(L,I+3) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+3,J+0) = BETA*C(I+3,J+0) + ALPHA*T30
                  END SELECT
               END DO
               SELECT CASE (M-I)
               CASE (3)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     T23 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T03 = T03 + A(L,I+0) * B(J+3,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T13 = T13 + A(L,I+1) * B(J+3,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T22 = T22 + A(L,I+2) * B(J+2,L)
                        T23 = T23 + A(L,I+2) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+0,J+3) = BETA*C(I+0,J+3) + ALPHA*T03
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+1,J+3) = BETA*C(I+1,J+3) + ALPHA*T13
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+2,J+2) = BETA*C(I+2,J+2) + ALPHA*T22
                     C(I+2,J+3) = BETA*C(I+2,J+3) + ALPHA*T23
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T20 = 0
                     T21 = 0
                     T22 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                        T22 = T22 + A(L,I+2) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                     C(I+2,J+2) = BETA*C(I+2,J+2) + ALPHA*T22
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     T20 = 0
                     T21 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                        T21 = T21 + A(L,I+2) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                     C(I+2,J+1) = BETA*C(I+2,J+1) + ALPHA*T21
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     T20 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T20 = T20 + A(L,I+2) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+2,J+0) = BETA*C(I+2,J+0) + ALPHA*T20
                  END SELECT
               CASE (2)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     T13 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T03 = T03 + A(L,I+0) * B(J+3,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                        T13 = T13 + A(L,I+1) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+0,J+3) = BETA*C(I+0,J+3) + ALPHA*T03
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                     C(I+1,J+3) = BETA*C(I+1,J+3) + ALPHA*T13
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T10 = 0
                     T11 = 0
                     T12 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                        T12 = T12 + A(L,I+1) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                     C(I+1,J+2) = BETA*C(I+1,J+2) + ALPHA*T12
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     T10 = 0
                     T11 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                        T11 = T11 + A(L,I+1) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                     C(I+1,J+1) = BETA*C(I+1,J+1) + ALPHA*T11
                  CASE (1)
                     T00 = 0
                     T10 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T10 = T10 + A(L,I+1) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+1,J+0) = BETA*C(I+1,J+0) + ALPHA*T10
                  END SELECT
               CASE (1)
                  DO J = 0,N-4,4
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     T03 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                        T03 = T03 + A(L,I+0) * B(J+3,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                     C(I+0,J+3) = BETA*C(I+0,J+3) + ALPHA*T03
                  END DO
                  SELECT CASE (N-J)
                  CASE (3)
                     T00 = 0
                     T01 = 0
                     T02 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                        T02 = T02 + A(L,I+0) * B(J+2,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                     C(I+0,J+2) = BETA*C(I+0,J+2) + ALPHA*T02
                  CASE (2)
                     T00 = 0
                     T01 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                        T01 = T01 + A(L,I+0) * B(J+1,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                     C(I+0,J+1) = BETA*C(I+0,J+1) + ALPHA*T01
                  CASE (1)
                     T00 = 0
                     DO L = 0,K-1
                        T00 = T00 + A(L,I+0) * B(J+0,L)
                     END DO
                     C(I+0,J+0) = BETA*C(I+0,J+0) + ALPHA*T00
                  END SELECT
               END SELECT
            END SELECT
         END SELECT
      END IF
      RETURN
      END SUBROUTINE RECSY_GEMM
