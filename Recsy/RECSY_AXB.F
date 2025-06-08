!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      SUBROUTINE RECSY_AXB(FA, FB, TA, TB, M1, M2, N1, N2, ALPHA, A, LDA, X, LDX, B, LDB, BETA, C, LDC, MACHINE, WORKSPACE)
!
!     Memory requirements
!
!     'U'/'V' and 'U'/'V'
!     MULT1:
!     M1*N1 + MAX(RECSY_MULT_LEFT(M1,N1,M2), RECSY_MULT_RIGHT(M1,N2,N1) =
!      M1*N1 + MAX(MIN(M1,M2)*N1, MIN(N1,N2)*M1)
!
!     MULT2:
!     M2*N2 + MAX(RECSY_MULT_RIGHT(M2,N2,N1), RECSY_MULT_LEFT(M1,N2,M2) =
!      M2*N2 + MAX(MIN(N1,N2)*M2, MIN(M1,M2)*N2)
!
!     M1        M2        N1        N2        WORKSPACE
!     M/2       M/2       N/2       N/2       M*N/2
!     M/2+1/2   M/2-1/2   N/2       N/2       M*N/2 + N/2
!     M/2-1/2   M/2+1/2   N/2       N/2       M*N/2 + N/2
!     M/2       M/2       N/2+1/2   N/2-1/2   M*N/2 + M/2
!     M/2+1/2   M/2-1/2   N/2+1/2   N/2-1/2   M*N/2 + max(M,N)/2
!     M/2       M/2       N         N         M*N
!     M         M         N/2       N/2       M*N
!     M         M         N         N         M*N*2
!     M         M         N/2+1/2   N/2-1/2   M*N/2 + M/2
!     M/2+1     M/2-1     N/2+1     N/2-1     M*N/2 + max(M,N)
!
!     M/2+1
!
!     'F' and 'U'/'V'
!      CALL RECSY_AXB('F','V','N','N', M1, M2,  N1, N2, -1D0, A(0,M1), LDA, C(M1,0), LDC, B(0,N1), LDB, SCALE, C(0,N1), LDC, MACHINE, WORKSPACE)
!
!     MULT1:
!     M1*N1 + MAX(RECSY_MULT_LEFT(M1,N1,M2), RECSY_MULT_RIGHT(M1,N2,N1) =
!      M1*N1 + MAX(MIN(M1,M2)*N1, MIN(N1,N2)*M1)
!
!     MULT2:
!     M2*N2 + MAX(RECSY_MULT_RIGHT(M2,N2,N1), RECSY_MULT_LEFT(M1,N2,M2) =
!      M2*N2 + MAX(MIN(N1,N2)*M2, 0)
!
!     M1        M2        N1        N2        WORKSPACE
!     M/2       M/2       N/2       N/2       M*N/2
!     M/2+1/2   M/2-1/2   N/2       N/2       M*N/2 + N/2
!     M/2       M/2       N/2+1/2   N/2-1/2   M*N/2 + M/2
!     M/2+1/2   M/2-1/2   N/2+1/2   N/2-1/2   M*N/2 + max(M,N)/2
!     M/2       M/2       N         N         M*N
!     M         M         N/2       N/2       M*N
!     M         M         N         N         M*N*2
!
!
!
!

      IMPLICIT NONE
      INTEGER M1,M2,N1,N2,LDA,LDX,LDB,LDC,I,J
      DOUBLE PRECISION A(0:lda-1, 0:*)
      DOUBLE PRECISION B(0:ldb-1, 0:*)
      DOUBLE PRECISION C(0:ldc-1, 0:*)
      DOUBLE PRECISION X(0:ldx-1, 0:*), ALPHA, BETA
      DOUBLE PRECISION WORKSPACE(0:*), MACHINE(0:*)
      CHARACTER*1 TA,TB,FA,FB
      DOUBLE PRECISION MULT1A,MULT1B,MULT2A,MULT2B
      INTEGER BIGGEMM

      BIGGEMM = MACHINE(4)

      IF (FA.eq.'F') THEN
         MULT1A = (2D0*m1)*m2*n1
         MULT2B = (2D0*m1)*m2*n2
      ELSE
         MULT1A = (1D0*m1)*m1*n1+(2D0*m1)*(m2-m1)*n1
         MULT2B = (1D0*m1)*m2*n2+(2D0*m1)*(m2-m1)*n2
      ENDIF

      IF (FB.eq.'F') THEN
         MULT1B = (2D0*m1)*n1*n2
         MULT2A = (2D0*m2)*n1*n2
      ELSE
         MULT1B = (1D0*m1)*n2*n2+(2D0*m1)*n2*(n1-n2) ! 5*10*10+2*5*10*0
         MULT2A = (1D0*m2)*n2*n2+(2D0*m2)*n2*(n1-n2) ! 5*10*10+2*5*10*0
      ENDIF

      IF ((MULT1A+MULT1B) .lt. (MULT2A+MULT2B)) THEN

         CALL RECSY_MULT_LEFT(FA, TA, M1, N1, M2, ALPHA, A, LDA, X, LDX, 0D0, WORKSPACE, M1, MACHINE, WORKSPACE(M1*N1))
         CALL RECSY_MULT_RIGHT(FB, TB, M1, N2, N1, 1D0, WORKSPACE, M1, B, LDB, BETA, C, LDC, MACHINE, WORKSPACE(M1*N1))

      ELSE

         CALL RECSY_MULT_RIGHT(FB, TB, M2, N2, N1, ALPHA, X, LDX, B, LDB, 0D0, WORKSPACE, M2, MACHINE, WORKSPACE(M2*N2))
         CALL RECSY_MULT_LEFT(FA, TA, M1, N2, M2, 1D0, A, LDA, WORKSPACE, M2, BETA, C, LDC, MACHINE, WORKSPACE(M2*N2))

      END IF
      END SUBROUTINE RECSY_AXB

      SUBROUTINE RECSY_TRMM_Q(SIDE, TRANS, M, N, ALPHA, A, LDA, B, LDB)

      IMPLICIT NONE
      INTEGER LDA,LDB,MODE,I,M,N
      DOUBLE PRECISION :: A(0:LDA-1, 0:*)
      DOUBLE PRECISION :: B(0:LDB-1, 0:*), ALPHA
      CHARACTER*1      :: SIDE
      CHARACTER*1      :: TRANS

      IF (SIDE.EQ.'L') THEN
         IF (TRANS.EQ.'N') THEN
            DO I = 0,M-2
               IF (A(I+1,I).NE.0D0) THEN
                  CALL DAXPY(N, ALPHA*A(I+1,I), B(I,0), LDB, B(I+1,0), LDB)
               END IF
            END DO
         ELSE
            DO I = 0,M-2
               IF (A(I+1,I).NE.0D0) THEN
                  CALL DAXPY(N, ALPHA*A(I+1,I), B(I+1,0), LDB, B(I,0), LDB)
               END IF
            END DO
         END IF
      ELSE
         IF (TRANS.EQ.'N') THEN
            DO I = 0,N-2
               IF (A(I+1,I).NE.0D0) THEN
                  CALL DAXPY(M, ALPHA*A(I+1,I), B(I,0), 1, B(I+1,0), 1)
               END IF
            END DO
         ELSE
            DO I = 0,N-2
               IF (A(I+1,I).NE.0D0) THEN
                  CALL DAXPY(M, ALPHA*A(I+1,I), B(I+1,0), 1, B(I,0), 1)
               END IF
            END DO
         END IF
      END IF

      END SUBROUTINE RECSY_TRMM_Q

      SUBROUTINE RECSY_TRMM_Q2(SIDE, TRANS, M, N, ALPHA, A, LDA, B, LDB, C, LDC)

      IMPLICIT NONE
      INTEGER LDA,LDB,LDC,MODE,I,M,N
      DOUBLE PRECISION :: A(0:LDA-1, 0:*)
      DOUBLE PRECISION :: B(0:LDB-1, 0:*), ALPHA
      DOUBLE PRECISION :: C(0:LDC-1, 0:*)
      CHARACTER*1      :: SIDE
      CHARACTER*1      :: TRANS

      IF (SIDE.EQ.'L') THEN
         IF (TRANS.EQ.'N') THEN
            DO I = 0,M-2
               IF (A(I+1,I).NE.0D0) THEN
                  CALL DAXPY(N, ALPHA*A(I+1,I), B(I,0), LDB, C(I+1,0), LDC)
               END IF
            END DO
         ELSE
            DO I = 0,M-2
               IF (A(I+1,I).NE.0D0) THEN
                  CALL DAXPY(N, ALPHA*A(I+1,I), B(I+1,0), LDB, C(I,0), LDC)
               END IF
            END DO
         END IF
      ELSE
         IF (TRANS.EQ.'N') THEN
            DO I = 0,N-2
               IF (A(I+1,I).NE.0D0) THEN
                  CALL DAXPY(M, ALPHA*A(I+1,I), B(0,I+1), 1, C(0,I), 1)
               END IF
            END DO
         ELSE
            DO I = 0,N-2
               IF (A(I+1,I).NE.0D0) THEN
                  CALL DAXPY(M, ALPHA*A(I+1,I), B(0,I), 1, C(0,I+1), 1)
               END IF
            END DO
         END IF
      END IF

      END SUBROUTINE RECSY_TRMM_Q2


      SUBROUTINE RECSY_MULT_LEFT(FA, TA, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, MACHINE, WORKSPACE)

!     C = op(A) * B
!     A = Upper(U), Quasi(V), Full(F)
!
!     Memory requirements:
!     BETA.EQ.0.OR.FA.EQ.'F'        0
!     OTHERWISE                    MIN(K,M)*N
!

      IMPLICIT NONE
      INTEGER M,N,K,LDA,LDB,LDC,I,J
      DOUBLE PRECISION A(0:LDA-1, 0:*)
      DOUBLE PRECISION B(0:LDB-1, 0:*)
      DOUBLE PRECISION C(0:LDC-1, 0:*)
      DOUBLE PRECISION MACHINE(0:*), WORKSPACE(0:*), ALPHA, BETA
      CHARACTER*1 TA,FA
      DOUBLE PRECISION MULT1A,MULT1B,MULT2A,MULT2B

      IF (FA.eq.'F'.or.MAX(MAX(M,N),K).lt.MACHINE(4)) THEN
         IF (TA.EQ.'N') THEN
            CALL RECSY_GEMM(0, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, MACHINE)
         ELSE
            CALL RECSY_GEMM(1, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, MACHINE)
         END IF
         RETURN
      END IF
      IF (TA.EQ.'N') THEN
         IF (K.GE.M) THEN
!     K     N
!   AAAaa BBBB CCCC
! M  AAaa BBBB CCCC
!     Aaa BBBB CCCC
!         bbbb
!         bbbb
            IF (BETA.EQ.0D0) THEN
               DO J=0,N-1
                  DO I=0,M-1
                     C(I,J) = B(I,J)
                  END DO
               END DO
               CALL DTRMM('L','U','N','N', M, N, ALPHA, A, LDA, C, LDC)
               IF (FA.eq.'V') THEN
                  CALL RECSY_TRMM_Q2('L','N', M, N, ALPHA, A, LDA, B, LDB, C, LDC)
               END IF
            ELSE
               DO J=0,N-1
                  DO I=0,M-1
                     WORKSPACE(I+J*M) = B(I,J)
                  END DO
               END DO
               CALL DTRMM('L','U','N','N', M, N, ALPHA, A, LDA, WORKSPACE, M)
               IF (FA.eq.'V') THEN
                  CALL RECSY_TRMM_Q2('L','N', M, N, ALPHA, A, LDA, B, LDB, WORKSPACE, M)
               END IF
               DO J=0,N-1
                  DO I=0,M-1
                     C(I,J) = BETA * C(I,J) + WORKSPACE(I+J*M)
                  END DO
               END DO
            END IF
            IF (K.GT.M) THEN
               CALL RECSY_GEMM(0, M, N, K-M, ALPHA, A(0,M), LDA, B(M,0), LDB, 1D0, C, LDC, MACHINE)
            END IF
         ELSE
!     K     N
!   aaa  BBBB cccc
! M aaa  BBBB cccc
!   AAA  BBBB CCCC
!    AA       CCCC
!     A       CCCC
            IF (BETA.EQ.0D0) THEN
               DO J=0,N-1
                  DO I=0,K-1
                     C(I+M-K,J) = B(I,J)
                  END DO
               END DO
               CALL DTRMM('L','U','N','N', K, N, ALPHA, A(M-K,0), LDA, C(M-K,0), LDC)
               IF (FA.eq.'V') THEN
                  CALL RECSY_TRMM_Q2('L','N', K, N, ALPHA, A(M-K,0), LDA, B, LDB, C(M-K,0), LDC)
               END IF
            ELSE
               DO J=0,N-1
                  DO I=0,K-1
                     WORKSPACE(I+J*K) = B(I,J)
                  END DO
               END DO
               CALL DTRMM('L','U','N','N', K, N, ALPHA, A(M-K,0), LDA, WORKSPACE, K)
               IF (FA.eq.'V') THEN
                  CALL RECSY_TRMM_Q2('L','N', K, N, ALPHA, A(M-K,0), LDA, B, LDB, WORKSPACE, K)
               END IF
               DO J=0,N-1
                  DO I=0,K-1
                     C(I+M-K,J) = BETA * C(I+M-K,J) + WORKSPACE(I+J*K)
                  END DO
               END DO
            END IF
            IF (M.GT.K) THEN
               CALL RECSY_GEMM(0, M-K, N, K, ALPHA, A(0,0), LDA, B(0,0), LDB, BETA, C(0,0), LDC, MACHINE)
            END IF
         END IF
      ELSE
         IF (M.GE.K) THEN
!   K      
!  A   BBBB  CCCC
!M AA  BBBB  CCCC
!  AAA BBBB  CCCC
!  aaa       cccc
!  aaa       cccc
            IF (BETA.EQ.0D0) THEN 
               DO J=0,N-1
                  DO I=0,K-1
                     C(I,J) = B(I,J)
                  END DO
               END DO
               CALL DTRMM('L','U','T', 'N', K, N, ALPHA, A, LDA, C, LDC)
               IF (FA.eq.'V') THEN
                  CALL RECSY_TRMM_Q2('L','T', K, N, ALPHA, A, LDA, B, LDB, C, LDC)
               END IF
            ELSE
               DO J=0,N-1
                  DO I=0,K-1
                     WORKSPACE(I+J*K) = B(I,J)
                  END DO
               END DO
               CALL DTRMM('L','U','T', 'N', K, N, ALPHA, A, LDA, WORKSPACE, K)
               IF (FA.eq.'V') THEN
                  CALL RECSY_TRMM_Q2('L','T', K, N, ALPHA, A, LDA, B, LDB, WORKSPACE, K)
               END IF
               DO J=0,N-1
                  DO I=0,K-1
                     C(I,J) = BETA * C(I,J) + WORKSPACE(I+J*K)
                  END DO
               END DO
            END IF
            IF (M.GT.K) THEN
               CALL RECSY_GEMM(1, M-K, N, K, ALPHA, A(0,K), LDA, B(0,0), LDB, BETA, C(K,0), LDC, MACHINE)
            END IF
         ELSE
!   K      N
!  aaA   bbbb  CCCC
!M aaAA  bbbb  CCCC
!  aaAAA BBBB  CCCC
!        BBBB  
!        BBBB  
            IF (BETA.EQ.0D0) THEN 
               DO J=0,N-1
                  DO I=0,M-1
                     C(I,J) = B(I+K-M,J)
                  END DO
               END DO
               CALL DTRMM('L','U','T', 'N', M, N, ALPHA, A(K-M,0), LDA, C, LDC)
               IF (FA.eq.'V') THEN
                  CALL RECSY_TRMM_Q2('L','T', M, N, ALPHA, A(K-M,0), LDA, B(K-M,0), LDB, C, LDC)
               END IF
            ELSE
               DO J=0,N-1
                  DO I=0,M-1
                     WORKSPACE(I+J*M) = B(I+K-M,J)
                  END DO
               END DO
               CALL DTRMM('L','U','T', 'N', M, N, ALPHA, A(K-M,0), LDA, WORKSPACE, M)
               IF (FA.eq.'V') THEN
                  CALL RECSY_TRMM_Q2('L','T', M, N, ALPHA, A(K-M,0), LDA, B(K-M,0), LDB, WORKSPACE, M)
               END IF
               DO J=0,N-1
                  DO I=0,M-1
                     C(I,J) = BETA * C(I,J) + WORKSPACE(I+J*M)
                  END DO
               END DO
            END IF
            IF (K.GT.M) THEN
               CALL RECSY_GEMM(1, M, N, K-M, ALPHA, A(0,0), LDA, B(0,0), LDB, 1D0, C(0,0), LDC, MACHINE)
            END IF
         END IF
      END IF
      END SUBROUTINE RECSY_MULT_LEFT
            
      SUBROUTINE RECSY_MULT_RIGHT(FB, TB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, MACHINE, WORKSPACE)

!     C = A * op(B)
!     B = Upper(U), Quasi(V), Full(F)
!
!     Memory requirements:
!     BETA.EQ.0.OR.FB.EQ.'F'        0
!     OTHERWISE                     MIN(K,N)*M
!
      IMPLICIT NONE
      INTEGER M,N,K,LDA,LDB,LDC,I,J
      DOUBLE PRECISION A(0:LDA-1, 0:*)
      DOUBLE PRECISION B(0:LDB-1, 0:*)
      DOUBLE PRECISION C(0:LDC-1, 0:N)
      DOUBLE PRECISION WORKSPACE(0:*), MACHINE(0:*), ALPHA, BETA
      CHARACTER*1 TB,FB
      DOUBLE PRECISION MULT1A,MULT1B,MULT2A,MULT2B

      IF (FB.eq.'F'.or.MAX(MAX(M,N),K).lt.MACHINE(4)) THEN
         IF (TB.EQ.'N') THEN
            CALL RECSY_GEMM(0, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, MACHINE)
         ELSE
            CALL RECSY_GEMM(2, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, MACHINE)
         END IF
         RETURN
      END IF
      IF (TB.EQ.'N') THEN
         IF (N.GE.K) THEN
!   K    N
!   AA BBbbb CCccc
! M AA  Bbbb CCccc
!   AA       CCccc
!   AA       CCccc
            IF (BETA.EQ.0D0) THEN
               DO J=0,K-1
                  DO I=0,M-1
                     C(I,J) = A(I,J)
                  END DO
               END DO
               CALL DTRMM('R','U','N','N', M, K, ALPHA, B, LDB, C, LDC)
               IF (FB.eq.'V') THEN
                  CALL RECSY_TRMM_Q2('R','N', M, K, ALPHA, B, LDB, A, LDA, C, LDC)
               END IF
            ELSE
               DO J=0,K-1
                  DO I=0,M-1
                     WORKSPACE(I+J*M) = A(I,J)
                  END DO
               END DO
               CALL DTRMM('R','U','N','N', M, K, ALPHA, B, LDB, WORKSPACE, M)
               IF (FB.eq.'V') THEN
                  CALL RECSY_TRMM_Q2('R','N', M, K, ALPHA, B, LDB, A, LDA, WORKSPACE, M)
               END IF
               DO J=0,K-1
                  DO I=0,M-1
                     C(I,J) = BETA * C(I,J) + WORKSPACE(I+J*M)
                  END DO
               END DO
            END IF
            IF (N.GT.K) THEN
               CALL RECSY_GEMM(0, M, N-K, K, ALPHA, A(0,0), LDA, B(0,K), LDB, BETA, C(0,K), LDC, MACHINE)
            END IF
         ELSE
!   K    N
!   aAAA bbb CCC
! M aAAA BBB CCC
!   aAAA  BB CCC
!   aAAA   B CCC
!   aAAA     CCC
            IF (BETA.EQ.0D0) THEN
               DO J=0,N-1
                  DO I=0,M-1
                     C(I,J) = A(I,J+K-N)
                  END DO
               END DO
               CALL DTRMM('R','U','N','N', M, N, ALPHA, B(K-N,0), LDB, C, LDC)
               IF (FB.eq.'V') THEN
                  CALL RECSY_TRMM_Q2('R','N', M, N, ALPHA, B(K-N,0), LDB, A(0,K-N), LDA, C, LDC)
               END IF
            ELSE
               DO J=0,N-1
                  DO I=0,M-1
                     WORKSPACE(I+J*M) = A(I,J+K-N)
                  END DO
               END DO
               CALL DTRMM('R','U','N','N', M, N, ALPHA, B(K-N,0), LDB, WORKSPACE, M)
               IF (FB.eq.'V') THEN
                  CALL RECSY_TRMM_Q2('R','N', M, N, ALPHA, B(K-N,0), LDB, A(0,K-N), LDA, WORKSPACE, M)
               END IF
               DO J=0,N-1
                  DO I=0,M-1
                     C(I,J) = BETA * C(I,J) + WORKSPACE(I+J*M)
                  END DO
               END DO
            END IF
            IF (K.GT.N) THEN
               CALL RECSY_GEMM(0, M, N, K-N, ALPHA, A(0,0), LDA, B(0,0), LDB, 1D0, C(0,0), LDC, MACHINE)
            END IF
         END IF
      ELSE
         IF (K.GE.N) THEN
!    K     N
! M AAaaa  B   CC
!   AAaaa  BB  CC
!   AAaaa  bb  CC
!          bb
!          bb
            IF (BETA.EQ.0D0) THEN
               DO J=0,N-1
                  DO I=0,M-1
                     C(I,J) = A(I,J)
                  END DO
               END DO
               CALL DTRMM('R','U','T', 'N', M, N, ALPHA, B, LDB, C, LDC)
               IF (FB.eq.'V') THEN
                  CALL RECSY_TRMM_Q2('R','T', M, N, ALPHA, B, LDB, A, LDA, C, LDC)
               END IF
            ELSE
               DO J=0,N-1
                  DO I=0,M-1
                     WORKSPACE(I+J*M) = A(I,J)
                  END DO
               END DO
               CALL DTRMM('R','U','T', 'N', M, N, ALPHA, B, LDB, WORKSPACE, M)
               IF (FB.eq.'V') THEN
                  CALL RECSY_TRMM_Q2('R','T', M, N, ALPHA, B, LDB, A, LDA, WORKSPACE, M)
               END IF
               DO J=0,N-1
                  DO I=0,M-1
                     C(I,J) = BETA * C(I,J) + WORKSPACE(I+J*M)
                  END DO
               END DO
            END IF
            IF (K.GT.N) THEN
               CALL RECSY_GEMM(2, M, N, K-N, ALPHA, A(0,N), LDA, B(0,N), LDB, 1D0, C(0,0), LDC, MACHINE)
            END IF
         ELSE
!    K     N
! M AAA  bbB    ccCCC
!   AAA  bbBB   ccCCC
!   AAA  bbBBB  ccCCC
!   AAA         ccCCC
            IF (BETA.EQ.0D0) THEN
               DO J=0,K-1
                  DO I=0,M-1
                     C(I,J+N-K) = A(I,J)
                  END DO
               END DO
               CALL DTRMM('R','U','T', 'N', M, K, ALPHA, B(N-K,0), LDB, C(0,N-K), LDC)
               IF (FB.eq.'V') THEN
                  CALL RECSY_TRMM_Q2('R','T', M, K, ALPHA, B(N-K,0), LDB, A, LDA, C(0,N-K), LDC)
               END IF
            ELSE
               DO J=0,K-1
                  DO I=0,M-1
                     WORKSPACE(I+J*M) = A(I,J)
                  END DO
               END DO
               CALL DTRMM('R','U','T', 'N', M, K, ALPHA, B(N-K,0), LDB, WORKSPACE, M)
               IF (FB.eq.'V') THEN
                  CALL RECSY_TRMM_Q2('R','T', M, K, ALPHA, B(N-K,0), LDB, A, LDA, WORKSPACE, M)
               END IF
               DO J=0,K-1
                  DO I=0,M-1
                     C(I,J+N-K) = BETA * C(I,J+N-K) + WORKSPACE(I+J*M)
                  END DO
               END DO
            END IF
            IF (N.GT.K) THEN
               CALL RECSY_GEMM(2, M, N-K, K, ALPHA, A(0,0), LDA, B(0,0), LDB, BETA, C(0,0), LDC, MACHINE)
            END IF
         END IF
      END IF
      END SUBROUTINE RECSY_MULT_RIGHT


