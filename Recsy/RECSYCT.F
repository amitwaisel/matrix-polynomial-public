!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      RECURSIVE SUBROUTINE RECSYCT(UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE )
!     Solves op(A)*X + sign*X*op(B) = scale*C
!     A (M*M) and B (N*N) quasi-triangular
!     
!     Iff close eigenvalues are detected, INFO is non-zero and elements 
!     of the C (M*N) matrix are no longer valid
!
      IMPLICIT NONE
!
      INTENT(IN) M,N,A,LDA,B,LDB,LDC,UPLOSIGN,MACHINE
      INTENT(OUT) INFO,SCALE
      INTENT(INOUT) C
!     .. Scalar Arguments ..
      INTEGER          :: M, N, LDA, LDB, LDC, INFO, UPLOSIGN
      DOUBLE PRECISION :: SCALE
!     .. Array Arguments ..
      DOUBLE PRECISION :: A(0:LDA-1,0:*), B(0:LDB-1,0:*), C(0:LDC-1,0:*)
      DOUBLE PRECISION, OPTIONAL :: MACHINE(0:*)
!     .. Local Scalars ..
      INTEGER          :: I, J, K, M1, M2, N1, N2, INF
      DOUBLE PRECISION :: SIGN, SC
!     .. Local Arrays ..
      DOUBLE PRECISION :: Z(0:31,0:31)
      DOUBLE PRECISION :: VEC(0:15), MACHINETEMP(0:9)
!     .. Intrinsic Functions ..
      INTRINSIC           MOD
!     .. Explicit Declarations
!     .. Executable Statements ..
!
      IF (.NOT.PRESENT(MACHINE)) THEN
         CALL RECSY_MACHINE(MACHINETEMP)
         CALL RECSYCT(UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINETEMP )
         RETURN
      END IF

      SCALE = 1D0
      SIGN = 1D0-(MOD(UPLOSIGN,2))*2D0
      IF (M.LE.4.AND.N.LE.4) THEN
         DO J=0,N-1
            DO I=0,M-1
               VEC(I+J*M) = C(I,J)
            END DO
         END DO
         SELECT CASE (UPLOSIGN)
         CASE (0)
            CALL RECSYCT_TR0(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
         CASE (1)
            CALL RECSYCT_TR1(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
         CASE (2)
            CALL RECSYCT_TR2(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
         CASE (3)
            CALL RECSYCT_TR3(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
         CASE (4)
            CALL RECSYCT_TR4(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
         CASE (5)
            CALL RECSYCT_TR5(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
         CASE (6)
            CALL RECSYCT_TR6(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
         CASE (7)
            CALL RECSYCT_TR7(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
         END SELECT
         IF (INFO.NE.0) THEN
            DO J=0,M*N-1
               DO I=0,M*N-1
                  Z(I,J) = 0
               END DO
            END DO
            DO K=0,M*N-1,M
                  DO J=0,M-1
                     DO I=0,M-1
                        IF (UPLOSIGN.LT.4) THEN
                        Z(K+I,K+J) = A(I,J)
                     ELSE
                        Z(K+I,K+J) = A(J,I)
                     END IF
                  END DO
               END DO
            END DO
            DO K=0,M-1
                  DO J=0,N-1
                     DO I=0,N-1
                        IF (MOD(UPLOSIGN/2,2).EQ.0) THEN
                        Z(K+I*M,K+J*M) = Z(K+I*M,K+J*M) + SIGN*B(J,I)
                     ELSE
                        Z(K+I*M,K+J*M) = Z(K+I*M,K+J*M) + SIGN*B(I,J)
                     END IF
                  END DO
               END DO
            END DO
            CALL RECSY_SOLVE_CP(M*N,32,VEC,Z,SCALE, INFO, MACHINE)
            DO J=0,N-1
               DO I=0,M-1
                  C(I,J) = VEC(I+J*M)
               END DO
            END DO
         END IF
         RETURN
      ENDIF
      INFO = N+1
      IF (M/2.GE.N.OR.N.LE.4) THEN ! Split M dimension
         M1 = M/2
         M2 = M - M1
            IF (A(M1,M1-1).NE.0D0) THEN
            M1 = M1+1
            M2 = M2-1
         END IF
         SELECT CASE (UPLOSIGN/2)
         CASE (0,1) ! (N,N) and (N,T)
            CALL RECSYCT(UPLOSIGN, SCALE, M2, N, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_GEMM(0, M1, N, M2, -1D0, A(0,M1), LDA, C(M1,0), LDC, SCALE, C(0,0), LDC, MACHINE)
            CALL RECSYCT(UPLOSIGN, SC, M1, N, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N, SC, C(M1,0), LDC)
               SCALE = SCALE * SC
            END IF
         CASE (2,3) ! (T,N) and (T,T)
            CALL RECSYCT(UPLOSIGN, SCALE, M1, N, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_GEMM(1, M2, N, M1, -1D0, A(0,M1), LDA, C(0,0), LDC, SCALE, C(M1,0), LDC, MACHINE)
            CALL RECSYCT(UPLOSIGN, SC, M2, N, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N, SC, C(0,0), LDC)
               SCALE = SCALE * SC
            END IF
         END SELECT
      ELSE IF (N/2.GE.M.OR.M.LE.4) THEN ! Split N dimension
         N1 = N/2
         N2 = N - N1
            IF (B(N1,N1-1).NE.0D0) THEN
            N1 = N1+1
            N2 = N2-1
         END IF
         SELECT CASE (UPLOSIGN/2)
         CASE (0,2) ! (N,N) and (T,N)
            CALL RECSYCT(UPLOSIGN, SCALE, M, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_GEMM(0, M, N2, N1, -SIGN, C(0,0), LDC, B(0,N1), LDB, SCALE, C(0,N1), LDC, MACHINE)
            CALL RECSYCT(UPLOSIGN, SC, M, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,N1), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, C(0,0), LDC)
               SCALE = SCALE * SC
            END IF
         CASE (1,3) ! (N,T) and (T,T)
            CALL RECSYCT(UPLOSIGN, SCALE, M, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,N1), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            CALL RECSY_GEMM(2, M, N1, N2, -SIGN, C(0,N1), LDC, B(0,N1), LDB, SCALE, C(0,0), LDC, MACHINE)
            CALL RECSYCT(UPLOSIGN, SC, M, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC, C(0,N1), LDC)
               SCALE = SCALE * SC
            END IF
         END SELECT
      ELSE
         M1 = M/2
         M2 = M - M1
         IF (A(M1,M1-1).NE.0D0) THEN
            M1 = M1+1
            M2 = M2-1
         END IF
         N1 = N/2
         N2 = N - N1
         IF (B(N1,N1-1).NE.0D0) THEN
            N1 = N1+1
            N2 = N2-1
         END IF
         SELECT CASE (UPLOSIGN/2)
         CASE (0) ! (N,N)
            CALL RECSYCT(UPLOSIGN, SCALE, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_GEMM(0, M1, N1, M2, -1D0, A(0,M1), LDA, C(M1,0), LDC, SCALE, C(0,0), LDC, MACHINE)
            CALL RECSYCT(UPLOSIGN, SC, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N1, SC, C(M1,0), LDC)
               SCALE = SCALE * SC
            END IF
            CALL RECSY_GEMM(0, M2, N2, N1, -SIGN, C(M1,0), LDC, B(0,N1), LDB, SCALE, C(M1,N1), LDC, MACHINE)
            CALL RECSYCT(UPLOSIGN, SC, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,N1), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, C(0,0), LDC)
               SCALE = SCALE * SC
            END IF
            CALL RECSY_GEMM(0, M1, N2, M2, -1D0, A(0,M1), LDA, C(M1,N1), LDC, SCALE, C(0,N1), LDC, MACHINE)
            CALL RECSY_GEMM(0, M1, N2, N1, -SIGN, C(0,0), LDC, B(0,N1), LDB, +1D0, C(0,N1), LDC, MACHINE)
            CALL RECSYCT(UPLOSIGN, SC, M1, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,N1), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, C(0,0), LDC)
               CALL RECSY_SCAL(M2, N2, SC, C(M1,N1), LDC)
               SCALE = SCALE * SC
            END IF
         CASE (1) ! (N,T)
            CALL RECSYCT(UPLOSIGN, SCALE, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,N1), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            CALL RECSY_GEMM(0, M1, N2, M2, -1D0, A(0,M1), LDA, C(M1,N1), LDC, SCALE, C(0,N1), LDC, MACHINE)
            CALL RECSYCT(UPLOSIGN, SC, M1, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,N1), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N2, SC, C(M1,N1), LDC)   
               SCALE = SCALE * SC
            END IF
            CALL RECSY_GEMM(2, M2, N1, N2, -SIGN, C(M1,N1), LDC, B(0,N1), LDB, SCALE, C(M1,0), LDC, MACHINE)
            CALL RECSYCT(UPLOSIGN, SC, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC, C(0,N1), LDC)
               SCALE = SCALE * SC
            END IF
            CALL RECSY_GEMM(0, M1, N1, M2, -1D0, A(0,M1), LDA, C(M1,0), LDC, SCALE, C(0,0), LDC, MACHINE)
            CALL RECSY_GEMM(2, M1, N1, N2, -SIGN, C(0,N1), LDC, B(0,N1), LDB, +1D0, C(0,0), LDC, MACHINE)
            CALL RECSYCT(UPLOSIGN, SC, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC, C(0,N1), LDC)
               CALL RECSY_SCAL(M2, N1, SC, C(M1,0), LDC)
               SCALE = SCALE * SC
            END IF
         CASE (2) ! (T,N)
            CALL RECSYCT(UPLOSIGN, SCALE, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_GEMM(1, M2, N1, M1, -1D0, A(0,M1), LDA, C(0,0), LDC, SCALE, C(M1,0), LDC, MACHINE)
            CALL RECSYCT(UPLOSIGN, SC, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N1, SC, C(0,0), LDC)
               SCALE = SCALE * SC
            END IF
            CALL RECSY_GEMM(0, M1, N2, N1, -SIGN, C(0,0), LDC, B(0,N1), LDB, SCALE, C(0,N1), LDC, MACHINE)
            CALL RECSYCT(UPLOSIGN, SC, M1, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,N1), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, C(0,0), LDC)
               SCALE = SCALE * SC
            END IF
            CALL RECSY_GEMM(1, M2, N2, M1, -1D0, A(0,M1), LDA, C(0,N1), LDC, SCALE, C(M1,N1), LDC, MACHINE)
            CALL RECSY_GEMM(0, M2, N2, N1, -SIGN, C(M1,0), LDC, B(0,N1), LDB, +1D0, C(M1,N1), LDC, MACHINE)
            CALL RECSYCT(UPLOSIGN, SC, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,N1), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, C(0,0), LDC)
               CALL RECSY_SCAL(M1, N2, SC, C(0,N1), LDC)
               SCALE = SCALE * SC
            END IF
         CASE (3) ! (T,T)
            CALL RECSYCT(UPLOSIGN, SCALE, M1, N2, A(0,0), LDA, B(N1,N1),LDB, C(0,N1), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            CALL RECSY_GEMM(1, M2, N2, M1, -1D0, A(0,M1), LDA, C(0,N1), LDC, SCALE, C(M1,N1), LDC, MACHINE)
            CALL RECSYCT(UPLOSIGN, SC, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,N1), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N2, SC, C(0,N1), LDC)
               SCALE = SCALE * SC
            END IF
            CALL RECSY_GEMM(2, M1, N1, N2, -SIGN, C(0,N1), LDC, B(0,N1), LDB, SCALE, C(0,0), LDC, MACHINE)
            CALL RECSYCT(UPLOSIGN, SC, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC, C(0,N1), LDC)
               SCALE = SCALE * SC
            END IF
            CALL RECSY_GEMM(1, M2, N1, M1, -1D0, A(0,M1), LDA, C(0,0), LDC, SCALE, C(M1,0), LDC, MACHINE)
            CALL RECSY_GEMM(2, M2, N1, N2, -SIGN, C(M1,N1), LDC, B(0,N1), LDB, +1D0, C(M1,0), LDC, MACHINE)
            CALL RECSYCT(UPLOSIGN, SC, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC, C(0,N1), LDC)
               CALL RECSY_SCAL(M1, N1, SC, C(0,0), LDC)
               SCALE = SCALE * SC
            END IF
         END SELECT
      END IF
      IF (INFO.GT.N) INFO = 0
      RETURN
      END SUBROUTINE RECSYCT

