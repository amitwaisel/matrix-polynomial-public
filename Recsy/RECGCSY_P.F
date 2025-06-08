!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      RECURSIVE SUBROUTINE RECGCSY_P(PROCS, UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, INFO, MACHINE )
!     Solves op(A)*X + sign*Y*op(B) = scale*C
!     Solves op(D)*X + sign*Y*op(E) = scale*F
!     A (M*M) and B (N*N) quasi-triangular, D (M*M) and E (N*N) quasi-triangular, 
!     
!     Iff close eigenvalues are detected, INFO is non-zero and elements 
!     of the C (M*N) and F (M*N) matrices are no longer valid
!
      IMPLICIT NONE
!
      INTENT(IN) PROCS,M,N,A,LDA,B,LDB,LDC,D,LDD,E,LDE,LDF,UPLOSIGN,MACHINE
      INTENT(OUT) INFO,SCALE
      INTENT(INOUT) C, F
!     .. Scalar Arguments ..
      INTEGER          :: PROCS, M, N, LDA, LDB, LDC, LDD, LDE, LDF, INFO, UPLOSIGN
      DOUBLE PRECISION :: SCALE
!     .. Array Arguments ..
      DOUBLE PRECISION :: A(0:LDA-1,0:*), B(0:LDB-1,0:*), C(0:LDC-1,0:*)
      DOUBLE PRECISION :: D(0:LDD-1,0:*), E(0:LDE-1,0:*), F(0:LDF-1,0:*)
      DOUBLE PRECISION, OPTIONAL :: MACHINE(0:*)
!     .. Local Scalars ..
      INTEGER          :: I, J, K, M1, M2, N1, N2, INF
      DOUBLE PRECISION :: SIGN, SC, SC2
!     .. Local Arrays ..
      DOUBLE PRECISION :: MACHINETEMP(0:9)
!     .. Intrinsic Functions ..
      INTRINSIC           MIN, MOD
!     .. Explicit Declarations
      INTERFACE
         RECURSIVE SUBROUTINE RECGCSY(UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, INFO, MACHINE )
            INTENT(IN) M,N,A,LDA,B,LDB,LDC,D,LDD,E,LDE,LDF,UPLOSIGN,MACHINE
            INTENT(OUT) INFO,SCALE
            INTENT(INOUT) C, F
            INTEGER          :: M, N, LDA, LDB, LDC, LDD, LDE, LDF, INFO, UPLOSIGN
            DOUBLE PRECISION :: SCALE
            DOUBLE PRECISION :: A(0:LDA-1,0:*), B(0:LDB-1,0:*), C(0:LDC-1,0:*)
            DOUBLE PRECISION :: D(0:LDD-1,0:*), E(0:LDE-1,0:*), F(0:LDF-1,0:*)
            DOUBLE PRECISION, OPTIONAL :: MACHINE(0:*)
         END SUBROUTINE RECGCSY
      END INTERFACE
!     .. Executable Statements ..
!
      IF (.NOT.PRESENT(MACHINE)) THEN
         CALL RECSY_MACHINE(MACHINETEMP)
         CALL RECGCSY_P(PROCS, UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, INFO, MACHINETEMP )
         RETURN
      END IF

      IF (PROCS.LE.1.OR.M*N*1D0*(M+N).LE.MACHINE(5)) THEN
         CALL RECGCSY(UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, INFO, MACHINE)
         RETURN
      ENDIF

      SIGN = 1D0-(MOD(UPLOSIGN,2))*2D0
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
            CALL RECGCSY_P(PROCS, UPLOSIGN, SCALE, M2, N, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, D(M1,M1), LDD, E(0,0), LDE, F(M1,0), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_GEMM(0, M1, N, M2, -1D0, A(0,M1), LDA, C(M1,0), LDC, SCALE, C(0,0), LDC, MACHINE)
            CALL RECSY_GEMM(0, M1, N, M2, -1D0, D(0,M1), LDD, C(M1,0), LDC, SCALE, F(0,0), LDF, MACHINE)
            CALL RECGCSY_P(PROCS, UPLOSIGN, SC, M1, N, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, F(0,0), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N, SC, C(M1,0), LDC)
               CALL RECSY_SCAL(M2, N, SC, F(M1,0), LDF)
               SCALE = SCALE * SC
            END IF
         CASE (2,3) ! (T,N) and (T,T)
            CALL RECGCSY_P(PROCS, UPLOSIGN, SCALE, M1, N, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, F(0,0), LDF, INF, MACHINE)
            IF (INF.NE.0) INFO = INF
            CALL RECSY_GEMM(1, M2, N, M1, -1D0, A(0,M1), LDA, C(0,0), LDC, SCALE, C(M1,0), LDC, MACHINE)
            CALL RECSY_GEMM(1, M2, N, M1, -1D0, D(0,M1), LDD, C(0,0), LDC, SCALE, F(M1,0), LDF, MACHINE)
            CALL RECGCSY_P(PROCS, UPLOSIGN, SC, M2, N, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, D(M1,M1), LDD, E(0,0), LDE, F(M1,0), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N, SC, C(0,0), LDC)
               CALL RECSY_SCAL(M1, N, SC, F(0,0), LDF)
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
            CALL RECGCSY_P(PROCS, UPLOSIGN, SCALE, M, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, F(0,0), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_GEMM(0, M, N2, N1, -SIGN, F(0,0), LDF, B(0,N1), LDB, SCALE, C(0,N1), LDC, MACHINE)
            CALL RECSY_GEMM(0, M, N2, N1, -SIGN, F(0,0), LDF, E(0,N1), LDE, SCALE, F(0,N1), LDF, MACHINE)
            CALL RECGCSY_P(PROCS, UPLOSIGN, SC, M, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,N1), LDC, D(0,0), LDD, E(N1,N1), LDE, F(0,N1), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, C(0,0), LDC)
               CALL RECSY_SCAL(M, N1, SC, F(0,0), LDF)
               SCALE = SCALE * SC
            END IF
         CASE (1,3) ! (N,T) and (T,T)
            CALL RECGCSY_P(PROCS, UPLOSIGN, SCALE, M, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,N1), LDC, D(0,0), LDD, E(N1,N1), LDE, F(0,N1), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            CALL RECSY_GEMM(2, M, N1, N2, -SIGN, F(0,N1), LDF, B(0,N1), LDB, SCALE, C(0,0), LDC, MACHINE)
            CALL RECSY_GEMM(2, M, N1, N2, -SIGN, F(0,N1), LDF, E(0,N1), LDE, SCALE, F(0,0), LDF, MACHINE)
            CALL RECGCSY_P(PROCS, UPLOSIGN, SC, M, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, F(0,0), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC, C(0,N1), LDC)
               CALL RECSY_SCAL(M, N2, SC, F(0,N1), LDF)
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
            CALL RECGCSY_P(PROCS, UPLOSIGN, SCALE, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, D(M1,M1), LDD, E(0,0), LDE, F(M1,0), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
!$OMP  PARALLEL SECTIONS PRIVATE(INF) REDUCTION(MIN:INFO)
!$OMP  SECTION
            CALL RECSY_GEMM(0, M1, N1, M2, -1D0, A(0,M1), LDA, C(M1,0), LDC, SCALE, C(0,0), LDC, MACHINE)
            CALL RECSY_GEMM(0, M1, N1, M2, -1D0, D(0,M1), LDD, C(M1,0), LDC, SCALE, F(0,0), LDF, MACHINE)
            CALL RECGCSY_P(PROCS / 2, UPLOSIGN, SC, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, F(0,0), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
!$OMP  SECTION
            CALL RECSY_GEMM(0, M2, N2, N1, -SIGN, F(M1,0), LDF, B(0,N1), LDB, SCALE, C(M1,N1), LDC, MACHINE)
            CALL RECSY_GEMM(0, M2, N2, N1, -SIGN, F(M1,0), LDF, E(0,N1), LDE, SCALE, F(M1,N1), LDF, MACHINE)
            CALL RECGCSY_P(PROCS / 2, UPLOSIGN, SC2, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,N1), LDC, D(M1,M1), LDD, E(N1,N1), LDE, F(M1,N1), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
!$OMP  END PARALLEL SECTIONS
            IF (SC.NE.1D0.AND.SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N2, SC, C(M1,N1), LDC)
               CALL RECSY_SCAL(M2, N2, SC, F(M1,N1), LDF)
               CALL RECSY_SCAL(M1, N1, SC2, C(0,0), LDC)
               CALL RECSY_SCAL(M1, N1, SC2, F(0,0), LDF)
               CALL RECSY_SCAL(M2, N1, SC*SC2, C(M1,0), LDC)
               CALL RECSY_SCAL(M2, N1, SC*SC2, F(M1,0), LDF)
               SCALE = SCALE * SC * SC2
            ELSE IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N, SC, C(M1,0), LDC)
               CALL RECSY_SCAL(M2, N, SC, F(M1,0), LDF)
               SCALE = SCALE * SC
            ELSE IF (SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC2, C(0,0), LDC)
               CALL RECSY_SCAL(M, N1, SC2, F(0,0), LDF)
               SCALE = SCALE * SC2
            END IF
            CALL RECSY_GEMM(0, M1, N2, M2, -1D0, A(0,M1), LDA, C(M1,N1), LDC, SCALE, C(0,N1), LDC, MACHINE)
            CALL RECSY_GEMM(0, M1, N2, M2, -1D0, D(0,M1), LDD, C(M1,N1), LDC, SCALE, F(0,N1), LDF, MACHINE)
            CALL RECSY_GEMM(0, M1, N2, N1, -SIGN, F(0,0), LDF, B(0,N1), LDB, +1D0, C(0,N1), LDC, MACHINE)
            CALL RECSY_GEMM(0, M1, N2, N1, -SIGN, F(0,0), LDF, E(0,N1), LDE, +1D0, F(0,N1), LDF, MACHINE)
            CALL RECGCSY_P(PROCS, UPLOSIGN, SC, M1, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,N1), LDC, D(0,0), LDD, E(N1,N1), LDE, F(0,N1), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, C(0,0), LDC)
               CALL RECSY_SCAL(M2, N2, SC, C(M1,N1), LDC)
               CALL RECSY_SCAL(M, N1, SC, F(0,0), LDF)
               CALL RECSY_SCAL(M2, N2, SC, F(M1,N1), LDF)
               SCALE = SCALE * SC
            END IF
         CASE (1) ! (N,T)
            CALL RECGCSY_P(PROCS, UPLOSIGN, SCALE, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,N1), LDC, D(M1,M1), LDD, E(N1,N1), LDE, F(M1,N1), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
!$OMP  PARALLEL SECTIONS PRIVATE(INF) REDUCTION(MIN:INFO)
!$OMP  SECTION
            CALL RECSY_GEMM(0, M1, N2, M2, -1D0, A(0,M1), LDA, C(M1,N1), LDC, SCALE, C(0,N1), LDC, MACHINE)
            CALL RECSY_GEMM(0, M1, N2, M2, -1D0, D(0,M1), LDD, C(M1,N1), LDC, SCALE, F(0,N1), LDF, MACHINE)
            CALL RECGCSY_P(PROCS / 2, UPLOSIGN, SC, M1, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,N1), LDC, D(0,0), LDD, E(N1,N1), LDE, F(0,N1), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
!$OMP  SECTION
            CALL RECSY_GEMM(2, M2, N1, N2, -SIGN, F(M1,N1), LDF, B(0,N1), LDB, SCALE, C(M1,0), LDC, MACHINE)
            CALL RECSY_GEMM(2, M2, N1, N2, -SIGN, F(M1,N1), LDF, E(0,N1), LDE, SCALE, F(M1,0), LDF, MACHINE)
            CALL RECGCSY_P(PROCS / 2, UPLOSIGN, SC2, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, D(M1,M1), LDD, E(0,0), LDE, F(M1,0), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
!$OMP  END PARALLEL SECTIONS
            IF (SC.NE.1D0.AND.SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N1, SC, C(M1,0), LDC)   
               CALL RECSY_SCAL(M2, N1, SC, F(M1,0), LDF)   
               CALL RECSY_SCAL(M1, N2, SC2, C(0,N1), LDC)
               CALL RECSY_SCAL(M1, N2, SC2, F(0,N1), LDF)
               CALL RECSY_SCAL(M2, N2, SC*SC2, C(M1,N1), LDC)
               CALL RECSY_SCAL(M2, N2, SC*SC2, F(M1,N1), LDF)
               SCALE = SCALE * SC * SC2
            ELSE IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N, SC, C(M1,0), LDC)   
               CALL RECSY_SCAL(M2, N, SC, F(M1,0), LDF)   
               SCALE = SCALE * SC
            ELSE IF (SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC2, C(0,N1), LDC)
               CALL RECSY_SCAL(M, N2, SC2, F(0,N1), LDF)
               SCALE = SCALE * SC2
            END IF
            CALL RECSY_GEMM(0, M1, N1, M2, -1D0, A(0,M1), LDA, C(M1,0), LDC, SCALE, C(0,0), LDC, MACHINE)
            CALL RECSY_GEMM(0, M1, N1, M2, -1D0, D(0,M1), LDD, C(M1,0), LDF, SCALE, F(0,0), LDF, MACHINE)
            CALL RECSY_GEMM(2, M1, N1, N2, -SIGN, F(0,N1), LDF, B(0,N1), LDB, +1D0, C(0,0), LDC, MACHINE)
            CALL RECSY_GEMM(2, M1, N1, N2, -SIGN, F(0,N1), LDF, E(0,N1), LDE, +1D0, F(0,0), LDF, MACHINE)
            CALL RECGCSY_P(PROCS, UPLOSIGN, SC, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, F(0,0), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC, C(0,N1), LDC)
               CALL RECSY_SCAL(M2, N1, SC, C(M1,0), LDC)
               CALL RECSY_SCAL(M, N2, SC, F(0,N1), LDF)
               CALL RECSY_SCAL(M2, N1, SC, F(M1,0), LDF)
               SCALE = SCALE * SC
            END IF
         CASE (2) ! (T,N)
            CALL RECGCSY_P(PROCS, UPLOSIGN, SCALE, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, F(0,0), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
!$OMP  PARALLEL SECTIONS PRIVATE(INF) REDUCTION(MIN:INFO)
!$OMP  SECTION
            CALL RECSY_GEMM(1, M2, N1, M1, -1D0, A(0,M1), LDA, C(0,0), LDC, SCALE, C(M1,0), LDC, MACHINE)
            CALL RECSY_GEMM(1, M2, N1, M1, -1D0, D(0,M1), LDD, C(0,0), LDC, SCALE, F(M1,0), LDF, MACHINE)
            CALL RECGCSY_P(PROCS / 2, UPLOSIGN, SC, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, D(M1,M1), LDD, E(0,0), LDE, F(M1,0), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
!$OMP  SECTION
            CALL RECSY_GEMM(0, M1, N2, N1, -SIGN, F(0,0), LDF, B(0,N1), LDB, SCALE, C(0,N1), LDC, MACHINE)
            CALL RECSY_GEMM(0, M1, N2, N1, -SIGN, F(0,0), LDF, E(0,N1), LDE, SCALE, F(0,N1), LDF, MACHINE)
            CALL RECGCSY_P(PROCS / 2, UPLOSIGN, SC2, M1, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,N1), LDC, D(0,0), LDD, E(N1,N1), LDE, F(0,N1), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
!$OMP  END PARALLEL SECTIONS
            IF (SC.NE.1D0.AND.SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N2, SC, C(0,N1), LDC)
               CALL RECSY_SCAL(M1, N2, SC, F(0,N1), LDF)
               CALL RECSY_SCAL(M2, N1, SC2, C(M1,0), LDC)
               CALL RECSY_SCAL(M2, N1, SC2, F(M1,0), LDF)
               CALL RECSY_SCAL(M1, N1, SC*SC2, C(0,0), LDC)
               CALL RECSY_SCAL(M1, N1, SC*SC2, F(0,0), LDF)
               SCALE = SCALE * SC * SC2
            ELSE IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N, SC, C(0,0), LDC)
               CALL RECSY_SCAL(M1, N, SC, F(0,0), LDF)
               SCALE = SCALE * SC
            ELSE IF (SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC2, C(0,0), LDC)
               CALL RECSY_SCAL(M, N1, SC2, F(0,0), LDF)
               SCALE = SCALE * SC2
            END IF
            CALL RECSY_GEMM(1, M2, N2, M1, -1D0, A(0,M1), LDA, C(0,N1), LDC, SCALE, C(M1,N1), LDC, MACHINE)
            CALL RECSY_GEMM(1, M2, N2, M1, -1D0, D(0,M1), LDD, C(0,N1), LDC, SCALE, F(M1,N1), LDF, MACHINE)
            CALL RECSY_GEMM(0, M2, N2, N1, -SIGN, F(M1,0), LDF, B(0,N1), LDB, +1D0, C(M1,N1), LDC, MACHINE)
            CALL RECSY_GEMM(0, M2, N2, N1, -SIGN, F(M1,0), LDF, E(0,N1), LDE, +1D0, F(M1,N1), LDF, MACHINE)
            CALL RECGCSY_P(PROCS, UPLOSIGN, SC, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,N1), LDC, D(M1,M1), LDD, E(N1,N1), LDE, F(M1,N1), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, C(0,0), LDC)
               CALL RECSY_SCAL(M1, N2, SC, C(0,N1), LDC)
               CALL RECSY_SCAL(M, N1, SC, F(0,0), LDF)
               CALL RECSY_SCAL(M1, N2, SC, F(0,N1), LDF)
               SCALE = SCALE * SC
            END IF
         CASE (3) ! (T,T)
            CALL RECGCSY_P(PROCS, UPLOSIGN, SCALE, M1, N2, A(0,0), LDA, B(N1,N1),LDB, C(0,N1), LDC, D(0,0), LDD, E(N1,N1), LDE, F(0,N1), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
!$OMP  PARALLEL SECTIONS PRIVATE(INF) REDUCTION(MIN:INFO)
!$OMP  SECTION
            CALL RECSY_GEMM(1, M2, N2, M1, -1D0, A(0,M1), LDA, C(0,N1), LDC, SCALE, C(M1,N1), LDC, MACHINE)
            CALL RECSY_GEMM(1, M2, N2, M1, -1D0, D(0,M1), LDD, C(0,N1), LDC, SCALE, F(M1,N1), LDF, MACHINE)
            CALL RECGCSY_P(PROCS / 2, UPLOSIGN, SC, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,N1), LDC, D(M1,M1), LDD, E(N1,N1), LDE, F(M1,N1), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
!$OMP  SECTION
            CALL RECSY_GEMM(2, M1, N1, N2, -SIGN, F(0,N1), LDF, B(0,N1), LDB, SCALE, C(0,0), LDC, MACHINE)
            CALL RECSY_GEMM(2, M1, N1, N2, -SIGN, F(0,N1), LDF, E(0,N1), LDE, SCALE, F(0,0), LDF, MACHINE)
            CALL RECGCSY_P(PROCS / 2, UPLOSIGN, SC2, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, F(0,0), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
!$OMP  END PARALLEL SECTIONS
            IF (SC.NE.1D0.AND.SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N1, SC, C(0,0), LDC)
               CALL RECSY_SCAL(M1, N1, SC, F(0,0), LDF)
               CALL RECSY_SCAL(M2, N2, SC2, C(M1,N1), LDC)
               CALL RECSY_SCAL(M2, N2, SC2, F(M1,N1), LDF)
               CALL RECSY_SCAL(M1, N2, SC*SC2, C(0,N1), LDC)
               CALL RECSY_SCAL(M1, N2, SC*SC2, F(0,N1), LDF)
               SCALE = SCALE * SC * SC2
            ELSE IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N, SC, C(0,0), LDC)
               CALL RECSY_SCAL(M1, N, SC, F(0,0), LDF)
               SCALE = SCALE * SC
            ELSE IF (SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC2, C(0,N1), LDC)
               CALL RECSY_SCAL(M, N2, SC2, F(0,N1), LDF)
               SCALE = SCALE * SC2
            END IF
            CALL RECSY_GEMM(1, M2, N1, M1, -1D0, A(0,M1), LDA, C(0,0), LDC, SCALE, C(M1,0), LDC, MACHINE)
            CALL RECSY_GEMM(1, M2, N1, M1, -1D0, D(0,M1), LDD, C(0,0), LDC, SCALE, F(M1,0), LDF, MACHINE)
            CALL RECSY_GEMM(2, M2, N1, N2, -SIGN, F(M1,N1), LDF, B(0,N1), LDB, +1D0, C(M1,0), LDC, MACHINE)
            CALL RECSY_GEMM(2, M2, N1, N2, -SIGN, F(M1,N1), LDF, E(0,N1), LDE, +1D0, F(M1,0), LDF, MACHINE)
            CALL RECGCSY_P(PROCS, UPLOSIGN, SC, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, D(M1,M1), LDD, E(0,0), LDE, F(M1,0), LDF, INF, MACHINE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC, C(0,N1), LDC)
               CALL RECSY_SCAL(M1, N1, SC, C(0,0), LDC)
               CALL RECSY_SCAL(M, N2, SC, F(0,N1), LDF)
               CALL RECSY_SCAL(M1, N1, SC, F(0,0), LDF)
               SCALE = SCALE * SC
            END IF
         END SELECT
      END IF
      IF (INFO.GT.N) INFO = 0
      RETURN
      END SUBROUTINE RECGCSY_P
