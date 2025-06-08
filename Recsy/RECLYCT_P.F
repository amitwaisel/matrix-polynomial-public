!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      RECURSIVE SUBROUTINE RECLYCT_P(PROCS, UPLO, SCALE, M, A, LDA, C, LDC, INFO, MACHINE )
!     Solves op(A)*X + sign*X*op(A) = scale*C
!     A (M*M) and B (N*N) quasi-triangular
!     
!     Iff close eigenvalues are detected, INFO is non-zero and elements 
!     of the C (M*N) matrix are no longer valid
!
      IMPLICIT NONE
!
      INTENT(IN) PROCS, M,A,LDA,LDC,UPLO,MACHINE
      INTENT(OUT) INFO,SCALE
      INTENT(INOUT) C
!     .. Scalar Arguments ..
      INTEGER          :: PROCS, M, LDA, LDC, INFO, UPLO
      DOUBLE PRECISION :: SCALE
!     .. Array Arguments ..
      DOUBLE PRECISION :: A(0:LDA-1,0:*), C(0:LDC-1,0:*)
      DOUBLE PRECISION, OPTIONAL :: MACHINE(0:*)
!     .. Local Scalars ..
      INTEGER          :: M1, M2, INF
      DOUBLE PRECISION :: SC
!     .. Local Arrays ..
      DOUBLE PRECISION :: TMPA(0:2,0:2), MACHINETEMP(0:9)
      DOUBLE PRECISION :: TMPC(0:2)
!     .. Intrinsic Functions ..
      INTRINSIC           MOD
!     .. Explicit Declarations
      INTERFACE
         RECURSIVE SUBROUTINE RECLYCT(UPLO, SCALE, M, A, LDA, C, LDC, INFO, MACHINE )
            INTENT(IN) M,A,LDA,LDC,UPLO,MACHINE
            INTENT(OUT) INFO,SCALE
            INTENT(INOUT) C
            INTEGER          :: M, LDA, LDC, INFO, UPLO
            DOUBLE PRECISION :: SCALE
            DOUBLE PRECISION :: A(0:LDA-1,0:*), C(0:LDC-1,0:*)
            DOUBLE PRECISION, OPTIONAL :: MACHINE(0:*)
         END SUBROUTINE RECLYCT
         RECURSIVE SUBROUTINE RECSYCT_P(PROCS, UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE )
            INTENT(IN) PROCS,M,N,A,LDA,B,LDB,LDC,UPLOSIGN,MACHINE
            INTENT(OUT) INFO,SCALE
            INTENT(INOUT) C
            INTEGER          :: PROCS, M, N, LDA, LDB, LDC, INFO, UPLOSIGN
            DOUBLE PRECISION :: SCALE
            DOUBLE PRECISION :: A(0:LDA-1,0:*), B(0:LDB-1,0:*), C(0:LDC-1,0:*)
            DOUBLE PRECISION, OPTIONAL :: MACHINE(0:*)
         END SUBROUTINE RECSYCT_P
      END INTERFACE
!     .. Executable Statements ..
!
      IF (.NOT.PRESENT(MACHINE)) THEN
         CALL RECSY_MACHINE(MACHINETEMP)
         CALL RECLYCT_P(PROCS, UPLO, SCALE, M, A, LDA, C, LDC, INFO, MACHINETEMP )
         RETURN
      END IF

      IF (M.LE.2.OR.M*M*M.LT.MACHINE(5)) THEN
         CALL RECLYCT(UPLO, SCALE, M, A, LDA, C, LDC, INFO, MACHINE)
         RETURN
      END IF
!
      M1 = M/2
      M2 = M - M1
      IF (A(M1,M1-1).NE.0D0) THEN
         M1 = M1+1
         M2 = M2-1
      END IF
      INFO = M+1
      SELECT CASE (UPLO)
      CASE (0)
         CALL RECLYCT_P(PROCS, UPLO, SCALE, M2, A(M1,M1), LDA, C(M1,M1), LDC, INF, MACHINE)
         IF (INF.NE.0.AND.INF+M1.LT.INFO) INFO = INF+M1
         CALL DSYMM('R', 'U', M1, M2, -1D0, C(M1,M1), LDC, A(0,M1), LDA, SCALE, C(0,M1), LDC)
         CALL RECSYCT_P(PROCS, 2, SC, M1, M2, A(0,0), LDA, A(M1,M1), LDA, C(0,M1), LDC, INF, MACHINE)      
         IF (INF.NE.0.AND.INF+M1.LT.INFO) INFO = INF+M1
         IF (SC.NE.1D0) THEN
            CALL RECSY_SCAL(M2, M2, SC, C(M1,M1), LDC)
            SCALE = SCALE * SC
         END IF
         CALL DSYR2K('U', 'N', M1, M2, -1D0, A(0,M1), LDA, C(0,M1), LDC, SCALE, C(0,0), LDC)
         CALL RECLYCT_P(PROCS, UPLO, SC, M1, A(0,0), LDA, C(0,0), LDC, INF, MACHINE)
         IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
         IF (SC.NE.1D0) THEN
            CALL RECSY_SCAL(M, M2, SC, C(0,M1), LDC)
            SCALE = SCALE * SC
         END IF
      CASE (1) ! (down,up) and (down,down)
         CALL RECLYCT_P(PROCS, UPLO, SCALE, M1, A(0,0), LDA, C(0,0), LDC, INF, MACHINE)
         IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
         CALL DSYMM('L', 'U', M1, M2, -1D0, C(0,0), LDC, A(0,M1), LDA, SCALE, C(0,M1), LDC)
         CALL RECSYCT_P(PROCS, 4, SC, M1, M2, A(0,0), LDA, A(M1,M1), LDA, C(0,M1), LDC, INF, MACHINE)      
         IF (INF.NE.0.AND.INF+M1.LT.INFO) INFO = INF+M1
         IF (SC.NE.1D0) THEN
            CALL RECSY_SCAL(M1, M1, SC, C(0,0), LDC)
            SCALE = SCALE * SC
         END IF
         CALL DSYR2K('U', 'T', M2, M1, -1D0, A(0,M1), LDA, C(0,M1), LDC, SCALE, C(M1,M1), LDC)
         CALL RECLYCT_P(PROCS, UPLO, SC, M2, A(M1,M1), LDA, C(M1,M1), LDC, INF, MACHINE)
         IF (INF.NE.0.AND.INF+M1.LT.INFO) INFO = INF+M1
         IF (SC.NE.1D0) THEN
            CALL RECSY_SCAL(M1, M, SC, C(0,0), LDC)
            SCALE = SCALE * SC
         END IF
      END SELECT
      IF (INFO.GT.M) INFO = 0
      RETURN
      END SUBROUTINE RECLYCT_P
