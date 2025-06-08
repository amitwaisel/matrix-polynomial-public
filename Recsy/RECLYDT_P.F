!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      RECURSIVE SUBROUTINE RECLYDT_P(PROCS, UPLO, SCALE, M, A, LDA, C, LDC, INFO, MACHINE, WORKSPACE, WKSIZE )
!     Solves op(A)*X*op(B) = scale*C
!     A (M*M) and B (N*N) quasi-triangular
!     
!     Iff close eigenvalues are detected, INFO is non-zero and elements 
!     of the C (M*N) matrix are no longer valid
!
      IMPLICIT NONE
!
      INTENT(IN) PROCS, M,A,LDA,LDC,UPLO,MACHINE
      INTENT(OUT) INFO,SCALE
      INTENT(INOUT) C, WORKSPACE,WKSIZE
!     .. Scalar Arguments ..
      INTEGER          :: PROCS, M, LDA, LDC, INFO, UPLO
      INTEGER, OPTIONAL :: WKSIZE
      DOUBLE PRECISION :: SCALE
!     .. Array Arguments ..
      DOUBLE PRECISION :: A(0:LDA-1,0:M-1), C(0:LDC-1,0:M-1)
      DOUBLE PRECISION, OPTIONAL :: MACHINE(0:*), WORKSPACE(0:*)
!     .. Local Scalars ..
      INTEGER          :: I, J, M1, M2, INF, WORKNEEDED, WKS
      DOUBLE PRECISION :: SC
!     .. Local Arrays ..
      DOUBLE PRECISION :: TMPA(0:2,0:2)
      DOUBLE PRECISION :: TMPC(0:2), MACHINETEMP(0:9)
      DOUBLE PRECISION, ALLOCATABLE :: WORKTEMP(:)
!     .. Intrinsic Functions ..
      INTRINSIC           MIN, MOD
!     .. Explicit Declarations
      INTERFACE
         RECURSIVE SUBROUTINE RECLYDT(UPLO, SCALE, M, A, LDA, C, LDC, INFO, MACHINE, WORKSPACE, WKSIZE )
            INTENT(IN) M,A,LDA,LDC,UPLO,MACHINE
            INTENT(OUT) INFO,SCALE
            INTENT(INOUT) C, WORKSPACE,WKSIZE
            INTEGER          :: M, LDA, LDC, INFO, UPLO
            INTEGER, OPTIONAL :: WKSIZE
            DOUBLE PRECISION :: SCALE
            DOUBLE PRECISION :: A(0:LDA-1,0:M-1), C(0:LDC-1,0:M-1)
            DOUBLE PRECISION, OPTIONAL :: MACHINE(0:*), WORKSPACE(0:*)
         END SUBROUTINE RECLYDT
         RECURSIVE SUBROUTINE RECSYDT_P(PROCS, UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE, WORKSPACE, WKSIZE )
            INTENT(IN) PROCS,M,N,A,LDA,B,LDB,LDC,UPLOSIGN,MACHINE
            INTENT(OUT) INFO,SCALE
            INTENT(INOUT) C,WORKSPACE,WKSIZE
            INTEGER          :: PROCS, M, N, LDA, LDB, LDC, INFO, UPLOSIGN
            INTEGER, OPTIONAL :: WKSIZE
            DOUBLE PRECISION :: SCALE
            DOUBLE PRECISION :: A(0:LDA-1,0:M-1), B(0:LDB-1,0:N-1), C(0:LDC-1,0:N-1)
            DOUBLE PRECISION, OPTIONAL :: WORKSPACE(0:*), MACHINE(0:*)
         END SUBROUTINE RECSYDT_P
      END INTERFACE
!     .. Executable Statements ..
!
      IF (.NOT.PRESENT(WKSIZE)) THEN
         WKS = 0
      ELSE
         WKS = WKSIZE
      END IF
      WORKNEEDED = M*M
      IF (WKS.LT.WORKNEEDED.OR..NOT.PRESENT(MACHINE)) THEN
         IF (WKS.EQ.-1) THEN
            WKSIZE = WORKNEEDED
            RETURN
         END IF
         IF (PRESENT(MACHINE)) THEN
            DO M1 = 0,9
               MACHINETEMP(M1) = MACHINE(M1)
            END DO
         ELSE
            CALL RECSY_MACHINE(MACHINETEMP)
         END IF
         IF (WKS.LT.WORKNEEDED) THEN
            ALLOCATE(WORKTEMP(WORKNEEDED), STAT=M1)
            IF (M1.NE.0) THEN
               INFO = -100
               RETURN
            END IF
            CALL RECLYDT_P(PROCS, UPLO, SCALE, M, A, LDA, C, LDC, INFO, MACHINETEMP, WORKTEMP, WORKNEEDED )
            DEALLOCATE(WORKTEMP)
         ELSE
            CALL RECLYDT_P(PROCS, UPLO, SCALE, M, A, LDA, C, LDC, INFO, MACHINETEMP, WORKSPACE, WKSIZE )
         END IF
         RETURN
      END IF

      IF (M.LE.2.OR.M*M*M.LT.MACHINE(5)) THEN
         CALL RECLYDT(UPLO, SCALE, M, A, LDA, C, LDC, INFO, MACHINE, WORKSPACE, WKSIZE)
         RETURN
      END IF

      M1 = M/2
      M2 = M - M1
      IF (A(M1,M1-1).NE.0D0) THEN
         M1 = M1+1
         M2 = M2-1
      END IF
      INFO = M+1
      SELECT CASE (UPLO)
      CASE (0)
         CALL RECLYDT_P(PROCS, UPLO, SCALE, M2, A(M1,M1), LDA, C(M1,M1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
         IF (INF.NE.0.AND.INF+M1.LT.INFO) INFO = INF+M1
         DO I=0,M-1
            DO J=I+1,M-1
              C(J,I)=C(I,J)
            END DO
         END DO
         CALL RECSY_AXB('F','V','N','T', M1, M2, M2, M2, -1D0, A(0,M1), LDA, C(M1,M1), LDC, A(M1,M1), LDA, SCALE, C(0,M1), LDC, MACHINE, WORKSPACE)
         CALL RECSYDT_P(PROCS, 3, SC, M1, M2, A(0,0), LDA, A(M1,M1), LDA, C(0,M1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
         IF (INF.NE.0.AND.INF+M1.LT.INFO) INFO = INF+M1
         IF (SC.NE.1D0) THEN
            CALL RECSY_SCAL(M2, M2, SC, C(M1,M1), LDC)
            SCALE = SCALE * SC
         END IF

         CALL MB01RD('U', 'N', M1, M2, SCALE, -1D0, C(0,0), LDC, A(0,M1), LDA, C(M1,M1), LDC, WORKSPACE, M1*M2, INFO)
         DO I=M1,M-1
            C(I,I) = C(I,I) * 2
         END DO

         CALL RECSY_MULT_LEFT('V','N', M1, M2, M1, 1D0, A(0,0), LDA, C(0,M1), LDC, 0D0, WORKSPACE, M1, MACHINE, WORKSPACE(M1*M2))
         CALL DSYR2K('U', 'N', M1, M2, -1D0, A(0,M1), LDA, WORKSPACE, M1, 1D0, C(0,0), LDC)
         CALL RECLYDT_P(PROCS, UPLO, SC, M1, A(0,0), LDA, C(0,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
         IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
         IF (SC.NE.1D0) THEN
            CALL RECSY_SCAL(M, M2, SC, C(0,M1), LDC)
            SCALE = SCALE * SC
         END IF
      CASE (1)
         CALL RECLYDT_P(PROCS, UPLO, SCALE, M1, A(0,0), LDA, C(0,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
         IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
         DO I=0,M-1
            DO J=I+1,M-1
              C(J,I)=C(I,J)
            END DO
         END DO
         CALL RECSY_AXB('V','F','T','N', M1, M1, M1, M2, -1D0, A(0,0), LDA, C(0,0), LDC, A(0,M1), LDA, SCALE, C(0,M1), LDC, MACHINE, WORKSPACE)
         CALL RECSYDT_P(PROCS, 5, SC, M1, M2, A(0,0), LDA, A(M1,M1), LDA, C(0,M1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
         IF (INF.NE.0.AND.INF+M1.LT.INFO) INFO = INF+M1
         IF (SC.NE.1D0) THEN
            CALL RECSY_SCAL(M1, M1, SC, C(0,0), LDC)
            SCALE = SCALE * SC
         END IF

         CALL MB01RD('U', 'T', M2, M1, SCALE, -1D0, C(M1,M1), LDC, A(0,M1), LDA, C(0,0), LDC, WORKSPACE, M1*M2, INFO)
         DO I=0,M1-1
            C(I,I) = C(I,I) * 2
         END DO

         CALL RECSY_MULT_RIGHT('V','N', M1, M2, M2, 1D0, C(0,M1), LDC, A(M1,M1), LDA, 0D0, WORKSPACE, M1, MACHINE, WORKSPACE(M1*M2))
         CALL DSYR2K('U', 'T', M2, M1, -1D0, A(0,M1), LDA, WORKSPACE, M1, 1D0, C(M1,M1), LDC)
         CALL RECLYDT_P(PROCS, UPLO, SC, M2, A(M1,M1), LDA, C(M1,M1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
         IF (INF.NE.0.AND.INF+M1.LT.INFO) INFO = INF+M1
         IF (SC.NE.1D0) THEN
            CALL RECSY_SCAL(M1, M, SC, C(0,0), LDC)
            SCALE = SCALE * SC
         END IF
      END SELECT
      IF (INFO.GT.M) INFO = 0
      RETURN
      END SUBROUTINE RECLYDT_P
