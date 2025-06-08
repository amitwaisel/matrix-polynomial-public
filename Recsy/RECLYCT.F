!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      RECURSIVE SUBROUTINE RECLYCT(UPLO, SCALE, M, A, LDA, C, LDC, INFO, MACHINE )
!     Solves op(A)*X + X*op(A) = scale*C
!     A (M*M) quasi-triangular
!     
!     Iff close eigenvalues are detected, INFO is non-zero and elements 
!     of the C (M*M) symmetric matrix are no longer valid
!
      IMPLICIT NONE
!
      INTENT(IN) M,A,LDA,LDC,UPLO,MACHINE
      INTENT(OUT) INFO,SCALE
      INTENT(INOUT) C
!     .. Scalar Arguments ..
      INTEGER          :: M, LDA, LDC, INFO, UPLO
      DOUBLE PRECISION :: SCALE
!     .. Array Arguments ..
      DOUBLE PRECISION :: A(0:LDA-1,0:*), C(0:LDC-1,0:*)
      DOUBLE PRECISION, OPTIONAL :: MACHINE(0:*)
!     .. Local Scalars ..
      INTEGER          :: M1, M2, INF
      DOUBLE PRECISION :: SC
!     .. Local Arrays ..
      DOUBLE PRECISION :: TMPA(0:2,0:2)
      DOUBLE PRECISION :: TMPC(0:2), MACHINETEMP(0:9)
!     .. Intrinsic Functions ..
      INTRINSIC           MOD
!     .. Explicit Declarations
      INTERFACE
         RECURSIVE SUBROUTINE RECSYCT(UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE )
            INTENT(IN) M,N,A,LDA,B,LDB,LDC,UPLOSIGN,MACHINE
            INTENT(OUT) INFO,SCALE
            INTENT(INOUT) C
            INTEGER          :: M, N, LDA, LDB, LDC, INFO, UPLOSIGN
            DOUBLE PRECISION :: SCALE
            DOUBLE PRECISION :: A(0:LDA-1,0:*), B(0:LDB-1,0:*), C(0:LDC-1,0:*)
            DOUBLE PRECISION, OPTIONAL :: MACHINE(0:*)
         END SUBROUTINE RECSYCT
      END INTERFACE
!     .. Executable Statements ..
!

      IF (.NOT.PRESENT(MACHINE)) THEN
         CALL RECSY_MACHINE(MACHINETEMP)
         CALL RECLYCT(UPLO, SCALE, M, A, LDA, C, LDC, INFO, MACHINETEMP )
         RETURN
      END IF

      SCALE = 1D0
      IF (M.EQ.1.OR.(M.EQ.2.AND.A(1,0).NE.0D0)) THEN
         IF (M.EQ.1) THEN  
            TMPC(0) = C(0,0)
            CALL RECSY_SOLVE_11(C(0,0), 2*A(0,0), INFO, MACHINE)
            IF (INFO.NE.0) THEN
               C(0,0) = TMPC(0)
               CALL RECSY_SOLVE_CP(1, 32, C(0,0), 2*A(0,0), SCALE, INFO, MACHINE)
            END IF
         ELSE
            IF (MOD(UPLO,2).EQ.0) THEN
               TMPA(0,0) = 2*A(0,0)
               TMPA(1,0) = A(1,0)
               TMPA(2,0) = 0
               TMPA(0,1) = 2*A(0,1)
               TMPA(1,1) = A(0,0)+A(1,1)
               TMPA(2,1) = 2*A(1,0)
               TMPA(0,2) = 0
               TMPA(1,2) = A(0,1)
               TMPA(2,2) = 2*A(1,1)
            ELSE
               TMPA(0,0) = 2*A(0,0)
               TMPA(1,0) = A(0,1)
               TMPA(2,0) = 0
               TMPA(0,1) = 2*A(1,0)
               TMPA(1,1) = A(0,0)+A(1,1)
               TMPA(2,1) = 2*A(0,1)
               TMPA(0,2) = 0
               TMPA(1,2) = A(1,0)
               TMPA(2,2) = 2*A(1,1)
            END IF
            TMPC(0) = C(0,0)
            IF (UPLO/2.EQ.0) THEN
               TMPC(1) = C(0,1)
            ELSE
               TMPC(1) = C(1,0)
            END IF
            TMPC(2) = C(1,1)
            CALL RECSY_SOLVE_33(TMPC, TMPA, INFO, MACHINE)
            INFO = 0
            IF (INFO.EQ.0) THEN
               C(0,0) = TMPC(0)
               IF (UPLO/2.EQ.0) THEN
                  C(0,1) = TMPC(1)
               ELSE
                  C(1,0) = TMPC(1)
               END IF
               C(1,1) = TMPC(2)
            ELSE
               IF (MOD(UPLO,2).EQ.0) THEN
                  TMPA(0,0) = 2*A(0,0)
                  TMPA(1,0) = A(1,0)
                  TMPA(2,0) = 0
                  TMPA(0,1) = 2*A(0,1)
                  TMPA(1,1) = A(0,0)+A(1,1)
                  TMPA(2,1) = 2*A(1,0)
                  TMPA(0,2) = 0
                  TMPA(1,2) = A(0,1)
                  TMPA(2,2) = 2*A(1,1)
               ELSE
                  TMPA(0,0) = 2*A(0,0)
                  TMPA(1,0) = A(0,1)
                  TMPA(2,0) = 0
                  TMPA(0,1) = 2*A(1,0)
                  TMPA(1,1) = A(0,0)+A(1,1)
                  TMPA(2,1) = 2*A(0,1)
                  TMPA(0,2) = 0
                  TMPA(1,2) = A(1,0)
                  TMPA(2,2) = 2*A(1,1)
               END IF
               TMPC(0) = C(0,0)
               IF (UPLO/2.EQ.0) THEN
                  TMPC(1) = C(0,1)
               ELSE
                  TMPC(1) = C(1,0)
               END IF
               TMPC(2) = C(1,1)
               CALL RECSY_SOLVE_CP(3, 3, TMPC, TMPA, SCALE, INFO, MACHINE)
               C(0,0) = TMPC(0)
               IF (UPLO/2.EQ.0) THEN
                  C(0,1) = TMPC(1)
               ELSE
                  C(1,0) = TMPC(1)
               END IF
               C(1,1) = TMPC(2)
            END IF
         END IF
         RETURN
      ENDIF
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
         CALL RECLYCT(UPLO, SCALE, M2, A(M1,M1), LDA, C(M1,M1), LDC, INF, MACHINE)
         IF (INF.NE.0.AND.INF+M1.LT.INFO) INFO = INF+M1
         CALL DSYMM('R', 'U', M1, M2, -1D0, C(M1,M1), LDC, A(0,M1), LDA, SCALE, C(0,M1), LDC)
         CALL RECSYCT(2, SC, M1, M2, A(0,0), LDA, A(M1,M1), LDA, C(0,M1), LDC, INF, MACHINE)      
         IF (INF.NE.0.AND.INF+M1.LT.INFO) INFO = INF+M1
         IF (SC.NE.1D0) THEN
            CALL RECSY_SCAL(M2, M2, SC, C(M1,M1), LDC)
            SCALE = SCALE * SC
         END IF
         CALL DSYR2K('U', 'N', M1, M2, -1D0, A(0,M1), LDA, C(0,M1), LDC, SCALE, C(0,0), LDC)
         CALL RECLYCT(UPLO, SC, M1, A(0,0), LDA, C(0,0), LDC, INF, MACHINE)
         IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
         IF (SC.NE.1D0) THEN
            CALL RECSY_SCAL(M, M2, SC, C(0,M1), LDC)
            SCALE = SCALE * SC
         END IF
      CASE (1) ! (down,up) and (down,down)
         CALL RECLYCT(UPLO, SCALE, M1, A(0,0), LDA, C(0,0), LDC, INF, MACHINE)
         IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
         CALL DSYMM('L', 'U', M1, M2, -1D0, C(0,0), LDC, A(0,M1), LDA, SCALE, C(0,M1), LDC)
         CALL RECSYCT(4, SC, M1, M2, A(0,0), LDA, A(M1,M1), LDA, C(0,M1), LDC, INF, MACHINE)      
         IF (INF.NE.0.AND.INF+M1.LT.INFO) INFO = INF+M1
         IF (SC.NE.1D0) THEN
            CALL RECSY_SCAL(M1, M1, SC, C(0,0), LDC)
            SCALE = SCALE * SC
         END IF
         CALL DSYR2K('U', 'T', M2, M1, -1D0, A(0,M1), LDA, C(0,M1), LDC, SCALE, C(M1,M1), LDC)
         CALL RECLYCT(UPLO, SC, M2, A(M1,M1), LDA, C(M1,M1), LDC, INF, MACHINE)
         IF (INF.NE.0.AND.INF+M1.LT.INFO) INFO = INF+M1
         IF (SC.NE.1D0) THEN
            CALL RECSY_SCAL(M1, M, SC, C(0,0), LDC)
            SCALE = SCALE * SC
         END IF
      END SELECT
      IF (INFO.GT.M) INFO = 0
      RETURN
      END SUBROUTINE RECLYCT
