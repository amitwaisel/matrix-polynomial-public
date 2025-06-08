!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      RECURSIVE SUBROUTINE RECGLYDT(UPLO, SCALE, M, A, LDA, E, LDE, C, LDC, INFO, MACHINE, WORKSPACE, WKSIZE )
!     Solves op(A)*X*op(A^T) + op(E)*X*op(E^T) = scale*C
!     A (M*M) quasi-triangular, E (M*M) triangular
!     
!     Iff close eigenvalues are detected, INFO is non-zero and elements 
!     of the C (M*M) matrix are no longer valid
!
      IMPLICIT NONE
!
      INTENT(IN) M,A,LDA,E,LDE,LDC,UPLO,MACHINE
      INTENT(OUT) INFO,SCALE
      INTENT(INOUT) C, WORKSPACE, WKSIZE
!     .. Scalar Arguments ..
      INTEGER          :: M, LDA, LDC, LDE, INFO, UPLO, I, J
      INTEGER, OPTIONAL :: WKSIZE
      DOUBLE PRECISION :: SCALE
!     .. Array Arguments ..
      DOUBLE PRECISION :: A(0:LDA-1,0:M-1), E(0:LDE-1,0:M-1), C(0:LDC-1,0:M-1)
      DOUBLE PRECISION, OPTIONAL :: MACHINE(0:*), WORKSPACE(0:*)
!     .. Local Scalars ..
      INTEGER          :: M1, M2, INF, WORKNEEDED, WKS
      DOUBLE PRECISION :: SC
!     .. Local Arrays ..
      DOUBLE PRECISION :: TMPA(0:2,0:2)
      DOUBLE PRECISION :: TMPC(0:2), MACHINETEMP(0:9)
      DOUBLE PRECISION, ALLOCATABLE :: WORKTEMP(:)
!     .. Intrinsic Functions ..
      INTRINSIC           MOD
!     .. Explicit Declarations
      INTERFACE
         RECURSIVE SUBROUTINE RECGSYL(UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, INFO, MACHINE, WORKSPACE, WKSIZE )
            INTENT(IN) M,N,A,LDA,B,LDB,C,LDC,D,LDD,LDE,UPLOSIGN,MACHINE
            INTENT(OUT) INFO,SCALE
            INTENT(INOUT) E,WORKSPACE,WKSIZE
            INTEGER          :: M, N, LDA, LDB, LDC, LDD, LDE, INFO, UPLOSIGN
            INTEGER, OPTIONAL :: WKSIZE
            DOUBLE PRECISION :: SCALE
            DOUBLE PRECISION :: A(0:LDA-1,0:M-1), B(0:LDB-1,0:N-1), C(0:LDC-1,0:M-1), D(0:LDD-1,0:N-1), E(0:LDE-1,0:N-1)
            DOUBLE PRECISION, OPTIONAL :: WORKSPACE(0:*), MACHINE(0:*)
         END SUBROUTINE RECGSYL
      END INTERFACE
!     .. Executable Statements ..
!
      IF (.NOT.PRESENT(WKSIZE)) THEN
         WKS = 0
      ELSE
         WKS = WKSIZE
      END IF
      WORKNEEDED = M*M/2
      IF (WKS.LT.WORKNEEDED.OR..NOT.PRESENT(MACHINE)) THEN
         IF (WKS.EQ.-1) THEN
            WKSIZE = WORKNEEDED
            RETURN
         END IF
         IF (PRESENT(MACHINE)) THEN
            DO I = 0,9
               MACHINETEMP(I) = MACHINE(I)
            END DO
         ELSE
            CALL RECSY_MACHINE(MACHINETEMP)
         END IF
         IF (WKS.LT.WORKNEEDED) THEN
            ALLOCATE(WORKTEMP(WORKNEEDED), STAT=I)
            IF (I.NE.0) THEN
               INFO = -100
               RETURN
            END IF
            CALL RECGLYDT(UPLO, SCALE, M, A, LDA, E, LDE, C, LDC, INFO, MACHINETEMP, WORKTEMP, WORKNEEDED )
            DEALLOCATE(WORKTEMP)
         ELSE
            CALL RECGLYDT(UPLO, SCALE, M, A, LDA, E, LDE, C, LDC, INFO, MACHINETEMP, WORKSPACE, WKS )
         END IF
         RETURN
      END IF

      SCALE = 1D0
      IF (M.EQ.1.OR.(M.EQ.2.AND.A(1,0).NE.0D0)) THEN
         IF (M.EQ.1) THEN  
            TMPC(0) = C(0,0)
            CALL RECSY_SOLVE_11(C(0,0), A(0,0)*A(0,0)-E(0,0)*E(0,0), INFO, MACHINE)
            IF (INFO.NE.0) THEN
               C(0,0) = TMPC(0)
               CALL RECSY_SOLVE_CP(1, 32, C(0,0), A(0,0)*A(0,0)-E(0,0)*E(0,0), SCALE, INFO, MACHINE)
            END IF
         ELSE
            IF (MOD(UPLO,2).EQ.0) THEN
               TMPA(0,0) = A(0,0)*A(0,0) - E(0,0)*E(0,0)
               TMPA(1,0) = A(0,0)*A(1,0)
               TMPA(2,0) = A(1,0)*A(1,0)
               TMPA(0,1) = 2*A(0,0)*A(0,1) - 2*E(0,0)*E(0,1)
               TMPA(1,1) = A(0,0)*A(1,1)+A(0,1)*A(1,0) - E(0,0)*E(1,1)
               TMPA(2,1) = 2*A(1,1)*A(1,0)
               TMPA(0,2) = A(0,1)*A(0,1) - E(0,1)*E(0,1)
               TMPA(1,2) = A(1,1)*A(0,1) - E(1,1)*E(0,1)
               TMPA(2,2) = A(1,1)*A(1,1) - E(1,1)*E(1,1)
            ELSE
               TMPA(0,0) = A(0,0)*A(0,0) - E(0,0)*E(0,0)
               TMPA(1,0) = A(0,0)*A(0,1) - E(0,0)*E(0,1)
               TMPA(2,0) = A(0,1)*A(0,1) - E(0,1)*E(0,1)
               TMPA(0,1) = 2*A(0,0)*A(1,0)
               TMPA(1,1) = A(0,0)*A(1,1)+A(1,0)*A(0,1) - E(0,0)*E(1,1)
               TMPA(2,1) = 2*A(1,1)*A(0,1) - 2*E(1,1)*E(0,1)
               TMPA(0,2) = A(1,0)*A(1,0)
               TMPA(1,2) = A(1,1)*A(1,0)
               TMPA(2,2) = A(1,1)*A(1,1) - E(1,1)*E(1,1)
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
                  TMPA(0,0) = A(0,0)*A(0,0) - E(0,0)*E(0,0)
                  TMPA(1,0) = A(0,0)*A(1,0)
                  TMPA(2,0) = A(1,0)*A(1,0)
                  TMPA(0,1) = 2*A(0,0)*A(0,1) - 2*E(0,0)*E(0,1)
                  TMPA(1,1) = A(0,0)*A(1,1)+A(0,1)*A(1,0) - E(0,0)*E(1,1)
                  TMPA(2,1) = 2*A(1,1)*A(1,0)
                  TMPA(0,2) = A(0,1)*A(0,1) - E(0,1)*E(0,1)
                  TMPA(1,2) = A(1,1)*A(0,1) - E(1,1)*E(0,1)
                  TMPA(2,2) = A(1,1)*A(1,1) - E(1,1)*E(1,1)
               ELSE
                  TMPA(0,0) = A(0,0)*A(0,0) - E(0,0)*E(0,0)
                  TMPA(1,0) = A(0,0)*A(0,1) - E(0,0)*E(0,1)
                  TMPA(2,0) = A(0,1)*A(0,1) - E(0,1)*E(0,1)
                  TMPA(0,1) = 2*A(0,0)*A(1,0)
                  TMPA(1,1) = A(0,0)*A(1,1)+A(1,0)*A(0,1) - E(0,0)*E(1,1)
                  TMPA(2,1) = 2*A(1,1)*A(0,1) - 2*E(1,1)*E(0,1)
                  TMPA(0,2) = A(1,0)*A(1,0)
                  TMPA(1,2) = A(1,1)*A(1,0)
                  TMPA(2,2) = A(1,1)*A(1,1) - E(1,1)*E(1,1)
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

      M1 = M/2
      M2 = M - M1
      IF (A(M1,M1-1).NE.0D0) THEN
         M1 = M1+1
         M2 = M2-1
      END IF
      INFO = M+1
      SELECT CASE (UPLO)
      CASE (0)
         CALL RECGLYDT(UPLO, SCALE, M2, A(M1,M1), LDA, E(M1,M1), LDE, C(M1,M1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
         IF (INF.NE.0.AND.INF+M1.LT.INFO) INFO = INF+M1
         DO I=0,M-1
            DO J=I+1,M-1
              C(J,I)=C(I,J)
            END DO
         END DO
         CALL RECSY_AXB('F','V','N','T', M1, M2, M2, M2, -1D0, A(0,M1), LDA, C(M1,M1), LDC, A(M1,M1), LDA, SCALE, C(0,M1), LDC, MACHINE, WORKSPACE)
         CALL RECSY_AXB('F','U','N','T', M1, M2, M2, M2, +1D0, E(0,M1), LDE, C(M1,M1), LDC, E(M1,M1), LDE,   1D0, C(0,M1), LDC, MACHINE, WORKSPACE)
         CALL RECGSYL(3, SC, M1, M2, A(0,0), LDA, A(M1,M1), LDA, E(0,0), LDE, E(M1,M1), LDE, C(0,M1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
         IF (INF.NE.0.AND.INF+M1.LT.INFO) INFO = INF+M1
         IF (SC.NE.1D0) THEN
            CALL RECSY_SCAL(M2, M2, SC, C(M1,M1), LDC)
            SCALE = SCALE * SC
         END IF

         CALL MB01RD('U', 'N', M1, M2, SCALE, -1D0, C(0,0), LDC, A(0,M1), LDA, C(M1,M1), LDC, WORKSPACE, M1*M2, INFO)
         DO I=M1,M-1
            C(I,I) = C(I,I) * 2
         END DO
         CALL MB01RD('U', 'N', M1, M2,   1D0, +1D0, C(0,0), LDC, E(0,M1), LDE, C(M1,M1), LDC, WORKSPACE, M1*M2, INFO)
         DO I=M1,M-1
            C(I,I) = C(I,I) * 2
         END DO

         CALL RECSY_MULT_LEFT('V','N', M1, M2, M1, 1D0, A(0,0), LDA, C(0,M1), LDC, 0D0, WORKSPACE, M1, MACHINE, WORKSPACE(M1*M2))
         CALL DSYR2K('U', 'N', M1, M2, -1D0, A(0,M1), LDA, WORKSPACE, M1, 1D0, C(0,0), LDC)
         CALL RECSY_MULT_LEFT('U','N', M1, M2, M1, 1D0, E(0,0), LDE, C(0,M1), LDC, 0D0, WORKSPACE, M1, MACHINE, WORKSPACE(M1*M2))
         CALL DSYR2K('U', 'N', M1, M2, +1D0, E(0,M1), LDE, WORKSPACE, M1, 1D0, C(0,0), LDC)
         CALL RECGLYDT(UPLO, SC, M1, A(0,0), LDA, E(0,0), LDE, C(0,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
         IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
         IF (SC.NE.1D0) THEN
            CALL RECSY_SCAL(M, M2, SC, C(0,M1), LDC)
            SCALE = SCALE * SC
         END IF
      CASE (1)
         CALL RECGLYDT(UPLO, SCALE, M1, A(0,0), LDA, E(0,0), LDE, C(0,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
         IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
         DO I=0,M-1
            DO J=I+1,M-1
              C(J,I)=C(I,J)
            END DO
         END DO
         CALL RECSY_AXB('V','F','T','N', M1, M1, M1, M2, -1D0, A(0,0), LDA, C(0,0), LDC, A(0,M1), LDA, SCALE, C(0,M1), LDC, MACHINE, WORKSPACE)
         CALL RECSY_AXB('U','F','T','N', M1, M1, M1, M2, +1D0, E(0,0), LDE, C(0,0), LDC, E(0,M1), LDE,   1D0, C(0,M1), LDC, MACHINE, WORKSPACE)
         CALL RECGSYL(5, SC, M1, M2, A(0,0), LDA, A(M1,M1), LDA, E(0,0), LDE, E(M1,M1), LDE, C(0,M1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
         IF (INF.NE.0.AND.INF+M1.LT.INFO) INFO = INF+M1
         IF (SC.NE.1D0) THEN
            CALL RECSY_SCAL(M1, M1, SC, C(0,0), LDC)
            SCALE = SCALE * SC
         END IF

         CALL MB01RD('U', 'T', M2, M1, SCALE, -1D0, C(M1,M1), LDC, A(0,M1), LDA, C(0,0), LDC, WORKSPACE, M1*M2, INFO)
         DO I=0,M1-1
            C(I,I) = C(I,I) * 2
         END DO
         CALL MB01RD('U', 'T', M2, M1,   1D0, +1D0, C(M1,M1), LDC, E(0,M1), LDE, C(0,0), LDC, WORKSPACE, M1*M2, INFO)
         DO I=0,M1-1
            C(I,I) = C(I,I) * 2
         END DO

         CALL RECSY_MULT_RIGHT('V','N', M1, M2, M2, 1D0, C(0,M1), LDC, A(M1,M1), LDA, 0D0, WORKSPACE, M1, MACHINE, WORKSPACE(M1*M2))
         CALL DSYR2K('U', 'T', M2, M1, -1D0, A(0,M1), LDA, WORKSPACE, M1, 1D0, C(M1,M1), LDC)
         CALL RECSY_MULT_RIGHT('U','N', M1, M2, M2, 1D0, C(0,M1), LDC, E(M1,M1), LDE, 0D0, WORKSPACE, M1, MACHINE, WORKSPACE(M1*M2))
         CALL DSYR2K('U', 'T', M2, M1, +1D0, E(0,M1), LDE, WORKSPACE, M1, 1D0, C(M1,M1), LDC)
         CALL RECGLYDT(UPLO, SC, M2, A(M1,M1), LDA, E(M1,M1), LDE, C(M1,M1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
         IF (INF.NE.0.AND.INF+M1.LT.INFO) INFO = INF+M1
         IF (SC.NE.1D0) THEN
            CALL RECSY_SCAL(M1, M, SC, C(0,0), LDC)
            SCALE = SCALE * SC
         END IF
      END SELECT
      IF (INFO.GT.M) INFO = 0
      RETURN
      END SUBROUTINE RECGLYDT
