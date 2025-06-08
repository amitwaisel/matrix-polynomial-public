!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      RECURSIVE SUBROUTINE RECGSYL_P(PROCS, UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, INFO, MACHINE, WORKSPACE, WKSIZE )
!     Solves op(A)*X*op(B) + sign*op(C)*X*op(D) = scale*E
!     A (M*M) and B (N*N) quasi-triangular, C and D triangular
!     
!     Iff close eigenvalues are detected, INFO is non-zero and elements 
!     of the E (M*N) matrix are no longer valid
!
      IMPLICIT NONE
!
      INTENT(IN) PROCS,M,N,A,LDA,B,LDB,C,LDC,D,LDD,LDE,UPLOSIGN,MACHINE
      INTENT(OUT) INFO,SCALE
      INTENT(INOUT) E,WORKSPACE,WKSIZE

!     .. Scalar Arguments ..
      INTEGER          :: PROCS, M, N, LDA, LDB, LDC, LDD, LDE, INFO, UPLOSIGN
      INTEGER, OPTIONAL :: WKSIZE
      DOUBLE PRECISION :: SCALE
!     .. Array Arguments ..
      DOUBLE PRECISION :: A(0:LDA-1,0:M-1), B(0:LDB-1,0:N-1), C(0:LDC-1,0:M-1), D(0:LDD-1,0:N-1), E(0:LDE-1,0:N-1)
      DOUBLE PRECISION, OPTIONAL :: WORKSPACE(0:*), MACHINE(0:*)
!     .. Local Scalars ..
      INTEGER          :: I, J, K, L, M1, M2, N1, N2, INF, WORKNEEDED, WKS
      DOUBLE PRECISION :: SIGN, SC, SC2, TMPA, TMPB, TMPC, TMPD
      CHARACTER*1      :: TA, TB
!     .. Local Arrays ..
      DOUBLE PRECISION :: Z(0:31,0:31)
      DOUBLE PRECISION :: VEC(0:15), MACHINETEMP(0:9)
      DOUBLE PRECISION, ALLOCATABLE :: WORKTEMP(:)
!     .. Intrinsic Functions ..
      INTRINSIC           MOD, MIN
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
      WORKNEEDED = M*N+4
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
            CALL RECGSYL_P(PROCS, UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, INFO, MACHINE, WORKTEMP, WORKNEEDED )
            DEALLOCATE(WORKTEMP)
         ELSE
            CALL RECGSYL_P(PROCS, UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, INFO, MACHINE, WORKSPACE, WKSIZE )
         END IF
         RETURN
      END IF

      IF (PROCS.LE.1.OR.M*N*1D0*(M+N).LE.MACHINE(5)) THEN
         CALL RECGSYL(UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, INFO, MACHINE, WORKSPACE, WKSIZE)
         RETURN
      ENDIF

      SIGN = 1D0-(MOD(UPLOSIGN,2))*2D0
      IF (MOD(UPLOSIGN,8)/4.EQ.0) THEN
         TA = 'N'
      ELSE
         TA = 'T'
      END IF
      IF (MOD(UPLOSIGN/2,2).EQ.0) THEN
         TB = 'N'
      ELSE
         TB = 'T'
      END IF
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
            CALL RECGSYL_P(PROCS, UPLOSIGN, SCALE, M2, N, A(M1,M1), LDA, B(0,0), LDB, C(M1,M1), LDC, D(0,0), LDD, E(M1,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_AXB('F','V','N',TB, M1, M2, N, N, -1D0, A(0,M1), LDA, E(M1,0), LDE, B, LDB, SCALE, E(0,0), LDC, MACHINE, WORKSPACE)
            CALL RECSY_AXB('F','U','N',TB, M1, M2, N, N, -SIGN,C(0,M1), LDC, E(M1,0), LDE, D, LDD, 1D0,   E(0,0), LDC, MACHINE, WORKSPACE)
            CALL RECGSYL_P(PROCS, UPLOSIGN, SC, M1, N, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N, SC, E(M1,0), LDE)
               SCALE = SCALE * SC
            END IF
         CASE (2,3) ! (T,N) and (T,T)
            CALL RECGSYL_P(PROCS, UPLOSIGN, SCALE, M1, N, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_AXB('F','V','T',TB, M2, M1, N, N, -1D0, A(0,M1), LDA, E(0,0), LDE, B, LDB, SCALE, E(M1,0), LDE, MACHINE, WORKSPACE)
            CALL RECSY_AXB('F','U','T',TB, M2, M1, N, N, -SIGN,C(0,M1), LDC, E(0,0), LDE, D, LDD, 1D0,   E(M1,0), LDE, MACHINE, WORKSPACE)
            CALL RECGSYL_P(PROCS, UPLOSIGN, SC, M2, N, A(M1,M1), LDA, B(0,0), LDB, C(M1,M1), LDC, D(0,0), LDD, E(M1,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N, SC, E(0,0), LDE)
               SCALE = SCALE * SC
            END IF
         CASE (5) 
            CALL RECGSYL_P(PROCS, UPLOSIGN, SCALE, M2, N, A(M1,M1), LDA, B(0,0), LDB, C(M1,M1), LDC, D(0,0), LDD, E(M1,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_AXB('F','U','N',TB, M1, M2, N, N, -1D0, A(0,M1), LDA, E(M1,0), LDE, B, LDB, SCALE, E(0,0), LDC, MACHINE, WORKSPACE)
            CALL RECSY_AXB('F','V','N',TB, M1, M2, N, N, -SIGN,C(0,M1), LDC, E(M1,0), LDE, D, LDD, 1D0,   E(0,0), LDC, MACHINE, WORKSPACE)
            CALL RECGSYL_P(PROCS, UPLOSIGN, SC, M1, N, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N, SC, E(M1,0), LDE)
               SCALE = SCALE * SC
            END IF
         CASE (6)
            CALL RECGSYL_P(PROCS, UPLOSIGN, SCALE, M1, N, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_AXB('F','U','T',TB, M2, M1, N, N, -1D0, A(0,M1), LDA, E(0,0), LDE, B, LDB, SCALE, E(M1,0), LDE, MACHINE, WORKSPACE)
            CALL RECSY_AXB('F','V','T',TB, M2, M1, N, N, -SIGN,C(0,M1), LDC, E(0,0), LDE, D, LDD, 1D0,   E(M1,0), LDE, MACHINE, WORKSPACE)
            CALL RECGSYL_P(PROCS, UPLOSIGN, SC, M2, N, A(M1,M1), LDA, B(0,0), LDB, C(M1,M1), LDC, D(0,0), LDD, E(M1,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N, SC, E(0,0), LDE)
               SCALE = SCALE * SC
            END IF
         END SELECT
      ELSE IF (N/2.GE.M.OR.M.LE.4) THEN ! Split N dimension
         N1 = N/2
         N2 = N - N1
         IF (UPLOSIGN.GE.8) THEN
            IF (D(N1,N1-1).NE.0D0) THEN
               N1 = N1+1
               N2 = N2-1
            END IF
         ELSE
            IF (B(N1,N1-1).NE.0D0) THEN
               N1 = N1+1
               N2 = N2-1
            END IF
         END IF
         SELECT CASE (UPLOSIGN/2)
         CASE (0,2,6) ! (N,N) and (T,N)
            CALL RECGSYL_P(PROCS, UPLOSIGN, SCALE, M, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_AXB('V','F',TA,'N', M, M, N1, N2, -1D0, A(0,0), LDA, E(0,0), LDE, B(0,N1), LDB, SCALE, E(0,N1), LDE, MACHINE, WORKSPACE)
            CALL RECSY_AXB('U','F',TA,'N', M, M, N1, N2, -SIGN,C(0,0), LDC, E(0,0), LDE, D(0,N1), LDD, 1D0,   E(0,N1), LDE, MACHINE, WORKSPACE)
            CALL RECGSYL_P(PROCS, UPLOSIGN, SC, M, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,0), LDC, D(N1,N1), LDD, E(0,N1), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, E(0,0), LDE)
               SCALE = SCALE * SC
            END IF
         CASE (1,3,5) ! (N,T) and (T,T)
            CALL RECGSYL_P(PROCS, UPLOSIGN, SCALE, M, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,0), LDC, D(N1,N1), LDD, E(0,N1), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            CALL RECSY_AXB('V','F',TA,'T', M, M, N2, N1, -1D0, A(0,0), LDA, E(0,N1), LDE, B(0,N1), LDB, SCALE, E(0,0), LDE, MACHINE, WORKSPACE)
            CALL RECSY_AXB('U','F',TA,'T', M, M, N2, N1, -SIGN,C(0,0), LDC, E(0,N1), LDE, D(0,N1), LDD, 1D0,   E(0,0), LDE, MACHINE, WORKSPACE)
            CALL RECGSYL_P(PROCS, UPLOSIGN, SC, M, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC, E(0,N1), LDE)
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
         IF (UPLOSIGN.GE.8) THEN
            IF (D(N1,N1-1).NE.0D0) THEN
               N1 = N1+1
               N2 = N2-1
            END IF
         ELSE
            IF (B(N1,N1-1).NE.0D0) THEN
               N1 = N1+1
               N2 = N2-1
            END IF
         END IF
         SELECT CASE (UPLOSIGN/2)
         CASE (0) ! (N,N)
            CALL RECGSYL_P(PROCS, UPLOSIGN, SCALE, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,M1), LDC, D(0,0), LDD, E(M1,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
!$OMP  PARALLEL SECTIONS PRIVATE(INF) REDUCTION(MIN:INFO)
!$OMP  SECTION
            CALL RECSY_AXB('F','V','N','N', M1, M2, N1, N1, -1D0, A(0,M1), LDA, E(M1,0), LDE, B(0,0), LDB, SCALE, E(0,0), LDE, MACHINE, WORKSPACE)
            CALL RECSY_AXB('F','U','N','N', M1, M2, N1, N1,-SIGN, C(0,M1), LDC, E(M1,0), LDE, D(0,0), LDD,   1D0, E(0,0), LDE, MACHINE, WORKSPACE)
            CALL RECGSYL_P(PROCS/2, UPLOSIGN, SC, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
!$OMP  SECTION
            CALL RECSY_AXB('V','F','N','N', M2, M2, N1, N2, -1D0, A(M1,M1), LDA, E(M1,0), LDE, B(0,N1), LDB, SCALE, E(M1,N1), LDE, MACHINE, WORKSPACE(2*M1*N1))
            CALL RECSY_AXB('U','F','N','N', M2, M2, N1, N2,-SIGN, C(M1,M1), LDC, E(M1,0), LDE, D(0,N1), LDD,   1D0, E(M1,N1), LDE, MACHINE, WORKSPACE(2*M1*N1))
            CALL RECGSYL_P(PROCS/2, UPLOSIGN, SC2, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,M1), LDC, D(N1,N1), LDD, E(M1,N1), LDE, INF, MACHINE, WORKSPACE(2*M1*N1), WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
!$OMP  END PARALLEL SECTIONS
            IF (SC.NE.1D0.AND.SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N2, SC, E(M1,N1), LDE)
               CALL RECSY_SCAL(M1, N1, SC2, E(0,0), LDE)
               CALL RECSY_SCAL(M2, N1, SC*SC2, E(M1,0), LDE)
               SCALE = SCALE * SC * SC2
            ELSE IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N, SC, E(M1,0), LDE)
               SCALE = SCALE * SC
            ELSE IF (SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC2, E(0,0), LDE)
               SCALE = SCALE * SC2
            END IF

            IF (M.GT.N) THEN
               CALL RECSY_AXB('F','V','N','N', M1, M2,  N, N2, -1D0, A(0,M1), LDA, E(M1,0), LDE, B(0,N1), LDB, SCALE, E(0,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','U','N','N', M1, M2,  N, N2,-SIGN, C(0,M1), LDC, E(M1,0), LDE, D(0,N1), LDD, 1D0, E(0,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('V','F','N','N', M1, M1, N1, N2, -1D0, A(0,0), LDA, E(0,0), LDE, B(0,N1), LDB, 1D0, E(0,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('U','F','N','N', M1, M1, N1, N2,-SIGN, C(0,0), LDC, E(0,0), LDE, D(0,N1), LDD, 1D0, E(0,N1), LDE, MACHINE, WORKSPACE)
            ELSE
               CALL RECSY_AXB('V','F','N','N', M1,  M, N1, N2, -1D0, A(0,0), LDA, E(0,0), LDE, B(0,N1), LDB, SCALE, E(0,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('U','F','N','N', M1,  M, N1, N2,-SIGN, C(0,0), LDC, E(0,0), LDE, D(0,N1), LDD, 1D0, E(0,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','V','N','N', M1, M2, N2, N2, -1D0, A(0,M1), LDA, E(M1,N1), LDE, B(N1,N1), LDB, 1D0, E(0,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','U','N','N', M1, M2, N2, N2,-SIGN, C(0,M1), LDC, E(M1,N1), LDE, D(N1,N1), LDD, 1D0, E(0,N1), LDE, MACHINE, WORKSPACE)
            END IF

            CALL RECGSYL_P(PROCS, UPLOSIGN, SC, M1, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,0), LDC, D(N1,N1), LDD, E(0,N1), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, E(0,0), LDE)
               CALL RECSY_SCAL(M2, N2, SC, E(M1,N1), LDE)
               SCALE = SCALE * SC
            END IF
         CASE (1) ! (N,T)
            CALL RECGSYL_P(PROCS, UPLOSIGN, SCALE, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,M1), LDC, D(N1,N1), LDD, E(M1,N1), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
!$OMP  PARALLEL SECTIONS PRIVATE(INF) REDUCTION(MIN:INFO)
!$OMP  SECTION
            CALL RECSY_AXB('F','V','N','T', M1, M2, N2, N2, -1D0, A(0,M1), LDA, E(M1,N1), LDE, B(N1,N1), LDB, SCALE, E(0,N1), LDE, MACHINE, WORKSPACE)
            CALL RECSY_AXB('F','U','N','T', M1, M2, N2, N2,-SIGN, C(0,M1), LDC, E(M1,N1), LDE, D(N1,N1), LDD,   1D0, E(0,N1), LDE, MACHINE, WORKSPACE)
            CALL RECGSYL_P(PROCS/2, UPLOSIGN, SC, M1, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,0), LDC, D(N1,N1), LDD, E(0,N1), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
!$OMP  SECTION
            CALL RECSY_AXB('V','F','N','T', M2, M2, N2, N1, -1D0, A(M1,M1), LDA, E(M1,N1), LDE, B(0,N1), LDB, SCALE, E(M1,0), LDE, MACHINE, WORKSPACE(2*M1*N2))
            CALL RECSY_AXB('V','F','N','T', M2, M2, N2, N1,-SIGN, C(M1,M1), LDC, E(M1,N1), LDE, D(0,N1), LDD,   1D0, E(M1,0), LDE, MACHINE, WORKSPACE(2*M1*N2))
            CALL RECGSYL_P(PROCS/2, UPLOSIGN, SC2, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,M1), LDC, D(0,0), LDD, E(M1,0), LDE, INF, MACHINE, WORKSPACE(2*M1*N2), WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
!$OMP  END PARALLEL SECTIONS
            IF (SC.NE.1D0.AND.SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N1, SC, E(M1,0), LDE)   
               CALL RECSY_SCAL(M1, N2, SC2, E(0,N1), LDE)
               CALL RECSY_SCAL(M2, N2, SC*SC2, E(M1,N1), LDE)
               SCALE = SCALE * SC * SC2
            ELSE IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N, SC, E(M1,0), LDE)   
               SCALE = SCALE * SC
            ELSE IF (SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC2, E(0,N1), LDE)
               SCALE = SCALE * SC2
            END IF

            IF (M.GT.N) THEN
               CALL RECSY_AXB('F','V','N','T', M1, M2, N, N1, -1D0, A(0,M1), LDA, E(M1,0), LDE, B(0,0), LDB, SCALE, E(0,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','U','N','T', M1, M2, N, N1, -SIGN, C(0,M1), LDC, E(M1,0), LDE, D(0,0), LDD, 1D0, E(0,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('V','F','N','T', M1, M1, N2, N1, -1D0, A(0,0), LDA, E(0,N1), LDE, B(0,N1), LDB, 1D0, E(0,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('U','F','N','T', M1, M1, N2, N1, -SIGN, C(0,0), LDC, E(0,N1), LDE, D(0,N1), LDD, 1D0, E(0,0), LDE, MACHINE, WORKSPACE)
            ELSE
               CALL RECSY_AXB('V','F','N','T', M1, M, N2, N1, -1D0, A(0,0), LDA, E(0,N1), LDE, B(0,N1), LDB, SCALE, E(0,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('U','F','N','T', M1, M, N2, N1, -SIGN, C(0,0), LDC, E(0,N1), LDE, D(0,N1), LDD, 1D0, E(0,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','V','N','T', M1, M2, N1, N1, -1D0, A(0,M1), LDA, E(M1,0), LDE, B(0,0), LDB, 1D0, E(0,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','U','N','T', M1, M2, N1, N1, -SIGN, C(0,M1), LDA, E(M1,0), LDE, D(0,0), LDD, 1D0, E(0,0), LDE, MACHINE, WORKSPACE)
            END IF
            CALL RECGSYL_P(PROCS, UPLOSIGN, SC, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC, E(0,N1), LDE)
               CALL RECSY_SCAL(M2, N1, SC, E(M1,0), LDE)
               SCALE = SCALE * SC
            END IF
         CASE (2) ! (T,N)
            CALL RECGSYL_P(PROCS, UPLOSIGN, SCALE, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
!$OMP  PARALLEL SECTIONS PRIVATE(INF) REDUCTION(MIN:INFO)
!$OMP  SECTION
            CALL RECSY_AXB('F','V','T','N', M2, M1, N1, N1, -1D0, A(0,M1), LDA, E(0,0), LDE, B(0,0), LDB, SCALE, E(M1,0), LDE, MACHINE, WORKSPACE)
            CALL RECSY_AXB('F','U','T','N', M2, M1, N1, N1, -SIGN, C(0,M1), LDC, E(0,0), LDE, D(0,0), LDD,   1D0, E(M1,0), LDE, MACHINE, WORKSPACE)
            CALL RECGSYL_P(PROCS/2, UPLOSIGN, SC, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,M1), LDC, D(0,0), LDD, E(M1,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
!$OMP  SECTION
            CALL RECSY_AXB('V','F','T','N', M1, M1, N1, N2, -1D0, A(0,0), LDA, E(0,0), LDE, B(0,N1), LDB, SCALE, E(0,N1), LDE, MACHINE, WORKSPACE(2*M2*N1))
            CALL RECSY_AXB('U','F','T','N', M1, M1, N1, N2, -SIGN, C(0,0), LDC, E(0,0), LDE, D(0,N1), LDD,   1D0, E(0,N1), LDE, MACHINE, WORKSPACE(2*M2*N1))
            CALL RECGSYL_P(PROCS/2, UPLOSIGN, SC2, M1, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,0), LDC, D(N1,N1), LDD, E(0,N1), LDE, INF, MACHINE, WORKSPACE(2*M2*N1), WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
!$OMP  END PARALLEL SECTIONS
            IF (SC.NE.1D0.AND.SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N2, SC, E(0,N1), LDE)
               CALL RECSY_SCAL(M2, N1, SC2, E(M1,0), LDE)
               CALL RECSY_SCAL(M1, N1, SC*SC2, E(0,0), LDE)
               SCALE = SCALE * SC * SC2
            ELSE IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N, SC, E(0,0), LDE)
               SCALE = SCALE * SC
            ELSE IF (SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC2, E(0,0), LDE)
               SCALE = SCALE * SC2
            END IF
            IF (M.GT.N) THEN
               CALL RECSY_AXB('F','V','T','N', M2, M1, N, N2, -1D0, A(0,M1), LDA, E(0,0), LDE, B(0,N1), LDB, SCALE, E(M1,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','U','T','N', M2, M1, N, N2, -SIGN, C(0,M1), LDC, E(0,0), LDE, D(0,N1), LDD,   1D0, E(M1,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('V','F','T','N', M2, M2, N1, N2, -1D0, A(M1,M1), LDA, E(M1,0), LDE, B(0,N1), LDB, 1D0, E(M1,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('U','F','T','N', M2, M2, N1, N2, -SIGN, C(M1,M1), LDC, E(M1,0), LDE, D(0,N1), LDD, 1D0, E(M1,N1), LDE, MACHINE, WORKSPACE)
            ELSE
               CALL RECSY_AXB('V','F','T','N', M2, M, N1, N2, -1D0, A(0,M1), LDA, E(0,0), LDE, B(0,N1), LDB, SCALE, E(M1,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('U','F','T','N', M2, M, N1, N2, -SIGN, C(0,M1), LDC, E(0,0), LDE, D(0,N1), LDD,   1D0, E(M1,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','V','T','N', M2, M1, N2, N2, -1D0, A(0,M1), LDA, E(0,N1), LDE, B(N1,N1), LDB, 1D0, E(M1,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','U','T','N', M2, M1, N2, N2, -SIGN, C(0,M1), LDC, E(0,N1), LDE, D(N1,N1), LDD, 1D0, E(M1,N1), LDE, MACHINE, WORKSPACE)
            END IF
            CALL RECGSYL_P(PROCS, UPLOSIGN, SC, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,M1), LDC, D(N1,N1), LDD, E(M1,N1), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, E(0,0), LDE)
               CALL RECSY_SCAL(M1, N2, SC, E(0,N1), LDE)
               SCALE = SCALE * SC
            END IF
         CASE (3) ! (T,T)
            CALL RECGSYL_P(PROCS, UPLOSIGN, SCALE, M1, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,0), LDC, D(N1,N1), LDD, E(0,N1), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
!$OMP  PARALLEL SECTIONS PRIVATE(INF) REDUCTION(MIN:INFO)
!$OMP  SECTION
            CALL RECSY_AXB('F','V','T','T', M2, M1, N2, N2, -1D0, A(0,M1), LDA, E(0,N1), LDE, B(N1,N1), LDB, SCALE, E(M1,N1), LDE, MACHINE, WORKSPACE)
            CALL RECSY_AXB('F','U','T','T', M2, M1, N2, N2, -SIGN, C(0,M1), LDC, E(0,N1), LDE, D(N1,N1), LDD,   1D0, E(M1,N1), LDE, MACHINE, WORKSPACE)
            CALL RECGSYL_P(PROCS/2, UPLOSIGN, SC, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,M1), LDC, D(N1,N1), LDD, E(M1,N1), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
!$OMP  SECTION
            CALL RECSY_AXB('V','F','T','T', M1, M1, N2, N1, -1D0, A(0,0), LDA, E(0,N1), LDE, B(0,N1), LDB, SCALE, E(0,0), LDE, MACHINE, WORKSPACE(2*M2*N2))
            CALL RECSY_AXB('U','F','T','T', M1, M1, N2, N1, -SIGN, C(0,0), LDC, E(0,N1), LDE, D(0,N1), LDD,   1D0, E(0,0), LDE, MACHINE, WORKSPACE(2*M2*N2))
            CALL RECGSYL_P(PROCS/2, UPLOSIGN, SC2, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, INF, MACHINE, WORKSPACE(2*M2*N2), WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
!$OMP  END PARALLEL SECTIONS
            IF (SC.NE.1D0.AND.SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N1, SC, E(0,0), LDE)
               CALL RECSY_SCAL(M2, N2, SC2, E(M1,N1), LDE)
               CALL RECSY_SCAL(M1, N2, SC*SC2, E(0,N1), LDE)
               SCALE = SCALE * SC * SC2
            ELSE IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N, SC, E(0,0), LDE)
               SCALE = SCALE * SC
            ELSE IF (SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC2, E(0,N1), LDE)
               SCALE = SCALE * SC2
            END IF
            IF (M.GT.N) THEN
               CALL RECSY_AXB('F','V','T','T', M2, M1, N, N1, -1D0, A(0,M1), LDA, E(0,0), LDE, B(0,0), LDB, SCALE, E(M1,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','U','T','T', M2, M1, N, N1, -SIGN, C(0,M1), LDC, E(0,0), LDE, D(0,0), LDD,   1D0, E(M1,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('V','F','T','T', M2, M2, N2, N1, -1D0, A(M1,M1), LDA, E(M1,N1), LDE, B(0,N1), LDB, 1D0, E(M1,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('U','F','T','T', M2, M2, N2, N1, -SIGN, C(M1,M1), LDC, E(M1,N1), LDE, D(0,N1), LDD, 1D0, E(M1,0), LDE, MACHINE, WORKSPACE)
            ELSE
               CALL RECSY_AXB('F','V','T','T', M2, M1, N1, N1, -1D0, A(0,M1), LDA, E(0,0), LDE, B(0,0), LDB, SCALE, E(M1,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','U','T','T', M2, M1, N1, N1, -SIGN, C(0,M1), LDC, E(0,0), LDE, D(0,0), LDD,   1D0, E(M1,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('V','F','T','T', M2, M, N2, N1, -1D0, A(0,M1), LDA, E(0,N1), LDE, B(0,N1), LDB, 1D0, E(M1,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('U','F','T','T', M2, M, N2, N1, -SIGN, C(0,M1), LDC, E(0,N1), LDE, D(0,N1), LDD, 1D0, E(M1,0), LDE, MACHINE, WORKSPACE)
            END IF
            CALL RECGSYL_P(PROCS, UPLOSIGN, SC, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,M1), LDC, D(0,0), LDD, E(M1,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC, E(0,N1), LDE)
               CALL RECSY_SCAL(M1, N1, SC, E(0,0), LDE)
               SCALE = SCALE * SC
            END IF
         CASE (5) ! (N,T)
            CALL RECGSYL_P(PROCS, UPLOSIGN, SCALE, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,M1), LDC, D(N1,N1), LDD, E(M1,N1), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
!$OMP  PARALLEL SECTIONS PRIVATE(INF) REDUCTION(MIN:INFO)
!$OMP  SECTION
            CALL RECSY_AXB('F','U','N','T', M1, M2, N2, N2, -1D0, A(0,M1), LDA, E(M1,N1), LDE, B(N1,N1), LDB, SCALE, E(0,N1), LDE, MACHINE, WORKSPACE)
            CALL RECSY_AXB('F','V','N','T', M1, M2, N2, N2,-SIGN, C(0,M1), LDC, E(M1,N1), LDE, D(N1,N1), LDD,   1D0, E(0,N1), LDE, MACHINE, WORKSPACE)
            CALL RECGSYL_P(PROCS/2, UPLOSIGN, SC, M1, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,0), LDC, D(N1,N1), LDD, E(0,N1), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
!$OMP  SECTION
            CALL RECSY_AXB('V','F','N','T', M2, M2, N2, N1, -1D0, A(M1,M1), LDA, E(M1,N1), LDE, B(0,N1), LDB, SCALE, E(M1,0), LDE, MACHINE, WORKSPACE(2*M1*N2))
            CALL RECSY_AXB('V','F','N','T', M2, M2, N2, N1,-SIGN, C(M1,M1), LDC, E(M1,N1), LDE, D(0,N1), LDD,   1D0, E(M1,0), LDE, MACHINE, WORKSPACE(2*M1*N2))
            CALL RECGSYL_P(PROCS/2, UPLOSIGN, SC2, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,M1), LDC, D(0,0), LDD, E(M1,0), LDE, INF, MACHINE, WORKSPACE(2*M1*N2), WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
!$OMP  END PARALLEL SECTIONS
            IF (SC.NE.1D0.AND.SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N1, SC, E(M1,0), LDE)   
               CALL RECSY_SCAL(M1, N2, SC2, E(0,N1), LDE)
               CALL RECSY_SCAL(M2, N2, SC*SC2, E(M1,N1), LDE)
               SCALE = SCALE * SC * SC2
            ELSE IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N, SC, E(M1,0), LDE)   
               SCALE = SCALE * SC
            ELSE IF (SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC2, E(0,N1), LDE)
               SCALE = SCALE * SC2
            END IF

            IF (M.GT.N) THEN
               CALL RECSY_AXB('F','U','N','T', M1, M2, N, N1, -1D0, A(0,M1), LDA, E(M1,0), LDE, B(0,0), LDB, SCALE, E(0,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','V','N','T', M1, M2, N, N1, -SIGN, C(0,M1), LDC, E(M1,0), LDE, D(0,0), LDD, 1D0, E(0,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('V','F','N','T', M1, M1, N2, N1, -1D0, A(0,0), LDA, E(0,N1), LDE, B(0,N1), LDB, 1D0, E(0,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('U','F','N','T', M1, M1, N2, N1, -SIGN, C(0,0), LDC, E(0,N1), LDE, D(0,N1), LDD, 1D0, E(0,0), LDE, MACHINE, WORKSPACE)
            ELSE
               CALL RECSY_AXB('V','F','N','T', M1, M, N2, N1, -1D0, A(0,0), LDA, E(0,N1), LDE, B(0,N1), LDB, SCALE, E(0,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('U','F','N','T', M1, M, N2, N1, -SIGN, C(0,0), LDC, E(0,N1), LDE, D(0,N1), LDD, 1D0, E(0,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','U','N','T', M1, M2, N1, N1, -1D0, A(0,M1), LDA, E(M1,0), LDE, B(0,0), LDB, 1D0, E(0,0), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','V','N','T', M1, M2, N1, N1, -SIGN, C(0,M1), LDA, E(M1,0), LDE, D(0,0), LDD, 1D0, E(0,0), LDE, MACHINE, WORKSPACE)
            END IF
            CALL RECGSYL_P(PROCS, UPLOSIGN, SC, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC, E(0,N1), LDE)
               CALL RECSY_SCAL(M2, N1, SC, E(M1,0), LDE)
               SCALE = SCALE * SC
            END IF
         CASE (6) ! (T,N)
            CALL RECGSYL_P(PROCS, UPLOSIGN, SCALE, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, D(0,0), LDD, E(0,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
!$OMP  PARALLEL SECTIONS PRIVATE(INF) REDUCTION(MIN:INFO)
!$OMP  SECTION
            CALL RECSY_AXB('F','U','T','N', M2, M1, N1, N1, -1D0, A(0,M1), LDA, E(0,0), LDE, B(0,0), LDB, SCALE, E(M1,0), LDE, MACHINE, WORKSPACE)
            CALL RECSY_AXB('F','V','T','N', M2, M1, N1, N1, -SIGN, C(0,M1), LDC, E(0,0), LDE, D(0,0), LDD,   1D0, E(M1,0), LDE, MACHINE, WORKSPACE)
            CALL RECGSYL_P(PROCS/2, UPLOSIGN, SC, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,M1), LDC, D(0,0), LDD, E(M1,0), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
!$OMP  SECTION
            CALL RECSY_AXB('V','F','T','N', M1, M1, N1, N2, -1D0, A(0,0), LDA, E(0,0), LDE, B(0,N1), LDB, SCALE, E(0,N1), LDE, MACHINE, WORKSPACE(2*M2*N1))
            CALL RECSY_AXB('U','F','T','N', M1, M1, N1, N2, -SIGN, C(0,0), LDC, E(0,0), LDE, D(0,N1), LDD,   1D0, E(0,N1), LDE, MACHINE, WORKSPACE(2*M2*N1))
            CALL RECGSYL_P(PROCS/2, UPLOSIGN, SC2, M1, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,0), LDC, D(N1,N1), LDD, E(0,N1), LDE, INF, MACHINE, WORKSPACE(2*M2*N1), WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
!$OMP  END PARALLEL SECTIONS
            IF (SC.NE.1D0.AND.SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N2, SC, E(0,N1), LDE)
               CALL RECSY_SCAL(M2, N1, SC2, E(M1,0), LDE)
               CALL RECSY_SCAL(M1, N1, SC*SC2, E(0,0), LDE)
               SCALE = SCALE * SC * SC2
            ELSE IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N, SC, E(0,0), LDE)
               SCALE = SCALE * SC
            ELSE IF (SC2.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC2, E(0,0), LDE)
               SCALE = SCALE * SC2
            END IF
            IF (M.GT.N) THEN
               CALL RECSY_AXB('F','U','T','N', M2, M1, N, N2, -1D0, A(0,M1), LDA, E(0,0), LDE, B(0,N1), LDB, SCALE, E(M1,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','V','T','N', M2, M1, N, N2, -SIGN, C(0,M1), LDC, E(0,0), LDE, D(0,N1), LDD,   1D0, E(M1,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('V','F','T','N', M2, M2, N1, N2, -1D0, A(M1,M1), LDA, E(M1,0), LDE, B(0,N1), LDB, 1D0, E(M1,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('U','F','T','N', M2, M2, N1, N2, -SIGN, C(M1,M1), LDC, E(M1,0), LDE, D(0,N1), LDD, 1D0, E(M1,N1), LDE, MACHINE, WORKSPACE)
            ELSE
               CALL RECSY_AXB('V','F','T','N', M2, M, N1, N2, -1D0, A(0,M1), LDA, E(0,0), LDE, B(0,N1), LDB, SCALE, E(M1,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('U','F','T','N', M2, M, N1, N2, -SIGN, C(0,M1), LDC, E(0,0), LDE, D(0,N1), LDD,   1D0, E(M1,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','U','T','N', M2, M1, N2, N2, -1D0, A(0,M1), LDA, E(0,N1), LDE, B(N1,N1), LDB, 1D0, E(M1,N1), LDE, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','V','T','N', M2, M1, N2, N2, -SIGN, C(0,M1), LDC, E(0,N1), LDE, D(N1,N1), LDD, 1D0, E(M1,N1), LDE, MACHINE, WORKSPACE)
            END IF
            CALL RECGSYL_P(PROCS, UPLOSIGN, SC, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,M1), LDC, D(N1,N1), LDD, E(M1,N1), LDE, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, E(0,0), LDE)
               CALL RECSY_SCAL(M1, N2, SC, E(0,N1), LDE)
               SCALE = SCALE * SC
            END IF
         END SELECT
      END IF
      IF (INFO.GT.N) INFO = 0
      RETURN
      END SUBROUTINE RECGSYL_P

