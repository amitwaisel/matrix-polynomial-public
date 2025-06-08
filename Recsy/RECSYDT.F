!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      RECURSIVE SUBROUTINE RECSYDT(UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE, WORKSPACE, WKSIZE )
!     Solves op(A)*X*op(B) + sign*X = scale*C
!     A (M*M) and B (N*N) quasi-triangular
!     
!     Iff close eigenvalues are detected, INFO is non-zero and elements 
!     of the C (M*N) matrix are no longer valid
!
      IMPLICIT NONE
!
      INTENT(IN) M,N,A,LDA,B,LDB,LDC,UPLOSIGN,MACHINE
      INTENT(OUT) INFO,SCALE
      INTENT(INOUT) C,WORKSPACE,WKSIZE

!     .. Scalar Arguments ..
      INTEGER          :: M, N, LDA, LDB, LDC, INFO, UPLOSIGN
      INTEGER, OPTIONAL :: WKSIZE
      DOUBLE PRECISION :: SCALE
!     .. Array Arguments ..
      DOUBLE PRECISION :: A(0:LDA-1,0:M-1), B(0:LDB-1,0:N-1), C(0:LDC-1,0:N-1)
      DOUBLE PRECISION, OPTIONAL :: WORKSPACE(0:*), MACHINE(0:*)
!     .. Local Scalars ..
      INTEGER          :: I, J, K, L, M1, M2, N1, N2, INF, WORKNEEDED, WKS
      DOUBLE PRECISION :: SIGN, SC, TMPA, TMPB, MACHINETEMP(0:9)
      CHARACTER*1      :: TA, TB
!     .. Local Arrays ..
      DOUBLE PRECISION :: Z(0:31,0:31)
      DOUBLE PRECISION :: VEC(0:15)
      DOUBLE PRECISION, ALLOCATABLE :: WORKTEMP(:)
!     .. Intrinsic Functions ..
      INTRINSIC           MOD
!     .. Explicit Declarations
!     .. Executable Statements ..
!
      IF (.NOT.PRESENT(WKSIZE)) THEN
         WKS = 0
      ELSE
         WKS = WKSIZE
      END IF
      WORKNEEDED = M*N
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
            CALL RECSYDT(UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINETEMP, WORKTEMP, WORKNEEDED )
            DEALLOCATE(WORKTEMP)
         ELSE
            CALL RECSYDT(UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINETEMP, WORKSPACE, WKS )
         END IF
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
            CALL RECSYDT_TR0(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
         CASE (1)
            CALL RECSYDT_TR1(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
         CASE (2)
            CALL RECSYDT_TR2(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
         CASE (3)
            CALL RECSYDT_TR3(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
         CASE (4)
            CALL RECSYDT_TR4(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
         CASE (5)
            CALL RECSYDT_TR5(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
         CASE (6)
            CALL RECSYDT_TR6(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
         CASE (7)
            CALL RECSYDT_TR7(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
         END SELECT
         IF (INFO.NE.0) THEN
            DO L=0,N-1
               DO K=0,N-1
                  IF (MOD(UPLOSIGN/2,2).EQ.0) THEN
                     TMPB = B(L,K)
                  ELSE
                     TMPB = B(K,L)
                  END IF
                  DO J=0,M-1
                     DO I=0,M-1
                        IF (UPLOSIGN.LT.4) THEN
                           TMPA = A(I,J)
                        ELSE
                           TMPA = A(J,I)
                        END IF
                        Z(K*M+I,L*M+J) = TMPB*TMPA
                     END DO
                  END DO
               END DO
            END DO
            DO I=0,M*N-1
               Z(I,I) = Z(I,I) + SIGN
            END DO
            CALL RECSY_SOLVE_CP(M*N,32,VEC,Z,SC, INFO, MACHINE)
            SCALE = SCALE * SC
            DO J=0,N-1
               DO I=0,M-1
                  C(I,J) = VEC(I+J*M)
               END DO
            END DO
         END IF
         RETURN
      ENDIF

      IF (UPLOSIGN/4.EQ.0) THEN
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
!      IF (M.GE.N.OR.N.LE.4) THEN ! Split M dimension
         M1 = M/2
         M2 = M - M1
            IF (A(M1,M1-1).NE.0D0) THEN
            M1 = M1+1
            M2 = M2-1
         END IF
         SELECT CASE (UPLOSIGN/2)
         CASE (0,1) ! (N,N) and (N,T)

! A*X*B+X = C
! [ A11 A12 ] [ X1 ] B + [ X1 ] = [ C1 ]
! [  0  A22 ] [ X2 ]     [ X2 ]   [ C2 ]
!
! A11*X1*B + A12*X2*B +/- X1 = C1
!            A22*X2*B +/- X2 = C2

            CALL RECSYDT(UPLOSIGN, SCALE, M2, N, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_AXB('F','V','N',TB, M1, M2, N, N, -1D0, A(0,M1), LDA, C(M1,0), LDC, B, LDB, SCALE, C(0,0), LDC, MACHINE, WORKSPACE)
            CALL RECSYDT(UPLOSIGN, SC, M1, N, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N, SC, C(M1,0), LDC)
               SCALE = SCALE * SC
            END IF
         CASE (2,3) ! (T,N) and (T,T)

! A'*X*B+X = C
! [ A11' 0    ] [ X1 ] B + [ X1 ] = [ C1 ]
! [ A12' A22' ] [ X2 ]   + [ X2 ] = [ C2 ]
!
! A11'*X1*B + X1 = C1
! A12'*X1*B + A22'*X2*B + X2 = C2

            CALL RECSYDT(UPLOSIGN, SCALE, M1, N, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_AXB('F','V','T',TB, M2, M1, N, N, -1D0, A(0,M1), LDA, C(0,0), LDC, B, LDB, SCALE, C(M1,0), LDC, MACHINE, WORKSPACE)
            CALL RECSYDT(UPLOSIGN, SC, M2, N, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N, SC, C(0,0), LDC)
               SCALE = SCALE * SC
            END IF
         END SELECT
      ELSE IF (N/2.GE.M.OR.M.LE.4) THEN ! Split N dimension
!      ELSE IF (.TRUE.) THEN ! Split N dimension
         N1 = N/2
         N2 = N - N1
            IF (B(N1,N1-1).NE.0D0) THEN
            N1 = N1+1
            N2 = N2-1
         END IF
         SELECT CASE (UPLOSIGN/2)
         CASE (0,2) ! (N,N) and (T,N)
! A*X*B+X = C
! A [ X1 X2 ] [ B11 B12 ] + [ X1 X2 ] = [ C1 C2 ]
!             [ 0   B22 ] 
! 
! A*X1*B11 + X1 = C1
! A*X1*B12 + A*X2*B22 + X2 = C2
!
            CALL RECSYDT(UPLOSIGN, SCALE, M, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_AXB('V','F',TA,'N', M, M, N1, N2, -1D0, A(0,0), LDA, C(0,0), LDC, B(0,N1), LDB, SCALE, C(0,N1), LDC, MACHINE, WORKSPACE)
            CALL RECSYDT(UPLOSIGN, SC, M, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,N1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, C(0,0), LDC)
               SCALE = SCALE * SC
            END IF
         CASE (1,3) ! (N,T) and (T,T)
! A*X*B'+X = C
! A [ X1 X2 ] [ B11' 0    ] + [ X1 X2 ] = [ C1 C2 ]
!             [ B12' B22' ] 
!
! A*X1*B11' + A*X2*B12' + X1 = C1
! A*X2*B22' + X2 = C2
!
            CALL RECSYDT(UPLOSIGN, SCALE, M, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,N1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            CALL RECSY_AXB('V','F',TA,'T', M, M, N2, N1, -1D0, A(0,0), LDA, C(0,N1), LDC, B(0,N1), LDB, SCALE, C(0,0), LDC, MACHINE, WORKSPACE)
            CALL RECSYDT(UPLOSIGN, SC, M, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
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
! [ A11 A12 ] [ X11 X12 ] [ B11 B12 ] + [ X11 X12 ] = [ C11 C12 ]
! [  0  A22 ] [ X21 X22 ] [  0  B22 ]   [ X21 X22 ]   [ C21 C22 ]
!
! [ A11*X11 + X12*X21   A11*X12 + A12*X22 ]
! [           A22*X21             A22*X22 ]
!
! A11*X11*B11 + A12*X21*B11   + X11 = C11
! A11*X11*B12 + A12*X21*B12 + A11*X12*B22 + A12*X22*B22 + X12 = C12
! A22*X21*B11 + X21 = C21
! A22*X21*B12 + A22*X22*B22 + X22 = C22
!
            CALL RECSYDT(UPLOSIGN, SCALE, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_AXB('F','V','N','N', M1, M2, N1, N1, -1D0, A(0,M1), LDA, C(M1,0), LDC, B(0,0), LDB, SCALE, C(0,0), LDC, MACHINE, WORKSPACE)
            CALL RECSYDT(UPLOSIGN, SC, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N1, SC, C(M1,0), LDC)
               SCALE = SCALE * SC
            END IF
            CALL RECSY_AXB('V','F','N','N', M2, M2, N1, N2, -1D0, A(M1,M1), LDA, C(M1,0), LDC, B(0,N1), LDB, SCALE, C(M1,N1), LDC, MACHINE, WORKSPACE)
            CALL RECSYDT(UPLOSIGN, SC, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,N1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, C(0,0), LDC)
               SCALE = SCALE * SC
            END IF
            IF (.FALSE.) THEN
               CALL RECSY_AXB('V','F','N','N', M1, M1, N1, N2, -1D0, A(0,0), LDA, C(0,0), LDC, B(0,N1), LDB, SCALE, C(0,N1), LDC, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','F','N','N', M1, M2, N1, N2, -1D0, A(0,M1), LDA, C(M1,0), LDC, B(0,N1), LDB, 1D0, C(0,N1), LDC, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','V','N','N', M1, M2, N2, N2, -1D0, A(0,M1), LDA, C(M1,N1), LDC, B(N1,N1), LDB, 1D0, C(0,N1), LDC, MACHINE, WORKSPACE)
            ELSE IF (M.GT.N) THEN
               CALL RECSY_AXB('F','V','N','N', M1, M2,  N, N2, -1D0, A(0,M1), LDA, C(M1,0), LDC, B(0,N1), LDB, SCALE, C(0,N1), LDC, MACHINE, WORKSPACE)
               CALL RECSY_AXB('V','F','N','N', M1, M1, N1, N2, -1D0, A(0,0), LDA, C(0,0), LDC, B(0,N1), LDB, 1D0, C(0,N1), LDC, MACHINE, WORKSPACE)
            ELSE
               CALL RECSY_AXB('V','F','N','N', M1,  M, N1, N2, -1D0, A(0,0), LDA, C(0,0), LDC, B(0,N1), LDB, SCALE, C(0,N1), LDC, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','V','N','N', M1, M2, N2, N2, -1D0, A(0,M1), LDA, C(M1,N1), LDC, B(N1,N1), LDB, 1D0, C(0,N1), LDC, MACHINE, WORKSPACE)
            END IF

            CALL RECSYDT(UPLOSIGN, SC, M1, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,N1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, C(0,0), LDC)
               CALL RECSY_SCAL(M2, N2, SC, C(M1,N1), LDC)
               SCALE = SCALE * SC
            END IF
         CASE (1) ! (N,T)
! [ A11 A12 ] [ X11 X12 ] [ B11'  0   ] + [ X11 X12 ] = [ C11 C12 ]
! [  0  A22 ] [ X21 X22 ] [ B12' B22' ]   [ X21 X22 ]   [ C21 C22 ]
!
! [ A11*X11 + X12*X21   A11*X12 + A12*X22 ]
! [           A22*X21             A22*X22 ]
!
! A11*X11*B11' + A12*X21*B11' + A11*X12*B12' + A12*X22*B12' + X11 = C11
! A11*X12*B22' + A12*X22*B22' + X12 = C12
! A22*X21*B11' + A22*X22*B12' + X21 = C21
! A22*X22*B22' + X22 = C22
!
            CALL RECSYDT(UPLOSIGN, SCALE, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,N1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            CALL RECSY_AXB('F','V','N','T', M1, M2, N2, N2, -1D0, A(0,M1), LDA, C(M1,N1), LDC, B(N1,N1), LDB, SCALE, C(0,N1), LDC, MACHINE, WORKSPACE)
            CALL RECSYDT(UPLOSIGN, SC, M1, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,N1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M2, N2, SC, C(M1,N1), LDC)   
               SCALE = SCALE * SC
            END IF
            CALL RECSY_AXB('V','F','N','T', M2, M2, N2, N1, -1D0, A(M1,M1), LDA, C(M1,N1), LDC, B(0,N1), LDB, SCALE, C(M1,0), LDC, MACHINE, WORKSPACE)
            CALL RECSYDT(UPLOSIGN, SC, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC, C(0,N1), LDC)
               SCALE = SCALE * SC
            END IF
            IF (.FALSE.) THEN
               CALL RECSY_AXB('F','V','N','T', M1, M2, N1, N1, -1D0, A(0,M1), LDA, C(M1,0), LDC, B(0,0), LDB, SCALE, C(0,0), LDC, MACHINE, WORKSPACE)
               CALL RECSY_AXB('V','F','N','T', M1, M1, N2, N1, -1D0, A(0,0), LDA, C(0,N1), LDC, B(0,N1), LDB, 1D0, C(0,0), LDC, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','F','N','T', M1, M2, N2, N1, -1D0, A(0,M1), LDA, C(M1,N1), LDC, B(0,N1), LDB, 1D0, C(0,0), LDC, MACHINE, WORKSPACE)
            ELSE IF (M.GT.N) THEN
               CALL RECSY_AXB('F','V','N','T', M1, M2, N, N1, -1D0, A(0,M1), LDA, C(M1,0), LDC, B(0,0), LDB, SCALE, C(0,0), LDC, MACHINE, WORKSPACE)
               CALL RECSY_AXB('V','F','N','T', M1, M1, N2, N1, -1D0, A(0,0), LDA, C(0,N1), LDC, B(0,N1), LDB, 1D0, C(0,0), LDC, MACHINE, WORKSPACE)
            ELSE
               CALL RECSY_AXB('V','F','N','T', M1, M, N2, N1, -1D0, A(0,0), LDA, C(0,N1), LDC, B(0,N1), LDB, SCALE, C(0,0), LDC, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','V','N','T', M1, M2, N1, N1, -1D0, A(0,M1), LDA, C(M1,0), LDC, B(0,0), LDB, 1D0, C(0,0), LDC, MACHINE, WORKSPACE)
            END IF
            CALL RECSYDT(UPLOSIGN, SC, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC, C(0,N1), LDC)
               CALL RECSY_SCAL(M2, N1, SC, C(M1,0), LDC)
               SCALE = SCALE * SC
            END IF
         CASE (2) ! (T,N)
! [ A11'   0  ] [ X11 X12 ] [ B11 B12 ] + [ X11 X12 ] = [ C11 C12 ]
! [ A12' A22' ] [ X21 X22 ] [  0  B22 ]   [ X21 X22 ]   [ C21 C22 ]
!
! [            A11'*X11             A11'*X12 ]
! [ A12'*X11 + A22'*X21  A12'*X12 + A22'*X22 ]
!
! A11'*X11*B11    + X11 = C11
! A11'*X11*B12 + A11'*X12*B22 + X12 = C12
! A12'*X11*B11 + A22'*X21*B11 + X21 = C21
! A12'*X11*B12 + A22'*X21*B12 + A12'*X12*B22 + A22'*X22*B22 + X22 = C22
!
            CALL RECSYDT(UPLOSIGN, SCALE, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            CALL RECSY_AXB('F','V','T','N', M2, M1, N1, N1, -1D0, A(0,M1), LDA, C(0,0), LDC, B(0,0), LDB, SCALE, C(M1,0), LDC, MACHINE, WORKSPACE)
            CALL RECSYDT(UPLOSIGN, SC, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N1, SC, C(0,0), LDC)
               SCALE = SCALE * SC
            END IF
            CALL RECSY_AXB('V','F','T','N', M1, M1, N1, N2, -1D0, A(0,0), LDA, C(0,0), LDC, B(0,N1), LDB, SCALE, C(0,N1), LDC, MACHINE, WORKSPACE)
            CALL RECSYDT(UPLOSIGN, SC, M1, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,N1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, C(0,0), LDC)
               SCALE = SCALE * SC
            END IF
            IF (.FALSE.) THEN
               CALL RECSY_AXB('F','F','T','N', M2, M1, N1, N2, -1D0, A(0,M1), LDA, C(0,0), LDC, B(0,N1), LDB, SCALE, C(M1,N1), LDC, MACHINE, WORKSPACE)
               CALL RECSY_AXB('V','F','T','N', M2, M2, N1, N2, -1D0, A(M1,M1), LDA, C(M1,0), LDC, B(0,N1), LDB, 1D0, C(M1,N1), LDC, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','V','T','N', M2, M1, N2, N2, -1D0, A(0,M1), LDA, C(0,N1), LDC, B(N1,N1), LDB, 1D0, C(M1,N1), LDC, MACHINE, WORKSPACE)
            ELSE IF (M.GT.N) THEN
               CALL RECSY_AXB('F','V','T','N', M2, M1, N, N2, -1D0, A(0,M1), LDA, C(0,0), LDC, B(0,N1), LDB, SCALE, C(M1,N1), LDC, MACHINE, WORKSPACE)
               CALL RECSY_AXB('V','F','T','N', M2, M2, N1, N2, -1D0, A(M1,M1), LDA, C(M1,0), LDC, B(0,N1), LDB, 1D0, C(M1,N1), LDC, MACHINE, WORKSPACE)
            ELSE
               CALL RECSY_AXB('V','F','T','N', M2, M, N1, N2, -1D0, A(0,M1), LDA, C(0,0), LDC, B(0,N1), LDB, SCALE, C(M1,N1), LDC, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','V','T','N', M2, M1, N2, N2, -1D0, A(0,M1), LDA, C(0,N1), LDC, B(N1,N1), LDB, 1D0, C(M1,N1), LDC, MACHINE, WORKSPACE)
            END IF
            CALL RECSYDT(UPLOSIGN, SC, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,N1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N1, SC, C(0,0), LDC)
               CALL RECSY_SCAL(M1, N2, SC, C(0,N1), LDC)
               SCALE = SCALE * SC
            END IF
         CASE (3) ! (T,T)
! [ A11'   0  ] [ X11 X12 ] [ B11'   0  ] + [ X11 X12 ] = [ C11 C12 ]
! [ A12' A22' ] [ X21 X22 ] [ B12' B22' ]   [ X21 X22 ]   [ C21 C22 ]
!
! [            A11'*X11             A11'*X12 ]
! [ A12'*X11 + A22'*X21  A12'*X12 + A22'*X22 ]
!
! A11'*X11*B11' + A11'*X12*B12' + X11 = C11
! A11'*X12*B22' + X12 = C12
! A12'*X11*B11' + A22'*X21*B11' + A12'*X12*B12' + A22'*X22*B12' + X21 = C21
! A12'*X12*B22' + A22'*X22*B22' + X22 = C22
!
            CALL RECSYDT(UPLOSIGN, SCALE, M1, N2, A(0,0), LDA, B(N1,N1), LDB, C(0,N1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            CALL RECSY_AXB('F','V','T','T', M2, M1, N2, N2, -1D0, A(0,M1), LDA, C(0,N1), LDC, B(N1,N1), LDB, SCALE, C(M1,N1), LDC, MACHINE, WORKSPACE)
            CALL RECSYDT(UPLOSIGN, SC, M2, N2, A(M1,M1), LDA, B(N1,N1), LDB, C(M1,N1), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF+N1.LT.INFO) INFO = INF+N1
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M1, N2, SC, C(0,N1), LDC)
               SCALE = SCALE * SC
            END IF
            CALL RECSY_AXB('V','F','T','T', M1, M1, N2, N1, -1D0, A(0,0), LDA, C(0,N1), LDC, B(0,N1), LDB, SCALE, C(0,0), LDC, MACHINE, WORKSPACE)
            CALL RECSYDT(UPLOSIGN, SC, M1, N1, A(0,0), LDA, B(0,0), LDB, C(0,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
            IF (INF.NE.0.AND.INF.LT.INFO) INFO = INF
            IF (SC.NE.1D0) THEN
               CALL RECSY_SCAL(M, N2, SC, C(0,N1), LDC)
               SCALE = SCALE * SC
            END IF
            IF (.FALSE.) THEN
               CALL RECSY_AXB('F','V','T','T', M2, M1, N1, N1, -1D0, A(0,M1), LDA, C(0,0), LDC, B(0,0), LDB, SCALE, C(M1,0), LDC, MACHINE, WORKSPACE)
               CALL RECSY_AXB('F','F','T','T', M2, M1, N2, N1, -1D0, A(0,M1), LDA, C(0,N1), LDC, B(0,N1), LDB, 1D0, C(M1,0), LDC, MACHINE, WORKSPACE)
               CALL RECSY_AXB('V','F','T','T', M2, M2, N2, N1, -1D0, A(M1,M1), LDA, C(M1,N1), LDC, B(0,N1), LDB, 1D0, C(M1,0), LDC, MACHINE, WORKSPACE)
            ELSE IF (M.GT.N) THEN
               CALL RECSY_AXB('F','V','T','T', M2, M1, N, N1, -1D0, A(0,M1), LDA, C(0,0), LDC, B(0,0), LDB, SCALE, C(M1,0), LDC, MACHINE, WORKSPACE)
               CALL RECSY_AXB('V','F','T','T', M2, M2, N2, N1, -1D0, A(M1,M1), LDA, C(M1,N1), LDC, B(0,N1), LDB, 1D0, C(M1,0), LDC, MACHINE, WORKSPACE)
            ELSE
               CALL RECSY_AXB('F','V','T','T', M2, M1, N1, N1, -1D0, A(0,M1), LDA, C(0,0), LDC, B(0,0), LDB, SCALE, C(M1,0), LDC, MACHINE, WORKSPACE)
               CALL RECSY_AXB('V','F','T','T', M2, M, N2, N1, -1D0, A(0,M1), LDA, C(0,N1), LDC, B(0,N1), LDB, 1D0, C(M1,0), LDC, MACHINE, WORKSPACE)
            END IF
            CALL RECSYDT(UPLOSIGN, SC, M2, N1, A(M1,M1), LDA, B(0,0), LDB, C(M1,0), LDC, INF, MACHINE, WORKSPACE, WKSIZE)
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
      END SUBROUTINE RECSYDT

