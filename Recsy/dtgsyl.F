!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!
      SUBROUTINE DTGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, &
          D, LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK, &
              IWORK, INFO )
!
!  -- LAPACK routine (MODIFIED version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, &
                         LWORK, M, N
      DOUBLE PRECISION   DIF, SCALE
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ), &
                         D( LDD, * ), E( LDE, * ), F( LDF, * ), &
                         WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DTGSYL solves the generalized Sylvester equation:
!
!              A * R - L * B = scale * C                 (1)
!              D * R - L * E = scale * F
!
!  where R and L are unknown m-by-n matrices, (A, D), (B, E) and
!  (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,
!  respectively, with real entries. (A, D) and (B, E) must be in
!  generalized (real) Schur canonical form, i.e. A, B are upper quasi
!  triangular and D, E are upper triangular.
!
!  The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output
!  scaling factor chosen to avoid overflow.
!
!  In matrix notation (1) is equivalent to solve  Zx = scale b, where
!  Z is defined as
!
!             Z = [ kron(In, A)  -kron(B', Im) ]         (2)
!                 [ kron(In, D)  -kron(E', Im) ].
!
!  Here Ik is the identity matrix of size k and X' is the transpose of
!  X. kron(X, Y) is the Kronecker product between the matrices X and Y.
!
!  If TRANS = 'T', DTGSYL solves the transposed system Z'*y = scale*b,
!  which is equivalent to solve for R and L in
!
!              A' * R  + D' * L   = scale *  C           (3)
!              R  * B' + L  * E'  = scale * (-F)
!
!  This case (TRANS = 'T') is used to compute an one-norm-based estimate
!  of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)
!  and (B,E), using DLACON.
!
!  If IJOB >= 1, DTGSYL computes a Frobenius norm-based estimate
!  of Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the
!  reciprocal of the smallest singular value of Z. See [1-2] for more
!  information.
!
!  This is a level 3 BLAS algorithm.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          = 'N', solve the generalized Sylvester equation (1).
!          = 'T', solve the 'transposed' system (3).
!
!  IJOB    (input) INTEGER
!          Specifies what kind of functionality to be performed.
!           =0: solve (1) only.
!           =1: The functionality of 0 and 3.
!           =2: The functionality of 0 and 4.
!           =3: Only an estimate of Dif[(A,D), (B,E)] is computed.
!               (look ahead strategy IJOB  = 1 is used).
!           =4: Only an estimate of Dif[(A,D), (B,E)] is computed.
!               ( DGECON on sub-systems is used ).
!          Not referenced if TRANS = 'T'.
!
!  M       (input) INTEGER
!          The order of the matrices A and D, and the row dimension of
!          the matrices C, F, R and L.
!
!  N       (input) INTEGER
!          The order of the matrices B and E, and the column dimension
!          of the matrices C, F, R and L.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA, M)
!          The upper quasi triangular matrix A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1, M).
!
!  B       (input) DOUBLE PRECISION array, dimension (LDB, N)
!          The upper quasi triangular matrix B.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B. LDB >= max(1, N).
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC, N)
!          On entry, C contains the right-hand-side of the first matrix
!          equation in (1) or (3).
!          On exit, if IJOB = 0, 1 or 2, C has been overwritten by
!          the solution R. If IJOB = 3 or 4 and TRANS = 'N', C holds R,
!          the solution achieved during the computation of the
!          Dif-estimate.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1, M).
!
!  D       (input) DOUBLE PRECISION array, dimension (LDD, M)
!          The upper triangular matrix D.
!
!  LDD     (input) INTEGER
!          The leading dimension of the array D. LDD >= max(1, M).
!
!  E       (input) DOUBLE PRECISION array, dimension (LDE, N)
!          The upper triangular matrix E.
!
!  LDE     (input) INTEGER
!          The leading dimension of the array E. LDE >= max(1, N).
!
!  F       (input/output) DOUBLE PRECISION array, dimension (LDF, N)
!          On entry, F contains the right-hand-side of the second matrix
!          equation in (1) or (3).
!          On exit, if IJOB = 0, 1 or 2, F has been overwritten by
!          the solution L. If IJOB = 3 or 4 and TRANS = 'N', F holds L,
!          the solution achieved during the computation of the
!          Dif-estimate.
!
!  LDF     (input) INTEGER
!          The leading dimension of the array F. LDF >= max(1, M).
!
!  DIF     (output) DOUBLE PRECISION
!          On exit DIF is the reciprocal of a lower bound of the
!          reciprocal of the Dif-function, i.e. DIF is an upper bound of
!          Dif[(A,D), (B,E)] = sigma_min(Z), where Z as in (2).
!          IF IJOB = 0 or TRANS = 'T', DIF is not touched.
!
!  SCALE   (output) DOUBLE PRECISION
!          On exit SCALE is the scaling factor in (1) or (3).
!          If 0 < SCALE < 1, C and F hold the solutions R and L, resp.,
!          to a slightly perturbed system but the input matrices A, B, D
!          and E have not been changed. If SCALE = 0, C and F hold the
!          solutions R and L, respectively, to the homogeneous system
!          with C = F = 0. Normally, SCALE = 1.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          If IJOB = 0, WORK is not referenced.  Otherwise,
!          on exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK > = 1.
!          If IJOB = 1 or 2 and TRANS = 'N', LWORK >= 2*M*N.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  IWORK   (workspace) INTEGER array, dimension (M+N+6)
!
!  INFO    (output) INTEGER
!            =0: successful exit
!            <0: If INFO = -i, the i-th argument had an illegal value.
!            >0: (A, D) and (B, E) have common or close eigenvalues.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
!     Umea University, S-901 87 Umea, Sweden.
!
!  [1] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software
!      for Solving the Generalized Sylvester Equation and Estimating the
!      Separation between Regular Matrix Pairs, Report UMINF - 93.23,
!      Department of Computing Science, Umea University, S-901 87 Umea,
!      Sweden, December 1993, Revised April 1994, Also as LAPACK Working
!      Note 75.  To appear in ACM Trans. on Math. Software, Vol 22,
!      No 1, 1996.
!
!  [2] B. Kagstrom, A Perturbation Analysis of the Generalized Sylvester
!      Equation (AR - LB, DR - LE ) = (C, F), SIAM J. Matrix Anal.
!      Appl., 15(4):1045-1060, 1994
!
!  [3] B. Kagstrom and L. Westin, Generalized Schur Methods with
!      Condition Estimators for Solving the Generalized Sylvester
!      Equation, IEEE Transactions on Automatic Control, Vol. 34, No. 7,
!      July 1989, pp 745-751.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, NOTRAN
      INTEGER            I, IE, IFUNC, IROUND, IS, ISOLVE, J, JE, JS, K, &
                         LINFO, LWMIN, MB, NB, P, PPQQ, PQ, Q
      DOUBLE PRECISION   DSCALE, DSUM, SCALE2, SCALOC
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLACPY, DSCAL, DTGSY2, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, SQRT
!     .. Explicit Declarations
      INCLUDE 'RECSY_DECL.F'
!     ..
!     .. Executable Statements ..
!
!     Decode and test input parameters
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!
      IF( ( IJOB.EQ.1 .OR. IJOB.EQ.2 ) .AND. NOTRAN ) THEN
         LWMIN = MAX( 1, 2*M*N )
      ELSE
         LWMIN = 1
      END IF
!
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -1
      ELSE IF( ( IJOB.LT.0 ) .OR. ( IJOB.GT.4 ) ) THEN
         INFO = -2
      ELSE IF( M.LE.0 ) THEN
         INFO = -3
      ELSE IF( N.LE.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LDD.LT.MAX( 1, M ) ) THEN
         INFO = -12
      ELSE IF( LDE.LT.MAX( 1, N ) ) THEN
         INFO = -14
      ELSE IF( LDF.LT.MAX( 1, M ) ) THEN
         INFO = -16
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -20
      END IF
!
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTGSYL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Determine optimal block sizes MB and NB
!
      MB = ILAENV( 2, 'DTGSYL', TRANS, M, N, -1, -1 )
      NB = ILAENV( 5, 'DTGSYL', TRANS, M, N, -1, -1 )
!
      ISOLVE = 1
      IFUNC = 0
      IF( IJOB.GE.3 .AND. NOTRAN ) THEN
         IFUNC = IJOB - 2
         DO 10 J = 1, N
            CALL DCOPY( M, ZERO, 0, C( 1, J ), 1 )
            CALL DCOPY( M, ZERO, 0, F( 1, J ), 1 )
   10    CONTINUE
      ELSE IF( IJOB.GE.1 .AND. NOTRAN ) THEN
         ISOLVE = 2
      END IF
!
      IF( ( MB.LE.1 .AND. NB.LE.1 ) .OR. ( MB.GE.M .AND. NB.GE.N ) ) &
           THEN
!
         DO 30 IROUND = 1, ISOLVE
!
!           Use unblocked Level 2 solver
!
            DSCALE = ZERO
            DSUM = ONE
            PQ = 0
            CALL DTGSY2( TRANS, IFUNC, M, N, A, LDA, B, LDB, C, LDC, D, &
                         LDD, E, LDE, F, LDF, SCALE, DSUM, DSCALE, &
                         IWORK, PQ, INFO )
            IF( DSCALE.NE.ZERO ) THEN
               IF( IJOB.EQ.1 .OR. IJOB.EQ.3 ) THEN
                  DIF = SQRT( DBLE( 2*M*N ) ) / ( DSCALE*SQRT( DSUM ) )
               ELSE
                  DIF = SQRT( DBLE( PQ ) ) / ( DSCALE*SQRT( DSUM ) )
               END IF
            END IF
!
            IF( ISOLVE.EQ.2 .AND. IROUND.EQ.1 ) THEN
               IFUNC = IJOB
               SCALE2 = SCALE
               CALL DLACPY( 'F', M, N, C, LDC, WORK, M )
               CALL DLACPY( 'F', M, N, F, LDF, WORK( M*N+1 ), M )
               DO 20 J = 1, N
                  CALL DCOPY( M, ZERO, 0, C( 1, J ), 1 )
                  CALL DCOPY( M, ZERO, 0, F( 1, J ), 1 )
   20          CONTINUE
            ELSE IF( ISOLVE.EQ.2 .AND. IROUND.EQ.2 ) THEN
               CALL DLACPY( 'F', M, N, WORK, M, C, LDC )
               CALL DLACPY( 'F', M, N, WORK( M*N+1 ), M, F, LDF )
               SCALE = SCALE2
            END IF
   30    CONTINUE
!
         RETURN
      END IF
!
!     Determine block structure of A
!
      P = 0
      I = 1
   40 CONTINUE
      IF( I.GT.M ) &
         GO TO 50
      P = P + 1
      IWORK( P ) = I
      I = I + MB
      IF( I.GE.M ) &
         GO TO 50
      IF( A( I, I-1 ).NE.ZERO ) &
         I = I + 1
      GO TO 40
   50 CONTINUE
!
      IWORK( P+1 ) = M + 1
      IF( IWORK( P ).EQ.IWORK( P+1 ) ) &
         P = P - 1
!
!     Determine block structure of B
!
      Q = P + 1
      J = 1
   60 CONTINUE
      IF( J.GT.N ) &
         GO TO 70
      Q = Q + 1
      IWORK( Q ) = J
      J = J + NB
      IF( J.GE.N ) &
         GO TO 70
      IF( B( J, J-1 ).NE.ZERO ) &
         J = J + 1
      GO TO 60
   70 CONTINUE
!
      IWORK( Q+1 ) = N + 1
      IF( IWORK( Q ).EQ.IWORK( Q+1 ) ) &
         Q = Q - 1
!
      IF( NOTRAN ) THEN
!
         IF (IJOB.EQ.0) THEN
            CALL RECGCSY(1, SCALE, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, INFO)
            RETURN
         ELSE
!
            DO 150 IROUND = 1, ISOLVE
!
!              Solve (I, J)-subsystem
!                  A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
!                  D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
!              for I = P, P - 1,..., 1; J = 1, 2,..., Q
!
               DSCALE = ZERO
               DSUM = ONE
               PQ = 0
               SCALE = ONE
               DO 130 J = P + 2, Q
                  JS = IWORK( J )
                  JE = IWORK( J+1 ) - 1
                  NB = JE - JS + 1
                  DO 120 I = P, 1, -1
                     IS = IWORK( I )
                     IE = IWORK( I+1 ) - 1
                     MB = IE - IS + 1
                     PPQQ = 0
                     CALL DTGSY2( TRANS, IFUNC, MB, NB, A( IS, IS ), LDA, &
                                  B( JS, JS ), LDB, C( IS, JS ), LDC, &
                                  D( IS, IS ), LDD, E( JS, JS ), LDE, &
                                  F( IS, JS ), LDF, SCALOC, DSUM, DSCALE, &
                                  IWORK( Q+2 ), PPQQ, LINFO )
                     IF( LINFO.GT.0 ) &
                        INFO = LINFO
!
                     PQ = PQ + PPQQ
                     IF( SCALOC.NE.ONE ) THEN
                        DO 80 K = 1, JS - 1
                           CALL DSCAL( M, SCALOC, C( 1, K ), 1 )
                           CALL DSCAL( M, SCALOC, F( 1, K ), 1 )
   80                   CONTINUE
                        DO 90 K = JS, JE
                           CALL DSCAL( IS-1, SCALOC, C( 1, K ), 1 )
                           CALL DSCAL( IS-1, SCALOC, F( 1, K ), 1 )
   90                   CONTINUE
                        DO 100 K = JS, JE
                           CALL DSCAL( M-IE, SCALOC, C( IE+1, K ), 1 )
                           CALL DSCAL( M-IE, SCALOC, F( IE+1, K ), 1 )
  100                   CONTINUE
                        DO 110 K = JE + 1, N
                           CALL DSCAL( M, SCALOC, C( 1, K ), 1 )
                           CALL DSCAL( M, SCALOC, F( 1, K ), 1 )
  110                   CONTINUE
                        SCALE = SCALE*SCALOC
                     END IF
!
!                    Substitute R(I, J) and L(I, J) into remaining
!                    equation.
!
                     IF( I.GT.1 ) THEN
                        CALL DGEMM( 'N', 'N', IS-1, NB, MB, -ONE, &
                                    A( 1, IS ), LDA, C( IS, JS ), LDC, ONE, &
                                    C( 1, JS ), LDC )
                        CALL DGEMM( 'N', 'N', IS-1, NB, MB, -ONE, &
                                    D( 1, IS ), LDD, C( IS, JS ), LDC, ONE, &
                                    F( 1, JS ), LDF )
                     END IF
                     IF( J.LT.Q ) THEN
                        CALL DGEMM( 'N', 'N', MB, N-JE, NB, ONE, &
                                    F( IS, JS ), LDF, B( JS, JE+1 ), LDB, &
                                    ONE, C( IS, JE+1 ), LDC )
                        CALL DGEMM( 'N', 'N', MB, N-JE, NB, ONE, &
                                    F( IS, JS ), LDF, E( JS, JE+1 ), LDE, &
                                    ONE, F( IS, JE+1 ), LDF )
                     END IF
  120             CONTINUE
  130          CONTINUE
               IF( DSCALE.NE.ZERO ) THEN
                  IF( IJOB.EQ.1 .OR. IJOB.EQ.3 ) THEN
                     DIF = SQRT( DBLE( 2*M*N ) ) / ( DSCALE*SQRT( DSUM ) )
                  ELSE
                     DIF = SQRT( DBLE( PQ ) ) / ( DSCALE*SQRT( DSUM ) )
                  END IF
               END IF
               IF( ISOLVE.EQ.2 .AND. IROUND.EQ.1 ) THEN
                  IFUNC = IJOB
                  SCALE2 = SCALE
                  CALL DLACPY( 'F', M, N, C, LDC, WORK, M )
                  CALL DLACPY( 'F', M, N, F, LDF, WORK( M*N+1 ), M )
                  DO 140 J = 1, N
                     CALL DCOPY( M, ZERO, 0, C( 1, J ), 1 )
                     CALL DCOPY( M, ZERO, 0, F( 1, J ), 1 )
  140             CONTINUE
               ELSE IF( ISOLVE.EQ.2 .AND. IROUND.EQ.2 ) THEN
                  CALL DLACPY( 'F', M, N, WORK, M, C, LDC )
                  CALL DLACPY( 'F', M, N, WORK( M*N+1 ), M, F, LDF )
                  SCALE = SCALE2
               END IF
  150       CONTINUE
!
         END IF
      ELSE
!
!        Solve transposed (I, J)-subsystem
!             A(I, I)' * R(I, J)  + D(I, I)' * L(I, J)  =  C(I, J)
!             R(I, J)  * B(J, J)' + L(I, J)  * E(J, J)' = -F(I, J)
!        for I = 1,2,..., P; J = Q, Q-1,..., 1
!
         SCALE = ONE
         DO 210 I = 1, P
            IS = IWORK( I )
            IE = IWORK( I+1 ) - 1
            MB = IE - IS + 1
            DO 200 J = Q, P + 2, -1
               JS = IWORK( J )
               JE = IWORK( J+1 ) - 1
               NB = JE - JS + 1
               CALL DTGSY2( TRANS, IFUNC, MB, NB, A( IS, IS ), LDA, &
                            B( JS, JS ), LDB, C( IS, JS ), LDC, &
                            D( IS, IS ), LDD, E( JS, JS ), LDE, &
                            F( IS, JS ), LDF, SCALOC, DSUM, DSCALE, &
                            IWORK( Q+2 ), PPQQ, LINFO )
               IF( LINFO.GT.0 ) &
                  INFO = LINFO
               IF( SCALOC.NE.ONE ) THEN
                  DO 160 K = 1, JS - 1
                     CALL DSCAL( M, SCALOC, C( 1, K ), 1 )
                     CALL DSCAL( M, SCALOC, F( 1, K ), 1 )
  160             CONTINUE
                  DO 170 K = JS, JE
                     CALL DSCAL( IS-1, SCALOC, C( 1, K ), 1 )
                     CALL DSCAL( IS-1, SCALOC, F( 1, K ), 1 )
  170             CONTINUE
                  DO 180 K = JS, JE
                     CALL DSCAL( M-IE, SCALOC, C( IE+1, K ), 1 )
                     CALL DSCAL( M-IE, SCALOC, F( IE+1, K ), 1 )
  180             CONTINUE
                  DO 190 K = JE + 1, N
                     CALL DSCAL( M, SCALOC, C( 1, K ), 1 )
                     CALL DSCAL( M, SCALOC, F( 1, K ), 1 )
  190             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
!
!              Substitute R(I, J) and L(I, J) into remaining equation.
!
               IF( J.GT.P+2 ) THEN
                  CALL DGEMM( 'N', 'T', MB, JS-1, NB, ONE, C( IS, JS ), &
                              LDC, B( 1, JS ), LDB, ONE, F( IS, 1 ), &
                              LDF )
                  CALL DGEMM( 'N', 'T', MB, JS-1, NB, ONE, F( IS, JS ), &
                              LDF, E( 1, JS ), LDE, ONE, F( IS, 1 ), &
                              LDF )
               END IF
               IF( I.LT.P ) THEN
                  CALL DGEMM( 'T', 'N', M-IE, NB, MB, -ONE, &
                              A( IS, IE+1 ), LDA, C( IS, JS ), LDC, ONE, &
                              C( IE+1, JS ), LDC )
                  CALL DGEMM( 'T', 'N', M-IE, NB, MB, -ONE, &
                              D( IS, IE+1 ), LDD, F( IS, JS ), LDF, ONE, &
                              C( IE+1, JS ), LDC )
               END IF
  200       CONTINUE
  210    CONTINUE
!
      END IF
!
      WORK( 1 ) = LWMIN
!
      RETURN
!
!     End of DTGSYL
!
      END
