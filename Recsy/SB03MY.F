!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!
      SUBROUTINE SB03MY( TRANA, N, A, LDA, C, LDC, SCALE, INFO )
!
!     MODIFIED RELEASE 4.0, WGS COPYRIGHT 1999.
!
!     PURPOSE
!
!     To solve the real Lyapunov matrix equation
!
!            op(A)'*X + X*op(A) = scale*C
!
!     where op(A) = A or A' (A**T), A is upper quasi-triangular and C is
!     symmetric (C = C'). (A' denotes the transpose of the matrix A.)
!     A is N-by-N, the right hand side C and the solution X are N-by-N,
!     and scale is an output scale factor, set less than or equal to 1
!     to avoid overflow in X. The solution matrix X is overwritten
!     onto C.
!
!     A must be in Schur canonical form (as returned by LAPACK routines
!     DGEES or DHSEQR), that is, block upper triangular with 1-by-1 and
!     2-by-2 diagonal blocks; each 2-by-2 diagonal block has its
!     diagonal elements equal and its off-diagonal elements of opposite
!     sign.
!
!     ARGUMENTS 
!
!     Mode Parameters
!
!     TRANA   CHARACTER*1
!             Specifies the form of op(A) to be used, as follows:
!             = 'N':  op(A) = A    (No transpose);
!             = 'T':  op(A) = A**T (Transpose);
!             = 'C':  op(A) = A**T (Conjugate transpose = Transpose).
!
!     Input/Output Parameters
!
!     N       (input) INTEGER
!             The order of the matrices A, X, and C.  N >= 0.
!
!     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!             The leading N-by-N part of this array must contain the
!             upper quasi-triangular matrix A, in Schur canonical form.
!             The part of A below the first sub-diagonal is not
!             referenced.
!
!     LDA     INTEGER
!             The leading dimension of array A.  LDA >= MAX(1,N).
!
!     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!             On entry, the leading N-by-N part of this array must  
!             contain the symmetric matrix C.
!             On exit, if INFO >= 0, the leading N-by-N part of this
!             array contains the symmetric solution matrix X.
!
!     LDC     INTEGER
!             The leading dimension of array C.  LDC >= MAX(1,N).
!
!     SCALE   (output) DOUBLE PRECISION
!             The scale factor, scale, set less than or equal to 1 to
!             prevent the solution overflowing.
!
!     Error Indicator
!
!     INFO    INTEGER
!             = 0:  successful exit;
!             < 0:  if INFO = -i, the i-th argument had an illegal
!                   value;
!             = 1:  if A and -A have common or very close eigenvalues;
!                   perturbed values were used to solve the equation
!                   (but the matrix A is unchanged).
!
!     METHOD
!
!     Bartels-Stewart algorithm is used. A set of equivalent linear
!     algebraic systems of equations of order at most four are formed
!     and solved using Gaussian elimination with complete pivoting.
!
!     REFERENCES
!
!     [1] Bartels, R.H. and Stewart, G.W.  T
!         Solution of the matrix equation A X + XB = C.
!         Comm. A.C.M., 15, pp. 820-826, 1972.
!
!     NUMERICAL ASPECTS
!                               3
!     The algorithm requires 0(N ) operations.
!
!     CONTRIBUTOR
!
!     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997.
!     Supersedes Release 2.0 routine SB03AY by Control Systems Research
!     Group, Kingston Polytechnic, United Kingdom, October 1982.
!     Based on DTRLYP by P. Petkov, Tech. University of Sofia, September
!     1993.
!
!     REVISIONS
!
!     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999.
!
!     KEYWORDS
!
!     Continuous-time system, Lyapunov equation, matrix algebra, real
!     Schur form.
!
!     ******************************************************************
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Scalar Arguments ..
      CHARACTER          TRANA
      INTEGER            INFO, LDA, LDC, N
      DOUBLE PRECISION   SCALE
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRNA
      INTEGER            UPLOSIGN
!     ..
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
!     .. Explicit Declarations
      INCLUDE 'RECSY_DECL.F'
!     ..
!     .. Executable Statements ..
!
!     Decode and Test input parameters.
!
      NOTRNA = LSAME( TRANA, 'N' )
!
      INFO = 0
      IF( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'T' ) .AND. &
                            .NOT.LSAME( TRANA, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDC.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB03MY', -INFO )
         RETURN
      END IF
!
      SCALE = ONE
!
!     Quick return if possible.
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF (.NOT.NOTRNA) THEN
         UPLOSIGN = 0
      ELSE
         UPLOSIGN = 1
      END IF
      CALL RECLYCT(UPLOSIGN, SCALE, N, A, LDA, C, LDC, INFO)
!
      RETURN
! *** Last line of SB03MY ***
      END
