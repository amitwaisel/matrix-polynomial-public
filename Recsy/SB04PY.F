!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!
      SUBROUTINE SB04PY( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, LDC, SCALE, DWORK, INFO )
!
!     MODIFIED RELEASE 4.0, WGS COPYRIGHT 2000.
!
!     PURPOSE
!
!     To solve for X the discrete-time Sylvester equation
!
!        op(A)*X*op(B) + ISGN*X = scale*C,
!
!     where op(A) = A or A**T, A and B are both upper quasi-triangular,
!     and ISGN = 1 or -1. A is M-by-M and B is N-by-N; the right hand 
!     side C and the solution X are M-by-N; and scale is an output scale
!     factor, set less than or equal to 1 to avoid overflow in X. The
!     solution matrix X is overwritten onto C.
!
!     A and B must be in Schur canonical form (as returned by LAPACK
!     Library routine DHSEQR), that is, block upper triangular with
!     1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block has
!     its diagonal elements equal and its off-diagonal elements of
!     opposite sign.
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
!     TRANB   CHARACTER*1
!             Specifies the form of op(B) to be used, as follows:
!             = 'N':  op(B) = B    (No transpose);
!             = 'T':  op(B) = B**T (Transpose);
!             = 'C':  op(B) = B**T (Conjugate transpose = Transpose).
!
!     ISGN    INTEGER
!             Specifies the sign of the equation as described before.
!             ISGN may only be 1 or -1.
!
!     Input/Output Parameters
!
!     M       (input) INTEGER
!             The order of the matrix A, and the number of rows in the
!             matrices X and C.  M >= 0.
!
!     N       (input) INTEGER
!             The order of the matrix B, and the number of columns in
!             the matrices X and C.  N >= 0.
!
!     A       (input) DOUBLE PRECISION array, dimension (LDA,M)
!             The leading M-by-M part of this array must contain the
!             upper quasi-triangular matrix A, in Schur canonical form.
!             The part of A below the first sub-diagonal is not
!             referenced.
!
!     LDA     INTEGER
!             The leading dimension of array A.  LDA >= MAX(1,M).
!
!     B       (input) DOUBLE PRECISION array, dimension (LDB,N)
!             The leading N-by-N part of this array must contain the
!             upper quasi-triangular matrix B, in Schur canonical form.
!             The part of B below the first sub-diagonal is not
!             referenced.
!
!     LDB     (input) INTEGER
!             The leading dimension of the array B.  LDB >= max(1,N).
!
!     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!             On entry, the leading M-by-N part of this array must  
!             contain the right hand side matrix C.
!             On exit, if INFO >= 0, the leading M-by-N part of this
!             array contains the solution matrix X.
!
!     LDC     INTEGER
!             The leading dimension of array C.  LDC >= MAX(1,M).
!
!     SCALE   (output) DOUBLE PRECISION
!             The scale factor, scale, set less than or equal to 1 to
!             prevent the solution overflowing.
!
!     Workspace
!
!     DWORK   DOUBLE PRECISION array, dimension (2*M)
!
!     Error Indicator
!
!     INFO    INTEGER
!             = 0:  successful exit;
!             < 0:  if INFO = -i, the i-th argument had an illegal
!                   value;
!             = 1:  A and -ISGN*B have almost reciprocal eigenvalues;
!                   perturbed values were used to solve the equation
!                   (but the matrices A and B are unchanged).
!
!     METHOD
!
!     The solution matrix X is computed column-wise via a back
!     substitution scheme, an extension and refinement of the algorithm
!     in [1], similar to that used in [2] for continuous-time Sylvester
!     equations. A set of equivalent linear algebraic systems of 
!     equations of order at most four are formed and solved using 
!     Gaussian elimination with complete pivoting.
!
!     REFERENCES
!
!     [1] Bartels, R.H. and Stewart, G.W.  T
!         Solution of the matrix equation A X + XB = C.
!         Comm. A.C.M., 15, pp. 820-826, 1972.
!
!     [2] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J.,
!         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A.,
!         Ostrouchov, S., and Sorensen, D.
!         LAPACK Users' Guide: Second Edition.
!         SIAM, Philadelphia, 1995.
!
!     NUMERICAL ASPECTS
!
!     The algorithm is stable and reliable, since Gaussian elimination
!     with complete pivoting is used.
!
!     CONTRIBUTORS
!
!     A. Varga, German Aerospace Center, Oberpfaffenhofen, March 2000.
!     D. Sima, University of Bucharest, April 2000.
!     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000.
!     Partly based on the routine SYLSV, A. Varga, 1992.
!
!     REVISIONS
!
!     -
!
!     KEYWORDS
!
!     Discrete-time system, matrix algebra, Sylvester equation.
!
!     ******************************************************************
!
      IMPLICIT NONE
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Scalar Arguments ..
      CHARACTER          TRANA, TRANB
      INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
      DOUBLE PRECISION   SCALE
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ), DWORK( * )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRNA, NOTRNB
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
!     Decode and Test input parameters
!
      NOTRNA = LSAME( TRANA, 'N' )
      NOTRNB = LSAME( TRANB, 'N' )
!
      INFO = 0
      IF( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'T' ) .AND. .NOT.LSAME( TRANA, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRNB .AND. .NOT.LSAME( TRANB, 'T' ) .AND. .NOT.LSAME( TRANB, 'C' ) ) THEN
         INFO = -2
      ELSE IF( ISGN.NE.1 .AND. ISGN.NE.-1 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB04PY', -INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
!
!     Set constants to control overflow.
!
      UPLOSIGN = (1-ISGN) / 2
      IF (.NOT.NOTRNA) UPLOSIGN = UPLOSIGN + 2
      IF (.NOT.NOTRNB) UPLOSIGN = UPLOSIGN + 4
      CALL RECSYDT(UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, INFO)
!
      RETURN
! *** Last line of SB04PY ***
      END
