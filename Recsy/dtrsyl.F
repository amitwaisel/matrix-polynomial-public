!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!
      SUBROUTINE RECSY_DTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, &
                         LDC, SCALE, INFO )
!
!  -- LAPACK routine (MODIFIED version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          TRANA, TRANB
      INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
      DOUBLE PRECISION   SCALE
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  DTRSYL solves the real Sylvester matrix equation:
!
!     op(A)*X + X*op(B) = scale*C or
!     op(A)*X - X*op(B) = scale*C,
!
!  where op(A) = A or A**T, and  A and B are both upper quasi-
!  triangular. A is M-by-M and B is N-by-N; the right hand side C and
!  the solution X are M-by-N; and scale is an output scale factor, set
!  <= 1 to avoid overflow in X.
!
!  A and B must be in Schur canonical form (as returned by DHSEQR), that
!  is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
!  each 2-by-2 diagonal block has its diagonal elements equal and its
!  off-diagonal elements of opposite sign.
!
!  Arguments
!  =========
!
!  TRANA   (input) CHARACTER*1
!          Specifies the option op(A):
!          = 'N': op(A) = A    (No transpose)
!          = 'T': op(A) = A**T (Transpose)
!          = 'C': op(A) = A**H (Conjugate transpose = Transpose)
!
!  TRANB   (input) CHARACTER*1
!          Specifies the option op(B):
!          = 'N': op(B) = B    (No transpose)
!          = 'T': op(B) = B**T (Transpose)
!          = 'C': op(B) = B**H (Conjugate transpose = Transpose)
!
!  ISGN    (input) INTEGER
!          Specifies the sign in the equation:
!          = +1: solve op(A)*X + X*op(B) = scale*C
!          = -1: solve op(A)*X - X*op(B) = scale*C
!
!  M       (input) INTEGER
!          The order of the matrix A, and the number of rows in the
!          matrices X and C. M >= 0.
!
!  N       (input) INTEGER
!          The order of the matrix B, and the number of columns in the
!          matrices X and C. N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,M)
!          The upper quasi-triangular matrix A, in Schur canonical form.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,M).
!
!  B       (input) DOUBLE PRECISION array, dimension (LDB,N)
!          The upper quasi-triangular matrix B, in Schur canonical form.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B. LDB >= max(1,N).
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N right hand side matrix C.
!          On exit, C is overwritten by the solution matrix X.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M)
!
!  SCALE   (output) DOUBLE PRECISION
!          The scale factor, scale, set <= 1 to avoid overflow in X.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          = 1: A and B have common or very close eigenvalues; perturbed
!               values were used to solve the equation (but the matrices
!               A and B are unchanged).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
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
      PRINT*,'DTRSYL',M,N
      NOTRNA = LSAME( TRANA, 'N' )
      NOTRNB = LSAME( TRANB, 'N' )
!
      INFO = 0
      IF( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'T' ) .AND. .NOT. &
          LSAME( TRANA, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRNB .AND. .NOT.LSAME( TRANB, 'T' ) .AND. .NOT. &
               LSAME( TRANB, 'C' ) ) THEN
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
         CALL XERBLA( 'DTRSYL', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) &
         RETURN
!
!     Set constants to control overflow
!
      UPLOSIGN = (1-ISGN) / 2
      IF (.NOT.NOTRNA) UPLOSIGN = UPLOSIGN + 2
      IF (.NOT.NOTRNB) UPLOSIGN = UPLOSIGN + 4
      CALL RECSYCT(UPLOSIGN, SCALE, M, N, A, LDA, B, LDB, C, LDC, INFO)
!
      RETURN
!
!     End of DTRSYL
!
      END
