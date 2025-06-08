!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!
      SUBROUTINE SG03AY( TRANS, N, A, LDA, E, LDE, X, LDX, SCALE, INFO )
!
!     MODIFIED RELEASE 4.0, WGS COPYRIGHT 1999.
!
!     PURPOSE
!
!     To solve for X either the reduced generalized continuous-time 
!     Lyapunov equation
!
!         T            T                     
!        A  * X * E + E  * X * A  =  SCALE * Y                       (1)
!
!     or
!
!                 T            T                      
!        A * X * E  + E * X * A   =  SCALE * Y                       (2)
!
!     where the right hand side Y is symmetric. A, E, Y, and the 
!     solution X are N-by-N matrices. The pencil A - lambda * E must be
!     in generalized Schur form (A upper quasitriangular, E upper 
!     triangular). SCALE is an output scale factor, set to avoid 
!     overflow in X.
!
!     ARGUMENTS
!
!     Mode Parameters
!    
!     TRANS   CHARACTER*1
!             Specifies whether the transposed equation is to be solved
!             or not:
!             = 'N':  Solve equation (1);    
!             = 'T':  Solve equation (2).
!
!     Input/Output Parameters
!   
!     N       (input) INTEGER
!             The order of the matrix A.  N >= 0.
!
!     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!             The leading N-by-N upper Hessenberg part of this array
!             must contain the quasitriangular matrix A.
!
!     LDA     INTEGER
!             The leading dimension of the array A.  LDA >= MAX(1,N).
!
!     E       (input) DOUBLE PRECISION array, dimension (LDE,N)
!             The leading N-by-N upper triangular part of this array
!             must contain the matrix E.
!
!     LDE     INTEGER
!             The leading dimension of the array E.  LDE >= MAX(1,N).
!
!     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N)
!             On entry, the leading N-by-N part of this array must
!             contain the right hand side matrix Y of the equation. Only
!             the upper triangular part of this matrix need be given.
!             On exit, the leading N-by-N part of this array contains
!             the solution matrix X of the equation.
!
!     LDX     INTEGER
!             The leading dimension of the array X.  LDX >= MAX(1,N).
!
!     SCALE   (output) DOUBLE PRECISION
!             The scale factor set to avoid overflow in X.
!             (0 < SCALE <= 1)
!
!     Error indicator1
! 
!     INFO    INTEGER
!             = 0:  successful exit;
!             < 0:  if INFO = -i, the i-th argument had an illegal
!                   value;
!             = 1:  equation is (almost) singular to working precision;
!                   perturbed values were used to solve the equation
!                   (but the matrices A and E are unchanged).
!
!     METHOD 
!
!     The solution X of (1) or (2) is computed via block back
!     substitution or block forward substitution, respectively. (See
!     [1] and [2] for details.)
!
!     REFERENCES
!
!     [1] Bartels, R.H., Stewart, G.W.
!         Solution of the equation A X + X B = C.
!         Comm. A.C.M., 15, pp. 820-826, 1972.
!       
!     [2] Penzl, T.
!         Numerical solution of generalized Lyapunov equations.
!         Advances in Comp. Math., vol. 8, pp. 33-48, 1998.
!
!     NUMERICAL ASPECTS
!
!     8/3 * N**3 flops are required by the routine. Note that we count a
!     single floating point arithmetic operation as one flop.
!
!     The algorithm is backward stable if the eigenvalues of the pencil 
!     A - lambda * E are real. Otherwise, linear systems of order at 
!     most 4 are involved into the computation. These systems are solved
!     by Gauss elimination with complete pivoting. The loss of stability
!     of the Gauss elimination with complete pivoting is rarely 
!     encountered in practice.
!
!     CONTRIBUTOR
!
!     T. Penzl, Technical University Chemnitz, Germany, Aug. 1998.
!
!     REVISIONS
!
!     Sep. 1998 (V. Sima).
!     Dec. 1998 (V. Sima).
!
!     KEYWORDS
!
!     Lyapunov equation
!
!     ******************************************************************
!
!     .. Parameters ..
      DOUBLE PRECISION  MONE, ONE, ZERO
      PARAMETER         ( MONE = -1.0D+0, ONE = 1.0D+0, ZERO = 0.0D+0 )
!     .. Scalar Arguments ..
      CHARACTER         TRANS
      DOUBLE PRECISION  SCALE
      INTEGER           INFO, LDA, LDE, LDX, N
!     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), E(LDE,*), X(LDX,*)
!     .. Local Scalars ..
      LOGICAL           NOTRNS
      INTEGER           UPLOSIGN
!     .. Local Arrays ..
!     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
!     .. External Subroutines ..
      EXTERNAL          XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC         MAX
!     .. Explicit Declarations
      INCLUDE 'RECSY_DECL.F'
!     .. Executable Statements ..
!
!     Decode input parameters.
!
      NOTRNS = LSAME( TRANS, 'N' )
!
!     Check the scalar input parameters.
!
      IF ( .NOT.( NOTRNS .OR. LSAME( TRANS, 'T' ) ) ) THEN
         INFO = -1
      ELSEIF ( N .LT. 0 ) THEN
         INFO = -2
      ELSEIF ( LDA .LT. MAX( 1, N ) ) THEN
         INFO = -4
      ELSEIF ( LDE .LT. MAX( 1, N ) ) THEN
         INFO = -6
      ELSEIF ( LDX .LT. MAX( 1, N ) ) THEN
         INFO = -8
      ELSE
         INFO = 0
      END IF
      IF ( INFO .NE. 0 ) THEN
         CALL XERBLA( 'SG03AY', -INFO )
         RETURN
      END IF 
!
      SCALE = ONE
!
!     Quick return if possible.
!
      IF ( N .EQ. 0 ) RETURN
!
      IF (.NOT.NOTRNS) THEN
         UPLOSIGN = 0
      ELSE
         UPLOSIGN = 1
      END IF
      CALL RECGLYCT(UPLOSIGN, SCALE, N, A, LDA, E, LDE, X, LDX, INFO)
      RETURN
! *** Last line of SG03AY ***
      END
