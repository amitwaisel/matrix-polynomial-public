!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Ume� University, Sweden
!
      SUBROUTINE RECSYDT_TR2(M, N, A, LDA, B, LDB, C, LDC, INFO, MACHINE)
!     Solves A*X*B^T + X = scale*C
!     A (M*M) and B (N*N) quasi-triangular
!     
!     Iff close eigenvalues are detected, INFO is non-zero and elements 
!     of the C (M*N) matrix are no longer valid
!
      IMPLICIT NONE
!
      INTENT(IN) M,N,A,LDA,B,LDB,LDC,MACHINE
      INTENT(OUT) INFO
      INTENT(INOUT) C
!     .. Scalar Arguments ..
      INTEGER          :: M, N, LDA, LDB, LDC, INFO
!     .. Array Arguments ..
      DOUBLE PRECISION :: A(0:LDA-1,0:*), B(0:LDB-1,0:*), C(0:LDC-1,0:*), MACHINE(0:*)
!     .. Local Scalars ..
      INTEGER          :: I, II, J, JJ
!     .. Local Arrays ..
      DOUBLE PRECISION :: D(0:3,0:3)
      DOUBLE PRECISION :: E(0:3)
      DOUBLE PRECISION :: T00,T01,T02,T03,T10,T11,T12,T13
!     .. Executable Statements ..
!
      INFO = 0
      I = M-1
      DO WHILE (I.GE.0) 
         IF (I.NE.0.AND.A(I,I-1).NE.0) THEN
            II = 2
         ELSE
            II = 1
         END IF
         I = I-II+1
         J = N-1
         DO WHILE (J.GE.0)
            IF (J.NE.0.AND.B(J,J-1).NE.0) THEN
               JJ = 2
            ELSE
               JJ = 1
            END IF
            J = J-JJ+1
            IF (II.EQ.1) THEN
               IF (JJ.EQ.1) THEN
                  CALL RECSY_SOLVE_11(C(I,J), A(I,I) * B(J,J) + 1, INFO, MACHINE)
                  IF (INFO.NE.0) RETURN
                  SELECT CASE (J)
                  CASE (1)
                     C(I,0) = C(I,0) - A(I,I) * C(I,1) * B(0,1)
                  CASE (2)                            
                     C(I,0) = C(I,0) - A(I,I) * C(I,2) * B(0,2)
                     C(I,1) = C(I,1) - A(I,I) * C(I,2) * B(1,2)
                  CASE (3)                            
                     C(I,0) = C(I,0) - A(I,I) * C(I,3) * B(0,3)
                     C(I,1) = C(I,1) - A(I,I) * C(I,3) * B(1,3)
                     C(I,2) = C(I,2) - A(I,I) * C(I,3) * B(2,3)
                  END SELECT
               ELSE
                  D(0,0) = A(I,I) * B(J,J) + 1
                  D(1,0) = A(I,I) * B(J+1,J)
                  D(2,0) = A(I,I) * B(J,J+1)
                  D(3,0) = A(I,I) * B(J+1,J+1) + 1
                  CALL RECSY_SOLVE_22(C(I,J), C(I,J+1), D, INFO, MACHINE)
                  IF (INFO.NE.0) RETURN
                  SELECT CASE (J)
                  CASE (1)
                     C(I,0) = C(I,0) - A(I,I) * (C(I,1) * B(0,1) + C(I,2) * B(0,2))
                  CASE (2)                            
                     C(I,0) = C(I,0) - A(I,I) * (C(I,2) * B(0,2) + C(I,3) * B(0,3))
                     C(I,1) = C(I,1) - A(I,I) * (C(I,2) * B(1,2) + C(I,3) * B(1,3))
                  END SELECT
               END IF
            ELSE
               IF (JJ.EQ.1) THEN
                  D(0,0) = A(I,I) * B(J,J) + 1
                  D(1,0) = A(I+1,I) * B(J,J)
                  D(2,0) = A(I,I+1) * B(J,J)
                  D(3,0) = A(I+1,I+1) * B(J,J) + 1
                  CALL RECSY_SOLVE_22(C(I,J), C(I+1,J), D, INFO, MACHINE)
                  IF (INFO.NE.0) RETURN
                  SELECT CASE (J)
                  CASE (1)
                     T00 = C(I,1) * B(0,1)
                     T10 = C(I+1,1) * B(0,1)
                     C(I,  0) = C(I,  0) - A(I,  I) * T00 - A(I,  I+1) * T10
                     C(I+1,0) = C(I+1,0) - A(I+1,I) * T00 - A(I+1,I+1) * T10
                  CASE (2)                            
                     T00 = C(I,2) * B(0,2)
                     T01 = C(I,2) * B(1,2)
                     T10 = C(I+1,2) * B(0,2)
                     T11 = C(I+1,2) * B(1,2)
                     C(I,  0) = C(I,  0) - A(I,  I) * T00 - A(I,  I+1) * T10
                     C(I+1,0) = C(I+1,0) - A(I+1,I) * T00 - A(I+1,I+1) * T10
                     C(I,  1) = C(I,  1) - A(I,  I) * T01 - A(I,  I+1) * T11
                     C(I+1,1) = C(I+1,1) - A(I+1,I) * T01 - A(I+1,I+1) * T11
                  CASE (3)                                      
                     T00 = C(I,3) * B(0,3)
                     T01 = C(I,3) * B(1,3)
                     T02 = C(I,3) * B(2,3)
                     T10 = C(I+1,3) * B(0,3)
                     T11 = C(I+1,3) * B(1,3)
                     T12 = C(I+1,3) * B(2,3)
                     C(I,  0) = C(I,  0) - A(I,  I) * T00 - A(I,  I+1) * T10
                     C(I+1,0) = C(I+1,0) - A(I+1,I) * T00 - A(I+1,I+1) * T10
                     C(I,  1) = C(I,  1) - A(I,  I) * T01 - A(I,  I+1) * T11
                     C(I+1,1) = C(I+1,1) - A(I+1,I) * T01 - A(I+1,I+1) * T11
                     C(I,  2) = C(I,  2) - A(I,  I) * T02 - A(I,  I+1) * T12
                     C(I+1,2) = C(I+1,2) - A(I+1,I) * T02 - A(I+1,I+1) * T12
                  END SELECT
               ELSE
                  D(0,0) = A(I,I)   * B(J,J) + 1
                  D(0,1) = A(I,I+1) * B(J,J)
                  D(0,2) = A(I,I)   * B(J,J+1)
                  D(0,3) = A(I,I+1) * B(J,J+1)
                  D(1,0) = A(I+1,I)   * B(J,J)
                  D(1,1) = A(I+1,I+1) * B(J,J) + 1
                  D(1,2) = A(I+1,I)   * B(J,J+1)
                  D(1,3) = A(I+1,I+1) * B(J,J+1)
                  D(2,0) = A(I,I)   * B(J+1,J)
                  D(2,1) = A(I,I+1) * B(J+1,J)
                  D(2,2) = A(I,I)   * B(J+1,J+1) + 1
                  D(2,3) = A(I,I+1) * B(J+1,J+1)
                  D(3,0) = A(I+1,I)   * B(J+1,J)
                  D(3,1) = A(I+1,I+1) * B(J+1,J)
                  D(3,2) = A(I+1,I)   * B(J+1,J+1)
                  D(3,3) = A(I+1,I+1) * B(J+1,J+1) + 1
                  E(0) = C(I,J)
                  E(1) = C(I+1,J)
                  E(2) = C(I,J+1)
                  E(3) = C(I+1,J+1)
                  CALL RECSY_SOLVE_44(E, D, INFO, MACHINE)
                  IF (INFO.EQ.1) RETURN
                  C(I,J) = E(0)
                  C(I+1,J) = E(1)
                  C(I,J+1) = E(2)
                  C(I+1,J+1) = E(3)
                  SELECT CASE (J)
                  CASE (1)
                     T00 = C(I  ,1) * B(0,1) + C(I  ,2) * B(0,2)
                     T10 = C(I+1,1) * B(0,1) + C(I+1,2) * B(0,2)
                     C(I,  0) = C(I,  0) - A(I,  I) * T00 - A(I,  I+1) * T10
                     C(I+1,0) = C(I+1,0) - A(I+1,I) * T00 - A(I+1,I+1) * T10
                  CASE (2)                            
                     T00 = C(I  ,2) * B(0,2) + C(I  ,3) * B(0,3)
                     T01 = C(I  ,2) * B(1,2) + C(I  ,3) * B(1,3)
                     T10 = C(I+1,2) * B(0,2) + C(I+1,3) * B(0,3)
                     T11 = C(I+1,2) * B(1,2) + C(I+1,3) * B(1,3)
                     C(I,  0) = C(I,  0) - A(I,  I) * T00 - A(I,  I+1) * T10
                     C(I+1,0) = C(I+1,0) - A(I+1,I) * T00 - A(I+1,I+1) * T10
                     C(I,  1) = C(I,  1) - A(I,  I) * T01 - A(I,  I+1) * T11
                     C(I+1,1) = C(I+1,1) - A(I+1,I) * T01 - A(I+1,I+1) * T11
                  END SELECT
               END IF
            END IF
            J = J-1
         END DO
         IF (II.EQ.1) THEN
            IF (I.GT.0) THEN
               SELECT CASE (N)
               CASE (1)
                  T00 = C(I,0) * B(0,0)
               CASE (2)
                  T00 = C(I,0) * B(0,0) + C(I,1) * B(0,1)
                  T01 = C(I,0) * B(1,0) + C(I,1) * B(1,1)
               CASE (3)
                  T00 = C(I,0) * B(0,0) + C(I,1) * B(0,1) + C(I,2) * B(0,2)
                  T01 = C(I,0) * B(1,0) + C(I,1) * B(1,1) + C(I,2) * B(1,2)
                  T02 =                   C(I,1) * B(2,1) + C(I,2) * B(2,2)
               CASE (4)
                  T00 = C(I,0) * B(0,0) + C(I,1) * B(0,1) + C(I,2) * B(0,2) + C(I,3) * B(0,3)
                  T01 = C(I,0) * B(1,0) + C(I,1) * B(1,1) + C(I,2) * B(1,2) + C(I,3) * B(1,3)
                  T02 =                   C(I,1) * B(2,1) + C(I,2) * B(2,2) + C(I,3) * B(2,3)
                  T03 =                                     C(I,2) * B(3,2) + C(I,3) * B(3,3)
               END SELECT
            END IF

            SELECT CASE (I)
            CASE (1)
               SELECT CASE (N)
               CASE (1)
                  C(0,0) = C(0,0) - A(0,1) * T00
               CASE (2)                            
                  C(0,0) = C(0,0) - A(0,1) * T00
                  C(0,1) = C(0,1) - A(0,1) * T01
               CASE (3)                            
                  C(0,0) = C(0,0) - A(0,1) * T00
                  C(0,1) = C(0,1) - A(0,1) * T01
                  C(0,2) = C(0,2) - A(0,1) * T02
               CASE (4)                            
                  C(0,0) = C(0,0) - A(0,1) * T00
                  C(0,1) = C(0,1) - A(0,1) * T01
                  C(0,2) = C(0,2) - A(0,1) * T02
                  C(0,3) = C(0,3) - A(0,1) * T03
               END SELECT
            CASE (2)
               SELECT CASE (N)
               CASE (1)
                  C(0,0) = C(0,0) - A(0,2) * T00
                  C(1,0) = C(1,0) - A(1,2) * T00
               CASE (2)                            
                  C(0,0) = C(0,0) - A(0,2) * T00
                  C(0,1) = C(0,1) - A(0,2) * T01
                  C(1,0) = C(1,0) - A(1,2) * T00
                  C(1,1) = C(1,1) - A(1,2) * T01
               CASE (3)                            
                  C(0,0) = C(0,0) - A(0,2) * T00
                  C(0,1) = C(0,1) - A(0,2) * T01
                  C(0,2) = C(0,2) - A(0,2) * T02
                  C(1,0) = C(1,0) - A(1,2) * T00
                  C(1,1) = C(1,1) - A(1,2) * T01
                  C(1,2) = C(1,2) - A(1,2) * T02
               CASE (4)                            
                  C(0,0) = C(0,0) - A(0,2) * T00
                  C(0,1) = C(0,1) - A(0,2) * T01
                  C(0,2) = C(0,2) - A(0,2) * T02
                  C(0,3) = C(0,3) - A(0,2) * T03
                  C(1,0) = C(1,0) - A(1,2) * T00
                  C(1,1) = C(1,1) - A(1,2) * T01
                  C(1,2) = C(1,2) - A(1,2) * T02
                  C(1,3) = C(1,3) - A(1,2) * T03
               END SELECT
            CASE (3)
               SELECT CASE (N)
               CASE (1)
                  C(0,0) = C(0,0) - A(0,3) * T00
                  C(1,0) = C(1,0) - A(1,3) * T00
                  C(2,0) = C(2,0) - A(2,3) * T00
               CASE (2)                            
                  C(0,0) = C(0,0) - A(0,3) * T00
                  C(0,1) = C(0,1) - A(0,3) * T01
                  C(1,0) = C(1,0) - A(1,3) * T00
                  C(1,1) = C(1,1) - A(1,3) * T01
                  C(2,0) = C(2,0) - A(2,3) * T00
                  C(2,1) = C(2,1) - A(2,3) * T01
               CASE (3)                            
                  C(0,0) = C(0,0) - A(0,3) * T00
                  C(0,1) = C(0,1) - A(0,3) * T01
                  C(0,2) = C(0,2) - A(0,3) * T02
                  C(1,0) = C(1,0) - A(1,3) * T00
                  C(1,1) = C(1,1) - A(1,3) * T01
                  C(1,2) = C(1,2) - A(1,3) * T02
                  C(2,0) = C(2,0) - A(2,3) * T00
                  C(2,1) = C(2,1) - A(2,3) * T01
                  C(2,2) = C(2,2) - A(2,3) * T02
               CASE (4)                            
                  C(0,0) = C(0,0) - A(0,3) * T00
                  C(0,1) = C(0,1) - A(0,3) * T01
                  C(0,2) = C(0,2) - A(0,3) * T02
                  C(0,3) = C(0,3) - A(0,3) * T03
                  C(1,0) = C(1,0) - A(1,3) * T00
                  C(1,1) = C(1,1) - A(1,3) * T01
                  C(1,2) = C(1,2) - A(1,3) * T02
                  C(1,3) = C(1,3) - A(1,3) * T03
                  C(2,0) = C(2,0) - A(2,3) * T00
                  C(2,1) = C(2,1) - A(2,3) * T01
                  C(2,2) = C(2,2) - A(2,3) * T02
                  C(2,3) = C(2,3) - A(2,3) * T03
               END SELECT
            END SELECT
         ELSE
            IF (I.GT.0) THEN
               SELECT CASE (N)
               CASE (1)
                  T00 = C(I,0)   * B(0,0)
                  T10 = C(I+1,0) * B(0,0)
               CASE (2)
                  T00 = C(I,0)   * B(0,0) + C(I,1)   * B(0,1)
                  T01 = C(I,0)   * B(1,0) + C(I,1)   * B(1,1)
                  T10 = C(I+1,0) * B(0,0) + C(I+1,1) * B(0,1)
                  T11 = C(I+1,0) * B(1,0) + C(I+1,1) * B(1,1)
               CASE (3)
                  T00 = C(I,0)   * B(0,0) + C(I,1)   * B(0,1) + C(I,2)   * B(0,2)
                  T01 = C(I,0)   * B(1,0) + C(I,1)   * B(1,1) + C(I,2)   * B(1,2)
                  T02 =                     C(I,1)   * B(2,1) + C(I,2)   * B(2,2)
                  T10 = C(I+1,0) * B(0,0) + C(I+1,1) * B(0,1) + C(I+1,2) * B(0,2)
                  T11 = C(I+1,0) * B(1,0) + C(I+1,1) * B(1,1) + C(I+1,2) * B(1,2)
                  T12 =                     C(I+1,1) * B(2,1) + C(I+1,2) * B(2,2)
               CASE (4)
                  T00 = C(I,0)   * B(0,0) + C(I,1)   * B(0,1) + C(I,2)   * B(0,2) + C(I,3)   * B(0,3)
                  T01 = C(I,0)   * B(1,0) + C(I,1)   * B(1,1) + C(I,2)   * B(1,2) + C(I,3)   * B(1,3)
                  T02 =                     C(I,1)   * B(2,1) + C(I,2)   * B(2,2) + C(I,3)   * B(2,3)
                  T03 =                                         C(I,2)   * B(3,2) + C(I,3)   * B(3,3)
                  T10 = C(I+1,0) * B(0,0) + C(I+1,1) * B(0,1) + C(I+1,2) * B(0,2) + C(I+1,3) * B(0,3)
                  T11 = C(I+1,0) * B(1,0) + C(I+1,1) * B(1,1) + C(I+1,2) * B(1,2) + C(I+1,3) * B(1,3)
                  T12 =                     C(I+1,1) * B(2,1) + C(I+1,2) * B(2,2) + C(I+1,3) * B(2,3)
                  T13 =                                         C(I+1,2) * B(3,2) + C(I+1,3) * B(3,3)
               END SELECT
            END IF

            SELECT CASE (I)
            CASE (1)
               SELECT CASE (N)
               CASE (1)
                  C(0,0) = C(0,0) - A(0,1) * T00 - A(0,2) * T10
               CASE (2)                            
                  C(0,0) = C(0,0) - A(0,1) * T00 - A(0,2) * T10
                  C(0,1) = C(0,1) - A(0,1) * T01 - A(0,2) * T11
               CASE (3)                            
                  C(0,0) = C(0,0) - A(0,1) * T00 - A(0,2) * T10
                  C(0,1) = C(0,1) - A(0,1) * T01 - A(0,2) * T11
                  C(0,2) = C(0,2) - A(0,1) * T02 - A(0,2) * T12
               CASE (4)                            
                  C(0,0) = C(0,0) - A(0,1) * T00 - A(0,2) * T10
                  C(0,1) = C(0,1) - A(0,1) * T01 - A(0,2) * T11
                  C(0,2) = C(0,2) - A(0,1) * T02 - A(0,2) * T12
                  C(0,3) = C(0,3) - A(0,1) * T03 - A(0,2) * T13
               END SELECT
            CASE (2)
               SELECT CASE (N)
               CASE (1)
                  C(0,0) = C(0,0) - A(0,2) * T00 - A(0,3) * T10
                  C(1,0) = C(1,0) - A(1,2) * T00 - A(1,3) * T10
               CASE (2)                                           
                  C(0,0) = C(0,0) - A(0,2) * T00 - A(0,3) * T10
                  C(0,1) = C(0,1) - A(0,2) * T01 - A(0,3) * T11
                  C(1,0) = C(1,0) - A(1,2) * T00 - A(1,3) * T10
                  C(1,1) = C(1,1) - A(1,2) * T01 - A(1,3) * T11
               CASE (3)                                           
                  C(0,0) = C(0,0) - A(0,2) * T00 - A(0,3) * T10
                  C(0,1) = C(0,1) - A(0,2) * T01 - A(0,3) * T11
                  C(0,2) = C(0,2) - A(0,2) * T02 - A(0,3) * T12
                  C(1,0) = C(1,0) - A(1,2) * T00 - A(1,3) * T10
                  C(1,1) = C(1,1) - A(1,2) * T01 - A(1,3) * T11
                  C(1,2) = C(1,2) - A(1,2) * T02 - A(1,3) * T12
               CASE (4)                                           
                  C(0,0) = C(0,0) - A(0,2) * T00 - A(0,3) * T10
                  C(0,1) = C(0,1) - A(0,2) * T01 - A(0,3) * T11
                  C(0,2) = C(0,2) - A(0,2) * T02 - A(0,3) * T12
                  C(0,3) = C(0,3) - A(0,2) * T03 - A(0,3) * T13
                  C(1,0) = C(1,0) - A(1,2) * T00 - A(1,3) * T10
                  C(1,1) = C(1,1) - A(1,2) * T01 - A(1,3) * T11
                  C(1,2) = C(1,2) - A(1,2) * T02 - A(1,3) * T12
                  C(1,3) = C(1,3) - A(1,2) * T03 - A(1,3) * T13
               END SELECT
            END SELECT
         END IF
         I = I-1
      END DO
      RETURN
      END SUBROUTINE


