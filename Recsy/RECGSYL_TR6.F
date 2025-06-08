!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      SUBROUTINE RECGSYL_TR6(M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, INFO, MACHINE)
!     Solves A^T*X*B^T + C^T*X*D^T=E
!     A (M*M) and B (N*N) quasi-triangular
!     C (M*M) and D (N*N) triangular
!     Iff close eigenvalues are detected, INFO is non-zero and elements 
!     of the E (M*N) matrix are no longer valid
!
      IMPLICIT NONE
!
      INTENT(IN) M,N,A,LDA,B,LDB,C,LDC,D,LDD,LDE,MACHINE
      INTENT(OUT) INFO
      INTENT(INOUT) E
!     .. Scalar Arguments ..
      INTEGER          :: M, N, LDA, LDB, LDC, LDD, LDE, INFO
!     .. Array Arguments ..
      DOUBLE PRECISION :: A(0:LDA-1,0:*), B(0:LDB-1,0:*), C(0:LDC-1,0:*), D(0:LDD-1,0:*), E(0:LDE-1,0:*), MACHINE(0:*)
!     .. Local Scalars ..
      INTEGER          :: I, II, J, JJ
!     .. Local Arrays ..
      DOUBLE PRECISION :: G(0:3,0:3)
      DOUBLE PRECISION :: H(0:3)
      DOUBLE PRECISION :: T00,T01,T02,T03,T10,T11,T12,T13
      DOUBLE PRECISION :: U00,U01,U02,U03,U10,U11,U12,U13
!     .. Executable Statements ..
!
      INFO = 0
      I = 0
      DO WHILE (I.LT.M) 
         IF (I.NE.M-1.AND.A(I+1,I).NE.0) THEN
            II = 2
         ELSE
            II = 1
         END IF
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
                  CALL RECSY_SOLVE_11(E(I,J), A(I,I) * B(J,J) + C(I,I) * D(J,J), INFO, MACHINE)
                  IF (INFO.NE.0) RETURN
                  SELECT CASE (J)
                  CASE (1)
                     E(I,0) = E(I,0) - A(I,I) * E(I,1) * B(0,1) - C(I,I) * E(I,1) * D(0,1)
                  CASE (2)                                              
                     E(I,0) = E(I,0) - A(I,I) * E(I,2) * B(0,2) - C(I,I) * E(I,2) * D(0,2)
                     E(I,1) = E(I,1) - A(I,I) * E(I,2) * B(1,2) - C(I,I) * E(I,2) * D(1,2)
                  CASE (3)                                              
                     E(I,0) = E(I,0) - A(I,I) * E(I,3) * B(0,3) - C(I,I) * E(I,3) * D(0,3)
                     E(I,1) = E(I,1) - A(I,I) * E(I,3) * B(1,3) - C(I,I) * E(I,3) * D(1,3)
                     E(I,2) = E(I,2) - A(I,I) * E(I,3) * B(2,3) - C(I,I) * E(I,3) * D(2,3)
                  END SELECT
               ELSE
                  G(0,0) = A(I,I) * B(J  ,J  ) + C(I,I) * D(J  ,J  )
                  G(1,0) = A(I,I) * B(J+1,J  )
                  G(2,0) = A(I,I) * B(J  ,J+1) + C(I,I) * D(J  ,J+1)
                  G(3,0) = A(I,I) * B(J+1,J+1) + C(I,I) * D(J+1,J+1)
                  CALL RECSY_SOLVE_22(E(I,J), E(I,J+1), G, INFO, MACHINE)
                  IF (INFO.NE.0) RETURN
                  SELECT CASE (J)
                  CASE (1)
                     E(I,0) = E(I,0) - A(I,I) * (E(I,1) * B(0,1) + E(I,2) * B(0,2)) - C(I,I) * (E(I,1) * D(0,1) + E(I,2) * D(0,2))
                  CASE (2)                                              
                     E(I,0) = E(I,0) - A(I,I) * (E(I,2) * B(0,2) + E(I,3) * B(0,3)) - C(I,I) * (E(I,2) * D(0,2) + E(I,3) * D(0,3))
                     E(I,1) = E(I,1) - A(I,I) * (E(I,2) * B(1,2) + E(I,3) * B(1,3)) - C(I,I) * (E(I,2) * D(1,2) + E(I,3) * D(1,3))
                  END SELECT
               END IF
            ELSE
               IF (JJ.EQ.1) THEN
                  G(0,0) = A(I  ,I  ) * B(J,J) + C(I  ,I  ) * D(J,J)
                  G(1,0) = A(I  ,I+1) * B(J,J) + C(I  ,I+1) * D(J,J)
                  G(2,0) = A(I+1,I  ) * B(J,J)
                  G(3,0) = A(I+1,I+1) * B(J,J) + C(I+1,I+1) * D(J,J)
                  CALL RECSY_SOLVE_22(E(I,J), E(I+1,J), G, INFO, MACHINE)
                  IF (INFO.NE.0) RETURN
                  SELECT CASE (J)
                  CASE (1)
                     T00 = E(I  ,1) * B(0,1)
                     T10 = E(I+1,1) * B(0,1)
                     U00 = E(I  ,1) * D(0,1)
                     U10 = E(I+1,1) * D(0,1)
                     E(I,  0) = E(I,  0) - A(I,I  ) * T00 - A(I+1,I  ) * T10 - C(I,I  ) * U00
                     E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - A(I+1,I+1) * T10 - C(I,I+1) * U00 - C(I+1,I+1) * U10
                  CASE (2)                            
                     T00 = E(I  ,2) * B(0,2)
                     T01 = E(I  ,2) * B(1,2)
                     T10 = E(I+1,2) * B(0,2)
                     T11 = E(I+1,2) * B(1,2)
                     U00 = E(I  ,2) * D(0,2)
                     U01 = E(I  ,2) * D(1,2)
                     U10 = E(I+1,2) * D(0,2)
                     U11 = E(I+1,2) * D(1,2)
                     E(I,  0) = E(I,  0) - A(I,I  ) * T00 - A(I+1,I  ) * T10 - C(I,I  ) * U00
                     E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - A(I+1,I+1) * T10 - C(I,I+1) * U00 - C(I+1,I+1) * U10
                     E(I,  1) = E(I,  1) - A(I,I  ) * T01 - A(I+1,I  ) * T11 - C(I,I  ) * U01
                     E(I+1,1) = E(I+1,1) - A(I,I+1) * T01 - A(I+1,I+1) * T11 - C(I,I+1) * U01 - C(I+1,I+1) * U11
                  CASE (3)                                      
                     T00 = E(I  ,3) * B(0,3)
                     T01 = E(I  ,3) * B(1,3)
                     T02 = E(I  ,3) * B(2,3)
                     T10 = E(I+1,3) * B(0,3)
                     T11 = E(I+1,3) * B(1,3)
                     T12 = E(I+1,3) * B(2,3)
                     U00 = E(I  ,3) * D(0,3)
                     U01 = E(I  ,3) * D(1,3)
                     U02 = E(I  ,3) * D(2,3)
                     U10 = E(I+1,3) * D(0,3)
                     U11 = E(I+1,3) * D(1,3)
                     U12 = E(I+1,3) * D(2,3)
                     E(I,  0) = E(I,  0) - A(I,I  ) * T00 - A(I+1,I  ) * T10 - C(I,I  ) * U00
                     E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - A(I+1,I+1) * T10 - C(I,I+1) * U00 - C(I+1,I+1) * U10
                     E(I,  1) = E(I,  1) - A(I,I  ) * T01 - A(I+1,I  ) * T11 - C(I,I  ) * U01
                     E(I+1,1) = E(I+1,1) - A(I,I+1) * T01 - A(I+1,I+1) * T11 - C(I,I+1) * U01 - C(I+1,I+1) * U11
                     E(I,  2) = E(I,  2) - A(I,I  ) * T02 - A(I+1,I  ) * T12 - C(I,I  ) * U02
                     E(I+1,2) = E(I+1,2) - A(I,I+1) * T02 - A(I+1,I+1) * T12 - C(I,I+1) * U02 - C(I+1,I+1) * U12
                  END SELECT
               ELSE
                  G(0,0) = A(I  ,I  ) * B(J  ,J  ) + C(I  ,I  ) * D(J  ,J  )
                  G(0,1) = A(I+1,I  ) * B(J  ,J  )
                  G(0,2) = A(I  ,I  ) * B(J  ,J+1) + C(I  ,I  ) * D(J  ,J+1)
                  G(0,3) = A(I+1,I  ) * B(J  ,J+1)
                  G(1,0) = A(I  ,I+1) * B(J  ,J  ) + C(I  ,I+1) * D(J  ,J  )
                  G(1,1) = A(I+1,I+1) * B(J  ,J  ) + C(I+1,I+1) * D(J  ,J  )
                  G(1,2) = A(I  ,I+1) * B(J  ,J+1) + C(I  ,I+1) * D(J  ,J+1)
                  G(1,3) = A(I+1,I+1) * B(J  ,J+1) + C(I+1,I+1) * D(J  ,J+1)
                  G(2,0) = A(I  ,I  ) * B(J+1,J  ) + C(I  ,I  ) * D(J+1,J  )
                  G(2,1) = A(I+1,I  ) * B(J+1,J  )
                  G(2,2) = A(I  ,I  ) * B(J+1,J+1) + C(I  ,I  ) * D(J+1,J+1)
                  G(2,3) = A(I+1,I  ) * B(J+1,J+1)
                  G(3,0) = A(I  ,I+1) * B(J+1,J  )
                  G(3,1) = A(I+1,I+1) * B(J+1,J  )
                  G(3,2) = A(I  ,I+1) * B(J+1,J+1) + C(I  ,I+1) * D(J+1,J+1)
                  G(3,3) = A(I+1,I+1) * B(J+1,J+1) + C(I+1,I+1) * D(J+1,J+1)
                  H(0) = E(I,J)
                  H(1) = E(I+1,J)
                  H(2) = E(I,J+1)
                  H(3) = E(I+1,J+1)
                  CALL RECSY_SOLVE_44(H, G, INFO, MACHINE)
                  IF (INFO.EQ.1) RETURN
                  E(I,J) = H(0)
                  E(I+1,J) = H(1)
                  E(I,J+1) = H(2)
                  E(I+1,J+1) = H(3)
                  SELECT CASE (J)
                  CASE (1)
                     T00 = E(I  ,1) * B(0,1) + E(I  ,2) * B(0,2)
                     T10 = E(I+1,1) * B(0,1) + E(I+1,2) * B(0,2)
                     U00 = E(I  ,1) * D(0,1) + E(I  ,2) * D(0,2)
                     U10 = E(I+1,1) * D(0,1) + E(I+1,2) * D(0,2)
                     E(I,  0) = E(I,  0) - A(I,I  ) * T00 - A(I+1,I  ) * T10 - C(I,I  ) * U00
                     E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - A(I+1,I+1) * T10 - C(I,I+1) * U00 - C(I+1,I+1) * U10
                  CASE (2)                            
                     T00 = E(I  ,2) * B(0,2) + E(I  ,3) * B(0,3)
                     T01 = E(I  ,2) * B(1,2) + E(I  ,3) * B(1,3)
                     T10 = E(I+1,2) * B(0,2) + E(I+1,3) * B(0,3)
                     T11 = E(I+1,2) * B(1,2) + E(I+1,3) * B(1,3)
                     U00 = E(I  ,2) * D(0,2) + E(I  ,3) * D(0,3)
                     U01 = E(I  ,2) * D(1,2) + E(I  ,3) * D(1,3)
                     U10 = E(I+1,2) * D(0,2) + E(I+1,3) * D(0,3)
                     U11 = E(I+1,2) * D(1,2) + E(I+1,3) * D(1,3)
                     E(I,  0) = E(I,  0) - A(I,I  ) * T00 - A(I+1,I  ) * T10 - C(I,I  ) * U00
                     E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - A(I+1,I+1) * T10 - C(I,I+1) * U00 - C(I+1,I+1) * U10
                     E(I,  1) = E(I,  1) - A(I,I  ) * T01 - A(I+1,I  ) * T11 - C(I,I  ) * U01
                     E(I+1,1) = E(I+1,1) - A(I,I+1) * T01 - A(I+1,I+1) * T11 - C(I,I+1) * U01 - C(I+1,I+1) * U11
                  END SELECT
               END IF
            END IF
            J = J-1
         END DO
         IF (II.EQ.1) THEN
            IF (M-I-1.GT.0) THEN
               SELECT CASE (N)
               CASE (1)
                  T00 = E(I,0) * B(0,0)
                  U00 = E(I,0) * D(0,0)
               CASE (2)
                  T00 = E(I,0) * B(0,0) + E(I,1) * B(0,1)
                  T01 = E(I,0) * B(1,0) + E(I,1) * B(1,1)
                  U00 = E(I,0) * D(0,0) + E(I,1) * D(0,1)
                  U01 =                   E(I,1) * D(1,1)
               CASE (3)
                  T00 = E(I,0) * B(0,0) + E(I,1) * B(0,1) + E(I,2) * B(0,2)
                  T01 = E(I,0) * B(1,0) + E(I,1) * B(1,1) + E(I,2) * B(1,2)
                  T02 =                   E(I,1) * B(2,1) + E(I,2) * B(2,2)
                  U00 = E(I,0) * D(0,0) + E(I,1) * D(0,1) + E(I,2) * D(0,2)
                  U01 =                   E(I,1) * D(1,1) + E(I,2) * D(1,2)
                  U02 =                                     E(I,2) * D(2,2)
               CASE (4)
                  T00 = E(I,0) * B(0,0) + E(I,1) * B(0,1) + E(I,2) * B(0,2) + E(I,3) * B(0,3)
                  T01 = E(I,0) * B(1,0) + E(I,1) * B(1,1) + E(I,2) * B(1,2) + E(I,3) * B(1,3)
                  T02 =                   E(I,1) * B(2,1) + E(I,2) * B(2,2) + E(I,3) * B(2,3)
                  T03 =                                     E(I,2) * B(3,2) + E(I,3) * B(3,3)
                  U00 = E(I,0) * D(0,0) + E(I,1) * D(0,1) + E(I,2) * D(0,2) + E(I,3) * D(0,3)
                  U01 =                   E(I,1) * D(1,1) + E(I,2) * D(1,2) + E(I,3) * D(1,3)
                  U02 =                                     E(I,2) * D(2,2) + E(I,3) * D(2,3)
                  U03 =                                                       E(I,3) * D(3,3)
               END SELECT
            END IF

            SELECT CASE (M-I-1)
            CASE (1)
               SELECT CASE (N)
               CASE (1)
                  E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - C(I,I+1) * U00
               CASE (2)                                          
                  E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - C(I,I+1) * U00
                  E(I+1,1) = E(I+1,1) - A(I,I+1) * T01 - C(I,I+1) * U01
               CASE (3)                                          
                  E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - C(I,I+1) * U00
                  E(I+1,1) = E(I+1,1) - A(I,I+1) * T01 - C(I,I+1) * U01
                  E(I+1,2) = E(I+1,2) - A(I,I+1) * T02 - C(I,I+1) * U02
               CASE (4)                                          
                  E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - C(I,I+1) * U00
                  E(I+1,1) = E(I+1,1) - A(I,I+1) * T01 - C(I,I+1) * U01
                  E(I+1,2) = E(I+1,2) - A(I,I+1) * T02 - C(I,I+1) * U02
                  E(I+1,3) = E(I+1,3) - A(I,I+1) * T03 - C(I,I+1) * U03
               END SELECT
            CASE (2)
               SELECT CASE (N)
               CASE (1)
                  E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - C(I,I+1) * U00
                  E(I+2,0) = E(I+2,0) - A(I,I+2) * T00 - C(I,I+2) * U00
               CASE (2)                                          
                  E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - C(I,I+1) * U00
                  E(I+1,1) = E(I+1,1) - A(I,I+1) * T01 - C(I,I+1) * U01
                  E(I+2,0) = E(I+2,0) - A(I,I+2) * T00 - C(I,I+2) * U00
                  E(I+2,1) = E(I+2,1) - A(I,I+2) * T01 - C(I,I+2) * U01
               CASE (3)                                          
                  E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - C(I,I+1) * U00
                  E(I+1,1) = E(I+1,1) - A(I,I+1) * T01 - C(I,I+1) * U01
                  E(I+1,2) = E(I+1,2) - A(I,I+1) * T02 - C(I,I+1) * U02
                  E(I+2,0) = E(I+2,0) - A(I,I+2) * T00 - C(I,I+2) * U00
                  E(I+2,1) = E(I+2,1) - A(I,I+2) * T01 - C(I,I+2) * U01
                  E(I+2,2) = E(I+2,2) - A(I,I+2) * T02 - C(I,I+2) * U02
               CASE (4)                                          
                  E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - C(I,I+1) * U00
                  E(I+1,1) = E(I+1,1) - A(I,I+1) * T01 - C(I,I+1) * U01
                  E(I+1,2) = E(I+1,2) - A(I,I+1) * T02 - C(I,I+1) * U02
                  E(I+1,3) = E(I+1,3) - A(I,I+1) * T03 - C(I,I+1) * U03
                  E(I+2,0) = E(I+2,0) - A(I,I+2) * T00 - C(I,I+2) * U00
                  E(I+2,1) = E(I+2,1) - A(I,I+2) * T01 - C(I,I+2) * U01
                  E(I+2,2) = E(I+2,2) - A(I,I+2) * T02 - C(I,I+2) * U02
                  E(I+2,3) = E(I+2,3) - A(I,I+2) * T03 - C(I,I+2) * U03
               END SELECT
            CASE (3)
               SELECT CASE (N)
               CASE (1)
                  E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - C(I,I+1) * U00
                  E(I+2,0) = E(I+2,0) - A(I,I+2) * T00 - C(I,I+2) * U00
                  E(I+3,0) = E(I+3,0) - A(I,I+3) * T00 - C(I,I+3) * U00
               CASE (2)                                          
                  E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - C(I,I+1) * U00
                  E(I+1,1) = E(I+1,1) - A(I,I+1) * T01 - C(I,I+1) * U01
                  E(I+2,0) = E(I+2,0) - A(I,I+2) * T00 - C(I,I+2) * U00
                  E(I+2,1) = E(I+2,1) - A(I,I+2) * T01 - C(I,I+2) * U01
                  E(I+3,0) = E(I+3,0) - A(I,I+3) * T00 - C(I,I+3) * U00
                  E(I+3,1) = E(I+3,1) - A(I,I+3) * T01 - C(I,I+3) * U01
               CASE (3)                                          
                  E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - C(I,I+1) * U00
                  E(I+1,1) = E(I+1,1) - A(I,I+1) * T01 - C(I,I+1) * U01
                  E(I+1,2) = E(I+1,2) - A(I,I+1) * T02 - C(I,I+1) * U02
                  E(I+2,0) = E(I+2,0) - A(I,I+2) * T00 - C(I,I+2) * U00
                  E(I+2,1) = E(I+2,1) - A(I,I+2) * T01 - C(I,I+2) * U01
                  E(I+2,2) = E(I+2,2) - A(I,I+2) * T02 - C(I,I+2) * U02
                  E(I+3,0) = E(I+3,0) - A(I,I+3) * T00 - C(I,I+3) * U00
                  E(I+3,1) = E(I+3,1) - A(I,I+3) * T01 - C(I,I+3) * U01
                  E(I+3,2) = E(I+3,2) - A(I,I+3) * T02 - C(I,I+3) * U02
               CASE (4)                                          
                  E(I+1,0) = E(I+1,0) - A(I,I+1) * T00 - C(I,I+1) * U00
                  E(I+1,1) = E(I+1,1) - A(I,I+1) * T01 - C(I,I+1) * U01
                  E(I+1,2) = E(I+1,2) - A(I,I+1) * T02 - C(I,I+1) * U02
                  E(I+1,3) = E(I+1,3) - A(I,I+1) * T03 - C(I,I+1) * U03
                  E(I+2,0) = E(I+2,0) - A(I,I+2) * T00 - C(I,I+2) * U00
                  E(I+2,1) = E(I+2,1) - A(I,I+2) * T01 - C(I,I+2) * U01
                  E(I+2,2) = E(I+2,2) - A(I,I+2) * T02 - C(I,I+2) * U02
                  E(I+2,3) = E(I+2,3) - A(I,I+2) * T03 - C(I,I+2) * U03
                  E(I+3,0) = E(I+3,0) - A(I,I+3) * T00 - C(I,I+3) * U00
                  E(I+3,1) = E(I+3,1) - A(I,I+3) * T01 - C(I,I+3) * U01
                  E(I+3,2) = E(I+3,2) - A(I,I+3) * T02 - C(I,I+3) * U02
                  E(I+3,3) = E(I+3,3) - A(I,I+3) * T03 - C(I,I+3) * U03
               END SELECT
            END SELECT
         ELSE
            IF (M-I-1.GT.0) THEN
               SELECT CASE (N)
               CASE (1)
                  T00 = E(I,0)   * B(0,0)
                  T10 = E(I+1,0) * B(0,0)
                  U00 = E(I,0)   * D(0,0)
                  U10 = E(I+1,0) * D(0,0)
               CASE (2)
                  T00 = E(I,0)   * B(0,0) + E(I,1)   * B(0,1)
                  T01 = E(I,0)   * B(1,0) + E(I,1)   * B(1,1)
                  T10 = E(I+1,0) * B(0,0) + E(I+1,1) * B(0,1)
                  T11 = E(I+1,0) * B(1,0) + E(I+1,1) * B(1,1)
                  U00 = E(I,0)   * D(0,0) + E(I,1)   * D(0,1)
                  U01 =                     E(I,1)   * D(1,1)
                  U10 = E(I+1,0) * D(0,0) + E(I+1,1) * D(0,1)
                  U11 =                     E(I+1,1) * D(1,1)
               CASE (3)
                  T00 = E(I,0)   * B(0,0) + E(I,1)   * B(0,1) + E(I,2)   * B(0,2)
                  T01 = E(I,0)   * B(1,0) + E(I,1)   * B(1,1) + E(I,2)   * B(1,2)
                  T02 =                     E(I,1)   * B(2,1) + E(I,2)   * B(2,2)
                  T10 = E(I+1,0) * B(0,0) + E(I+1,1) * B(0,1) + E(I+1,2) * B(0,2)
                  T11 = E(I+1,0) * B(1,0) + E(I+1,1) * B(1,1) + E(I+1,2) * B(1,2)
                  T12 =                     E(I+1,1) * B(2,1) + E(I+1,2) * B(2,2)
                  U00 = E(I,0)   * D(0,0) + E(I,1)   * D(0,1) + E(I,2)   * D(0,2)
                  U01 =                     E(I,1)   * D(1,1) + E(I,2)   * D(1,2)
                  U02 =                                         E(I,2)   * D(2,2)
                  U10 = E(I+1,0) * D(0,0) + E(I+1,1) * D(0,1) + E(I+1,2) * D(0,2)
                  U11 =                     E(I+1,1) * D(1,1) + E(I+1,2) * D(1,2)
                  U12 =                                         E(I+1,2) * D(2,2)
               CASE (4)
                  T00 = E(I,0)   * B(0,0) + E(I,1)   * B(0,1) + E(I,2)   * B(0,2) + E(I,3)   * B(0,3)
                  T01 = E(I,0)   * B(1,0) + E(I,1)   * B(1,1) + E(I,2)   * B(1,2) + E(I,3)   * B(1,3)
                  T02 =                     E(I,1)   * B(2,1) + E(I,2)   * B(2,2) + E(I,3)   * B(2,3)
                  T03 =                                         E(I,2)   * B(3,2) + E(I,3)   * B(3,3)
                  T10 = E(I+1,0) * B(0,0) + E(I+1,1) * B(0,1) + E(I+1,2) * B(0,2) + E(I+1,3) * B(0,3)
                  T11 = E(I+1,0) * B(1,0) + E(I+1,1) * B(1,1) + E(I+1,2) * B(1,2) + E(I+1,3) * B(1,3)
                  T12 =                     E(I+1,1) * B(2,1) + E(I+1,2) * B(2,2) + E(I+1,3) * B(2,3)
                  T13 =                                         E(I+1,2) * B(3,2) + E(I+1,3) * B(3,3)
                  U00 = E(I,0)   * D(0,0) + E(I,1)   * D(0,1) + E(I,2)   * D(0,2) + E(I,3)   * D(0,3)
                  U01 =                     E(I,1)   * D(1,1) + E(I,2)   * D(1,2) + E(I,3)   * D(1,3)
                  U02 =                                         E(I,2)   * D(2,2) + E(I,3)   * D(2,3)
                  U03 =                                                             E(I,3)   * D(3,3)
                  U10 = E(I+1,0) * D(0,0) + E(I+1,1) * D(0,1) + E(I+1,2) * D(0,2) + E(I+1,3) * D(0,3)
                  U11 =                     E(I+1,1) * D(1,1) + E(I+1,2) * D(1,2) + E(I+1,3) * D(1,3)
                  U12 =                                         E(I+1,2) * D(2,2) + E(I+1,3) * D(2,3)
                  U13 =                                                             E(I+1,3) * D(3,3)
               END SELECT
            END IF

            SELECT CASE (M-I-2)
            CASE (1)
               SELECT CASE (N)
               CASE (1)
                  E(I+2,0) = E(I+2,0) - A(I,I+2) * T00 - A(I+1,I+2) * T10 - C(I,I+2) * U00 - C(I+1,I+2) * U10
               CASE (2)                                          
                  E(I+2,0) = E(I+2,0) - A(I,I+2) * T00 - A(I+1,I+2) * T10 - C(I,I+2) * U00 - C(I+1,I+2) * U10
                  E(I+2,1) = E(I+2,1) - A(I,I+2) * T01 - A(I+1,I+2) * T11 - C(I,I+2) * U01 - C(I+1,I+2) * U11
               CASE (3)                                          
                  E(I+2,0) = E(I+2,0) - A(I,I+2) * T00 - A(I+1,I+2) * T10 - C(I,I+2) * U00 - C(I+1,I+2) * U10
                  E(I+2,1) = E(I+2,1) - A(I,I+2) * T01 - A(I+1,I+2) * T11 - C(I,I+2) * U01 - C(I+1,I+2) * U11
                  E(I+2,2) = E(I+2,2) - A(I,I+2) * T02 - A(I+1,I+2) * T12 - C(I,I+2) * U02 - C(I+1,I+2) * U12
               CASE (4)                                          
                  E(I+2,0) = E(I+2,0) - A(I,I+2) * T00 - A(I+1,I+2) * T10 - C(I,I+2) * U00 - C(I+1,I+2) * U10
                  E(I+2,1) = E(I+2,1) - A(I,I+2) * T01 - A(I+1,I+2) * T11 - C(I,I+2) * U01 - C(I+1,I+2) * U11
                  E(I+2,2) = E(I+2,2) - A(I,I+2) * T02 - A(I+1,I+2) * T12 - C(I,I+2) * U02 - C(I+1,I+2) * U12
                  E(I+2,3) = E(I+2,3) - A(I,I+2) * T03 - A(I+1,I+2) * T13 - C(I,I+2) * U03 - C(I+1,I+2) * U13
               END SELECT
            CASE (2)
               SELECT CASE (N)
               CASE (1)
                  E(I+2,0) = E(I+2,0) - A(I,I+2) * T00 - A(I+1,I+2) * T10 - C(I,I+2) * U00 - C(I+1,I+2) * U10
                  E(I+3,0) = E(I+3,0) - A(I,I+3) * T00 - A(I+1,I+3) * T10 - C(I,I+3) * U00 - C(I+1,I+3) * U10
               CASE (2)                                          
                  E(I+2,0) = E(I+2,0) - A(I,I+2) * T00 - A(I+1,I+2) * T10 - C(I,I+2) * U00 - C(I+1,I+2) * U10
                  E(I+2,1) = E(I+2,1) - A(I,I+2) * T01 - A(I+1,I+2) * T11 - C(I,I+2) * U01 - C(I+1,I+2) * U11
                  E(I+3,0) = E(I+3,0) - A(I,I+3) * T00 - A(I+1,I+3) * T10 - C(I,I+3) * U00 - C(I+1,I+3) * U10
                  E(I+3,1) = E(I+3,1) - A(I,I+3) * T01 - A(I+1,I+3) * T11 - C(I,I+3) * U01 - C(I+1,I+3) * U11
               CASE (3)                                          
                  E(I+2,0) = E(I+2,0) - A(I,I+2) * T00 - A(I+1,I+2) * T10 - C(I,I+2) * U00 - C(I+1,I+2) * U10
                  E(I+2,1) = E(I+2,1) - A(I,I+2) * T01 - A(I+1,I+2) * T11 - C(I,I+2) * U01 - C(I+1,I+2) * U11
                  E(I+2,2) = E(I+2,2) - A(I,I+2) * T02 - A(I+1,I+2) * T12 - C(I,I+2) * U02 - C(I+1,I+2) * U12
                  E(I+3,0) = E(I+3,0) - A(I,I+3) * T00 - A(I+1,I+3) * T10 - C(I,I+3) * U00 - C(I+1,I+3) * U10
                  E(I+3,1) = E(I+3,1) - A(I,I+3) * T01 - A(I+1,I+3) * T11 - C(I,I+3) * U01 - C(I+1,I+3) * U11
                  E(I+3,2) = E(I+3,2) - A(I,I+3) * T02 - A(I+1,I+3) * T12 - C(I,I+3) * U02 - C(I+1,I+3) * U12
               CASE (4)                                          
                  E(I+2,0) = E(I+2,0) - A(I,I+2) * T00 - A(I+1,I+2) * T10 - C(I,I+2) * U00 - C(I+1,I+2) * U10
                  E(I+2,1) = E(I+2,1) - A(I,I+2) * T01 - A(I+1,I+2) * T11 - C(I,I+2) * U01 - C(I+1,I+2) * U11
                  E(I+2,2) = E(I+2,2) - A(I,I+2) * T02 - A(I+1,I+2) * T12 - C(I,I+2) * U02 - C(I+1,I+2) * U12
                  E(I+2,3) = E(I+2,3) - A(I,I+2) * T03 - A(I+1,I+2) * T13 - C(I,I+2) * U03 - C(I+1,I+2) * U13
                  E(I+3,0) = E(I+3,0) - A(I,I+3) * T00 - A(I+1,I+3) * T10 - C(I,I+3) * U00 - C(I+1,I+3) * U10
                  E(I+3,1) = E(I+3,1) - A(I,I+3) * T01 - A(I+1,I+3) * T11 - C(I,I+3) * U01 - C(I+1,I+3) * U11
                  E(I+3,2) = E(I+3,2) - A(I,I+3) * T02 - A(I+1,I+3) * T12 - C(I,I+3) * U02 - C(I+1,I+3) * U12
                  E(I+3,3) = E(I+3,3) - A(I,I+3) * T03 - A(I+1,I+3) * T13 - C(I,I+3) * U03 - C(I+1,I+3) * U13
               END SELECT
            END SELECT
         END IF
         I = I+II
      END DO
      RETURN
      END SUBROUTINE
