!  This file is part of the RECSY library. See http://www.cs.umu.se/~isak/recsy
!  Copyright 2003 Umeå University, Sweden
!
      SUBROUTINE RECSY_MACHINE(MACHINE)

      IMPLICIT NONE

      DOUBLE PRECISION MACHINE(0:*)
      DOUBLE PRECISION DLAMCH

      MACHINE(2) = DLAMCH( 'P' )
      MACHINE(0) = DLAMCH( 'S' ) / MACHINE(2)
      MACHINE(1) = 1D0 / MACHINE(0)
      MACHINE(3) = 40
      MACHINE(4) = 40
      MACHINE(5) = 3E6

      END SUBROUTINE RECSY_MACHINE
