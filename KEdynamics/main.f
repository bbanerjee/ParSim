!******************************************************************************!
!      Written by Kumar Mithraratne                                            !
!      University of Auckland                                                  !
!                                                                              !
!      01/04/2014                                                              !
!*******************************************************************************

	PROGRAM main
	
	IMPLICIT NONE

	REAL*8 XP(3,8),UP(2,4,0:2),C(6,6),roh,FN(3,8),UN(3,8,0:1),DUNDT(3,8,0:1)
	REAL*8 M(24,24),K(24,24)
	INTEGER itt,converged,FIX(25)


	CALL undeformed(XP)

       CALL material(C,roh)

       CALL initial_boundary(FN,UN,DUNDT)

   	CALL mass_matrix(M,XP,roh)	
	
	CALL stiffness_matrix(K,XP,C)










	STOP
	END
