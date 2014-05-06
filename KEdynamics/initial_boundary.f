!******************************************************************************!
!      Written by Kumar Mithraratne                                            !
!      University of Auckland                                                  !
!                                                                              !
!      01/04/2014                                                              !
!*******************************************************************************


	SUBROUTINE initial_boundary(FN,UN,DUNDT)
	
	
	IMPLICIT NONE
	
	INTEGER nj,nn,nt
	REAL*8 FN(3,8),UN(3,8,0:1),DUNDT(3,8,0:1)	
	
	DO nj=1,3,1
	   DO nn=1,8,1
	      FN(nj,nn)=0.0d0
	      DO nt=0,1,1
	         UN(nj,nn,nt)=0.0d0
	         DUNDT(nj,nn,nt)=0.0d0
	      ENDDO
	   ENDDO
	ENDDO
	      


	
	RETURN
	END
