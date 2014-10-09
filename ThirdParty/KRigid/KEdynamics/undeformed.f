!******************************************************************************!
!      Written by Kumar Mithraratne                                            !
!      University of Auckland                                                  !
!                                                                              !
!      01/04/2014                                                              !
!*******************************************************************************

	SUBROUTINE undeformed(XP)
	
	IMPLICIT NONE
	
	INTEGER nj,nn
	REAL*8 XP(3,8)	
	
	DO nj=1,3,1
	  DO nn=1,8,1
	    XP(nj,nn)=0.0d0
	  ENDDO
	ENDDO
	
	XP(1,1)=0.0d0
	XP(1,2)=1.0d0
	XP(1,3)=0.0d0
	XP(1,4)=1.0d0
	XP(1,5)=0.0d0
	XP(1,6)=1.0d0
	XP(1,7)=0.0d0
	XP(1,8)=1.0d0
		
	XP(2,1)=0.0d0
	XP(2,2)=0.0d0
	XP(2,3)=1.0d0
	XP(2,4)=1.0d0	
	XP(2,5)=0.0d0
	XP(2,6)=0.0d0
	XP(2,7)=1.0d0
	XP(2,8)=1.0d0	

	XP(3,1)=0.0d0
	XP(3,2)=0.0d0
	XP(3,3)=0.0d0
	XP(3,4)=0.0d0	
	XP(3,5)=1.0d0
	XP(3,6)=1.0d0
	XP(3,7)=1.0d0
	XP(3,8)=1.0d0	
	
	RETURN
	END
