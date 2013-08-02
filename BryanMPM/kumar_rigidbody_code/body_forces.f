	SUBROUTINE body_forces(bforce)
	
	
!***********************************************************************
!	Written by Kumar Mithraratne
!	© Auckland Bioengineering Institute
!	V1.1 July 2013
!***********************************************************************
	
	IMPLICIT NONE
	
	INCLUDE 'rigidbody1.cmn'
	
	
	REAL*8 bforce(nbm,3)
	!INTEGER
	
	!REAL*8 
	INTEGER nb

        
	DO nb=1,nbm,1
           bforce(nb,1)=0.0d0
	   bforce(nb,2)=0.0d0
	   bforce(nb,3)=-g
	ENDDO 


	RETURN
	END


