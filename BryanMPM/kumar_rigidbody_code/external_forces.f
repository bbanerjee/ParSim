	SUBROUTINE external_forces(eforce,etorque)
	
	
!***********************************************************************
!	Written by Kumar Mithraratne
!	© Auckland Bioengineering Institute
!	V1.1 July 2013
!***********************************************************************
	
	IMPLICIT NONE
	
	INCLUDE 'rigidbody1.cmn'
	
	
	REAL*8 eforce(nbm,3),etorque(nbm,3)
	!INTEGER
	
	!REAL*8 
	INTEGER nb,nj

        
	DO nb=1,nbm,1
	   DO nj=1,3,1
	      eforce(nb,nj)=0.0d0
	      etorque(nb,nj)=0.0d0
	   ENDDO
	ENDDO 


	RETURN
	END


