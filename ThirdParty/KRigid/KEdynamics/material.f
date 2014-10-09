!******************************************************************************!
!      Written by Kumar Mithraratne                                            !
!      University of Auckland                                                  !
!                                                                              !
!      01/04/2014                                                              !
!*******************************************************************************


	SUBROUTINE material(C,roh)
	
		
	IMPLICIT NONE
	
	REAL*8 C(6,6),roh,poisson,young,factor
	INTEGER i,j
	
	roh=1500.0d0             ! density (kg/m3)
	poisson=0.47d0           ! poisson ratio
	young=0.001d0*1.0d9     ! Young's modulus (Pa)
	
	DO i=1,6,1
	  DO j=1,6,1
	     C(i,j)=0.0d0
	  ENDDO
	ENDDO
	 	 
	factor=young/((1.0d0+poisson)*(1.0d0-poisson))
	
	C(1,1)=factor*(1.0d0-poisson)
	C(1,2)=factor*poisson
	C(1,3)=factor*poisson
	
	C(2,1)=factor*poisson
	C(2,2)=factor*(1.0d0-poisson)
	C(2,3)=factor*poisson
	
	C(3,1)=factor*poisson
	C(3,2)=factor*poisson
	C(3,3)=factor*(1.0d0-poisson)
	
	C(4,4)=factor*(1.0d0-2.0d0*poisson)/2.0d0
	C(5,5)=factor*(1.0d0-2.0d0*poisson)/2.0d0
	C(6,6)=factor*(1.0d0-2.0d0*poisson)/2.0d0
	
	roh=1500.0d0

!	DO i=1,6,1
!	   WRITE(6,85) (C(i,j),j=1,6)
!	ENDDO  
!85     FORMAT (6(1x,f10.3))
!       WRITE(6,'(2x,f15.2)') factor
	
	RETURN
	END
