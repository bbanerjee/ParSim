	SUBROUTINE initial_conditions(posni,lveli,lacci,
     '             ornti,aveli,aacci)	
	
	
!***********************************************************************
!	Written by Kumar Mithraratne
!	© Auckland Bioengineering Institute
!	V1.1 July 2013
!***********************************************************************
	
	IMPLICIT NONE
	
	INCLUDE 'rigidbody1.cmn'
	
	
	REAL*8 posni(nbm,3),lveli(nbm,3),lacci(nbm,3),
     '         ornti(nbm,3,3),aveli(nbm,3),aacci(nbm,3)
	!INTEGER
	
	!REAL*8 
	INTEGER n,nb,nj

         
	!initialise
	DO nb=1,nbm,1
	   DO nj=1,3,1
	      posni(nb,nj)=0.0d0
	      lveli(nb,nj)=0.0d0
	      lacci(nb,nj)=0.0d0
	      
	      aveli(nb,nj)=0.0d0
	      aacci(nb,nj)=0.0d0
	      DO n=1,3,1
	         ornti(nb,nj,n)=0.0d0
	      ENDDO
           ENDDO
	ENDDO

	!!!!!!!This information must be read in from an input file - nb=1..nbt
	DO nb=1,nbt,1
           posni(nb,1)=200.0d0/mm2m  !mm
	   posni(nb,2)=250.0d0/mm2m  !mm
	   posni(nb,3)=375.0d0/mm2m  !mm
        
           lveli(nb,1)=2000.0d0/mm2m !mm/s
           lveli(nb,2)=3000.0d0/mm2m !mm/s
           lveli(nb,3)=5000.0d0/mm2m !mm/s
	
	   lacci(nb,1)=0.0d0/mm2m    !mm/s2
	   lacci(nb,2)=0.0d0/mm2m    !mm/s2
	   lacci(nb,3)=0.0d0/mm2m    !mm/s2
			
	   ornti(nb,1,1)=1.0d0  !normalised vector 1
	   ornti(nb,2,1)=0.0d0	
	   ornti(nb,3,1)=0.0d0	
	   ornti(nb,1,2)=0.0d0  !normalised vector 2
	   ornti(nb,2,2)=1.0d0	
	   ornti(nb,3,2)=0.0d0	
	   ornti(nb,1,3)=0.0d0  !normalised vector 3
	   ornti(nb,2,3)=0.0d0	
	   ornti(nb,3,3)=1.0d0	
	
           aveli(nb,1)=20.0d0    !deg/s
           aveli(nb,2)=0.0d0    !deg/s
           aveli(nb,3)=10.0d0   !deg/s
		
           aacci(nb,1)=0.0d0    !deg/s2
           aacci(nb,2)=0.0d0    !deg/s2
           aacci(nb,3)=0.0d0    !deg/s2
	ENDDO   
		

	RETURN
	END


