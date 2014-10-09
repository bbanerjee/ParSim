!******************************************************************************!
!      Written by Kumar Mithraratne                                            !
!      University of Auckland                                                  !
!                                                                              !
!      01/04/2014                                                              !
!*******************************************************************************


	SUBROUTINE dywrtdxi(Jb,dy_dxi,xi1,xi2,xi3,YP)  
	
	
	IMPLICIT NONE
	
	REAL*8 Jb,dy_dxi(3,3),xi1,xi2,xi3,YP(3,8)
	REAL*8 df_dxi(8,3)
	INTEGER i,j,n
	
	DO i=1,3,1
	  DO j=1,3,1
	    dy_dxi(i,j)=0.0d0
          ENDDO
	ENDDO
	
	df_dxi(1,1)=-(1.0d0-xi2)*(1.0d0-xi3)
	df_dxi(1,2)=-(1.0d0-xi1)*(1.0d0-xi3)
	df_dxi(1,3)=-(1.0d0-xi1)*(1.0d0-xi2)	
	
	df_dxi(2,1)=(1.0d0-xi2)*(1.0d0-xi3)
	df_dxi(2,2)=-(xi1)*(1.0d0-xi3)
	df_dxi(2,3)=-(xi1)*(1.0d0-xi2)
		
	df_dxi(3,1)=-(xi2)*(1.0d0-xi3)
	df_dxi(3,2)=(1.0d0-xi1)*(1.0d0-xi3)
	df_dxi(3,3)=-(1.0d0-xi1)*(xi2)	
	
	df_dxi(4,1)=(xi2)*(1.0d0-xi3)
	df_dxi(4,2)=(xi1)*(1.0d0-xi3)
	df_dxi(4,3)=-(xi1)*(xi2)	
	
	df_dxi(5,1)=-(1.0d0-xi2)*(xi3)
	df_dxi(5,2)=-(1.0d0-xi1)*(xi3)
	df_dxi(5,3)=(1.0d0-xi1)*(1.0d0-xi2)
	
	df_dxi(6,1)=(1.0d0-xi2)*(xi3)
	df_dxi(6,2)=-(xi1)*(xi3)
	df_dxi(6,3)=(xi1)*(1.0d0-xi2)
		
	df_dxi(7,1)=-(xi2)*(xi3)
	df_dxi(7,2)=(1.0d0-xi1)*(xi3)
	df_dxi(7,3)=(1.0d0-xi1)*(xi2)
	
	df_dxi(8,1)=(xi2)*(xi3)
	df_dxi(8,2)=(xi1)*(xi3)
	df_dxi(8,3)=(xi1)*(xi2)	
	
		
	DO j=1,3,1
          DO i=1,3,1
            dy_dxi(j,i)=0.0d0
            DO n=1,8,1                       
	        dy_dxi(j,i)=dy_dxi(j,i)+YP(j,n)*df_dxi(n,i)               
	     ENDDO
          ENDDO    
	ENDDO      
	  
       Jb=dy_dxi(1,1)*(dy_dxi(2,2)*dy_dxi(3,3)-dy_dxi(2,3)*dy_dxi(3,2))+
     1    dy_dxi(1,2)*(dy_dxi(2,3)*dy_dxi(3,1)-dy_dxi(2,1)*dy_dxi(3,3))+
     2    dy_dxi(1,3)*(dy_dxi(2,1)*dy_dxi(3,2)-dy_dxi(2,2)*dy_dxi(3,1))

!	DO i=1,3,1
!	   WRITE(6,100) (dy_dxi(i,j),j=1,3)
!	ENDDO
!100    FORMAT (24(1x,f7.5))
!	print*,'*********************************'
	
	RETURN
	END
