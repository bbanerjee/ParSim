!******************************************************************************!
!      Written by Kumar Mithraratne                                            !
!      University of Auckland                                                  !
!                                                                              !
!      01/04/2014                                                              !
!*******************************************************************************


	SUBROUTINE stiffness_matrix(STIFFNESS,XP,C)
	
	
	IMPLICIT NONE
	
	INTEGER i,j,g1,g2,g3,NGT
	INTEGER INFO,IPIV(3)
	REAL*8 STIFFNESS(24,24),XP(3,8),C(6,6)
	REAL*8 XI(3,3),w(3),xi1,xi2,xi3,Jxxi,df_dxi(8,3),dx_dxi(3,3)	
	REAL*8 B(6,24),BT(24,6),temp(24,6)
	REAL*8 LOCALSTIFFNESS(24,24),dxi_dx(3,3),df_dx(8,3)

	XI(1,1)=0.112701665d0
	XI(1,2)=0.500000000d0
       XI(1,3)=0.887298334d0		
       XI(2,1)=0.112701665d0
       XI(2,2)=0.500000000d0
       XI(2,3)=0.887298334d0	
       XI(3,1)=0.112701665d0
       XI(3,2)=0.500000000d0
       XI(3,3)=0.887298334d0	

       w(1)=0.277778d0
       w(2)=0.444444d0	
       w(3)=0.277778d0	


       DO i=1,8,1
	   DO j=1,3,1
	      df_dxi(i,j)=0.0d0
	   ENDDO
	ENDDO

       DO i=1,6,1
	   DO j=1,24,1
	      B(i,j)=0.0d0
	      BT(j,i)=0.0d0
	      temp(j,i)=0.0d0
	   ENDDO
	ENDDO
      
       DO i=1,24,1
	   DO j=1,24,1
	      LOCALSTIFFNESS(i,j)=0.0d0
	      STIFFNESS(i,j)=0.0d0
	   ENDDO
	ENDDO

       DO i=1,3,1
	   DO j=1,3,1
	      dxi_dx(i,j)=0.0d0
	      IF (i.EQ.j) dxi_dx(i,j)=1.0d0
          ENDDO
	ENDDO
	   
       NGT=0
	DO g2=1,3,1
	   DO g1=1,3,1
	      DO g3=1,3,1
	    
	         NGT=NGT+1
	       
	         xi1=XI(1,g1)
	         xi2=XI(2,g2)
		  xi3=XI(3,g3)
	      
                CALL dywrtdxi(Jxxi,dx_dxi,xi1,xi2,xi3,XP) 

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
		  
		  !test
		  !dx_dxi(1,1)=1.0d0; dx_dxi(1,2)=2.0d0; dx_dxi(1,3)=1.0d0; 
		  !dx_dxi(2,1)=2.0d0; dx_dxi(2,2)=2.0d0; dx_dxi(2,3)=3.0d0; 
		  !dx_dxi(3,1)=1.0d0; dx_dxi(3,2)=3.0d0; dx_dxi(3,3)=1.0d0; 

                CALL DGETRF(3,3,dx_dxi,3,IPIV,INFO)
                CALL DGETRS('N',3,3,dx_dxi,3,IPIV,dxi_dx,3,INFO) 

!	         DO i=1,3,1
!	            WRITE(6,100) (dx_dxi(i,j),j=1,3)
!	         ENDDO
!	         DO i=1,3,1
!	            WRITE(6,100) (dxi_dx(i,j),j=1,3)
!	         ENDDO		  
!100             FORMAT (3(1x,f7.3))

                CALL dgemm('N','N',8,3,3,1.0d0,df_dxi,8,dxi_dx,
	1                   3,0.0d0,df_dx,8)
               
		  DO j=1,8,1
		     B(1,j)=df_dx(j,1)
                   B(2,8+j)=df_dx(j,2)
                   B(3,16+j)=df_dx(j,3)
                   
		     B(4,j)=df_dx(j,2)
                   B(4,8+j)=df_dx(j,1)
                   
		     B(5,8+j)=df_dx(j,3)
		     B(5,16+j)=df_dx(j,2)
		     
		     B(6,j)=df_dx(j,3)
		     B(6,16+j)=df_dx(j,1)
		  ENDDO
		  
!	         DO i=1,6,1
!	            WRITE(6,80) (B(i,j),j=1,24)
!	         ENDDO  
!80              FORMAT (24(1x,f7.3))

!	         DO i=1,6,1
!	            WRITE(6,85) (C(i,j),j=1,6)
!	         ENDDO  
!85              FORMAT (6(1x,f15.2))

!		  STOP		  
		  
		  DO i=1,6,1
		     DO j=1,24,1
		        BT(j,i)=B(i,j)
		     ENDDO
		  ENDDO
		  
                CALL dgemm('N','N',24,6,6,1.0d0,BT,24,C,
	1                   6,0.0d0,temp,24)
!	         DO i=1,6,1
!	            WRITE(6,90) (temp(i,j),j=1,24)
!	         ENDDO  
!90              FORMAT (24(1x,f7.3))
!		  STOP		  
	
                CALL dgemm('N','N',24,24,6,1.0d0,temp,24,B,
	1                   6,0.0d0,LOCALSTIFFNESS,24)

!	         DO i=1,24,1
!	            WRITE(6,90) (LOCALSTIFFNESS(i,j)*1.0d-05,j=1,24)
!	         ENDDO  
!90              FORMAT (24(1x,f7.4))
!		  STOP		  
			  
                DO i=1,24,1
		     DO j=1,24,1
		        STIFFNESS(i,j)=STIFFNESS(i,j)+
	1                 Jxxi*w(g1)*w(g2)*w(g3)*LOCALSTIFFNESS(i,j)
	            ENDDO
		  ENDDO   


	      ENDDO
	   ENDDO  		
	ENDDO

!	DO i=1,24,1
!	   WRITE(6,100) (STIFFNESS(i,j)*1.0d-05,j=1,24)
!	ENDDO  
!100    FORMAT (24(1x,f7.4))


	RETURN
	END
