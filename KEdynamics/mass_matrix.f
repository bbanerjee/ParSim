!******************************************************************************!
!      Written by Kumar Mithraratne                                            !
!      University of Auckland                                                  !
!                                                                              !
!      01/04/2014                                                              !
!*******************************************************************************


	SUBROUTINE mass_matrix(MASS,XP,roh)
	
	
	IMPLICIT NONE
	
	INTEGER i,j,g1,g2,g3,NGT
	REAL*8 MASS(24,24),XP(3,8),roh
	REAL*8 XI(3,3),w(3),xi1,xi2,xi3,Jxxi,fi(8),dx_dxi(3,3)	
        
	 
	REAL*8 FIMATRIX(3,24),FIMATRIXT(24,3),LOCALMASS(24,24) 

!******Matrix Multiplication Example********************************************
	INTEGER M,N,K,LDA,LDB,LDC
	REAL*8 ALPHA,BETA
	CHARACTER*1 TRANSA,TRANSB
       REAL*8 A(3,2),B(2,1),C(3,1)
	
 
       A(1,1)=1.0d0; A(1,2)=3.0d0
	A(2,1)=4.0d0; A(2,2)=2.0d0
	A(3,1)=2.0d0; A(3,2)=1.0d0
	
	B(1,1)=2.0
	B(2,1)=1.0
	
       C(1,1)=0.0d0
	C(2,1)=0.0d0
	C(3,1)=0.0d0
	
	TRANSA='N'  !no transpose, transpose or conjugate transpose
	TRANSB='N'  !no transpose, transpose or conjugate transpose
	M=3 	     ! number of rows of A
	N=1         ! number of columns of B
	K=2         ! no. of columns of A or no. of rows of B (should be equal)
	ALPHA=1.0d0 ! scalar multiplier of A*B
	LDA=3       ! first dimension of A (rows)
	LDB=2       ! first dimension of B (rows)
	BETA=0.0d0  ! scalar multiplier of C
	LDC=3       ! first dimension of C (rows)
	
	!C = ALPHA*A*B + BETA*C

       CALL dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
 !******Matrix Multiplication Example********************************************

	 
!*******************************************************************************
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
	
!	DO i=1,8,1
!	   fi(i)=0.0d0
!	ENDDO
	
!	DO i=1,24,1
!	   DO j=1,24,1
!	      MASS(i,j)=0.0d0
!	   ENDDO
!	ENDDO

!       NGT=0
!	DO g2=1,3,1
!	   DO g1=1,3,1
!	      DO g3=1,3,1
	    
!	         NGT=NGT+1
	       
!	         xi1=XI(1,g1)
!	         xi2=XI(2,g2)
!		  xi3=XI(3,g3)
	      
!                CALL dywrtdxi(Jxxi,dx_dxi,xi1,xi2,xi3,XP) 

!	         fi(1)=(1.0d0-xi1)*(1.0d0-xi2)*(1.0d0-xi3)
!	         fi(2)=(xi1)*(1.0d0-xi2)*(1.0d0-xi3)
!	         fi(3)=(1.0d0-xi1)*(xi2)*(1.0d0-xi3)
!	         fi(4)=(xi1)*(xi2)*(1.0d0-xi3)
!	         fi(5)=(1.0d0-xi1)*(1.0d0-xi2)*(xi3)
!	         fi(6)=(xi1)*(1.0d0-xi2)*(xi3)
!	         fi(7)=(1.0d0-xi1)*(xi2)*(xi3)
!	         fi(8)=(xi1)*(xi2)*(xi3)
	      	      		      
!		  MASS(1,1)=MASS(1,1)+Jxxi*w(g1)*w(g2)*w(g3)*fi(1)*fi(1)
!	         MASS(1,2)=MASS(1,2)+Jxxi*w(g1)*w(g2)*w(g3)*fi(1)*fi(2)
!	         MASS(1,3)=MASS(1,3)+Jxxi*w(g1)*w(g2)*w(g3)*fi(1)*fi(3)
!	         MASS(1,4)=MASS(1,4)+Jxxi*w(g1)*w(g2)*w(g3)*fi(1)*fi(4)
!		  MASS(1,5)=MASS(1,5)+Jxxi*w(g1)*w(g2)*w(g3)*fi(1)*fi(5)
!	         MASS(1,6)=MASS(1,6)+Jxxi*w(g1)*w(g2)*w(g3)*fi(1)*fi(6)
!	         MASS(1,7)=MASS(1,7)+Jxxi*w(g1)*w(g2)*w(g3)*fi(1)*fi(7)
!	         MASS(1,8)=MASS(1,8)+Jxxi*w(g1)*w(g2)*w(g3)*fi(1)*fi(8)
	     
!	         MASS(2,2)=MASS(2,2)+Jxxi*w(g1)*w(g2)*w(g3)*fi(2)*fi(2)
!	         MASS(2,3)=MASS(2,3)+Jxxi*w(g1)*w(g2)*w(g3)*fi(2)*fi(3)
!	         MASS(2,4)=MASS(2,4)+Jxxi*w(g1)*w(g2)*w(g3)*fi(2)*fi(4)
!	         MASS(2,5)=MASS(2,5)+Jxxi*w(g1)*w(g2)*w(g3)*fi(2)*fi(5)
!	         MASS(2,6)=MASS(2,6)+Jxxi*w(g1)*w(g2)*w(g3)*fi(2)*fi(6)
!	         MASS(2,7)=MASS(2,7)+Jxxi*w(g1)*w(g2)*w(g3)*fi(2)*fi(7)
!	         MASS(2,8)=MASS(2,8)+Jxxi*w(g1)*w(g2)*w(g3)*fi(2)*fi(8)
		  	         
!	         MASS(3,3)=MASS(3,3)+Jxxi*w(g1)*w(g2)*w(g3)*fi(3)*fi(3)
!	         MASS(3,4)=MASS(3,4)+Jxxi*w(g1)*w(g2)*w(g3)*fi(3)*fi(4)
!	         MASS(3,5)=MASS(3,5)+Jxxi*w(g1)*w(g2)*w(g3)*fi(3)*fi(5)
!	         MASS(3,6)=MASS(3,6)+Jxxi*w(g1)*w(g2)*w(g3)*fi(3)*fi(6)		  
!	         MASS(3,7)=MASS(3,7)+Jxxi*w(g1)*w(g2)*w(g3)*fi(3)*fi(7)
!	         MASS(3,8)=MASS(3,8)+Jxxi*w(g1)*w(g2)*w(g3)*fi(3)*fi(8)
		  	
!	         MASS(4,4)=MASS(4,4)+Jxxi*w(g1)*w(g2)*w(g3)*fi(4)*fi(4)
!	         MASS(4,5)=MASS(4,5)+Jxxi*w(g1)*w(g2)*w(g3)*fi(4)*fi(5)
!	         MASS(4,6)=MASS(4,6)+Jxxi*w(g1)*w(g2)*w(g3)*fi(4)*fi(6)
!	         MASS(4,7)=MASS(4,7)+Jxxi*w(g1)*w(g2)*w(g3)*fi(4)*fi(7)
!	         MASS(4,8)=MASS(4,8)+Jxxi*w(g1)*w(g2)*w(g3)*fi(4)*fi(8)
		  
!	         MASS(5,5)=MASS(5,5)+Jxxi*w(g1)*w(g2)*w(g3)*fi(5)*fi(5)
!	         MASS(5,6)=MASS(5,6)+Jxxi*w(g1)*w(g2)*w(g3)*fi(5)*fi(6)
!	         MASS(5,7)=MASS(5,7)+Jxxi*w(g1)*w(g2)*w(g3)*fi(5)*fi(7)
!	         MASS(5,8)=MASS(5,8)+Jxxi*w(g1)*w(g2)*w(g3)*fi(5)*fi(8)
		  
!	         MASS(6,6)=MASS(6,6)+Jxxi*w(g1)*w(g2)*w(g3)*fi(6)*fi(6)
!	         MASS(6,7)=MASS(6,7)+Jxxi*w(g1)*w(g2)*w(g3)*fi(6)*fi(7)
!	         MASS(6,8)=MASS(6,8)+Jxxi*w(g1)*w(g2)*w(g3)*fi(6)*fi(8)
	         
!		  MASS(7,7)=MASS(7,7)+Jxxi*w(g1)*w(g2)*w(g3)*fi(7)*fi(7)
!	         MASS(7,8)=MASS(7,8)+Jxxi*w(g1)*w(g2)*w(g3)*fi(7)*fi(8)
		  
!	         MASS(8,8)=MASS(8,8)+Jxxi*w(g1)*w(g2)*w(g3)*fi(8)*fi(8)
		  
!	      ENDDO
!	   ENDDO  		
!	ENDDO
	 
!       DO i=2,8,1
!	   MASS(i,1)=MASS(1,i)
!	ENDDO
!       DO i=3,8,1
!	   MASS(i,2)=MASS(2,i)
!	ENDDO
!       DO i=4,8,1
!	   MASS(i,3)=MASS(3,i)
!	ENDDO
!       DO i=5,8,1
!	   MASS(i,4)=MASS(4,i)
!	ENDDO
!       DO i=6,8,1
!	   MASS(i,5)=MASS(5,i)
!	ENDDO
!       DO i=7,8,1
!	   MASS(i,6)=MASS(6,i)
!	ENDDO			
!       MASS(8,7)=MASS(7,8)
		
!	DO i=9,16,1
!	   DO j=9,16,1
!	      MASS(i,j)=MASS(i-8,j-8)
!	   ENDDO
!	ENDDO
	
!	DO i=17,24,1
!	   DO j=17,24,1
!	      MASS(i,j)=MASS(i-8,j-8)
!	   ENDDO
!	ENDDO

 !      DO i=1,24,1
!	   DO j=1,24,1
!	      MASS(i,j)=roh*MASS(i,j)
!	   ENDDO
!	ENDDO	
	
!	DO i=1,24,1
!	   WRITE(6,100) (MASS(i,j),j=1,24)
!	ENDDO
!100    FORMAT (24(1x,f7.3))
!       WRITE(6,'("************************************************")')
!*******************************************************************************
       DO i=1,8,1
	   fi(i)=0.0d0
	ENDDO

       DO i=1,3,1
	   DO j=1,24,1
	      FIMATRIX(i,j)=0.0d0
	      FIMATRIXT(j,i)=0.0d0
	   ENDDO
	ENDDO
       DO i=1,24,1
	   DO j=1,24,1
	      LOCALMASS(i,j)=0.0d0
	      MASS(i,j)=0.0d0
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

	         fi(1)=(1.0d0-xi1)*(1.0d0-xi2)*(1.0d0-xi3)
	         fi(2)=(xi1)*(1.0d0-xi2)*(1.0d0-xi3)
	         fi(3)=(1.0d0-xi1)*(xi2)*(1.0d0-xi3)
	         fi(4)=(xi1)*(xi2)*(1.0d0-xi3)
	         fi(5)=(1.0d0-xi1)*(1.0d0-xi2)*(xi3)
	         fi(6)=(xi1)*(1.0d0-xi2)*(xi3)
	         fi(7)=(1.0d0-xi1)*(xi2)*(xi3)
	         fi(8)=(xi1)*(xi2)*(xi3)

                DO j=1,8,1
		     FIMATRIX(1,j)=fi(j)    
		     FIMATRIXT(j,1)=fi(j)
		     FIMATRIX(2,8+j)=fi(j)  
		     FIMATRIXT(8+j,2)=fi(j)
		     FIMATRIX(3,16+j)=fi(j) 
		     FIMATRIXT(16+j,3)=fi(j)
		  ENDDO

                CALL dgemm('N','N',24,24,3,roh,FIMATRIXT,24,FIMATRIX,
	1                   3,0.0d0,LOCALMASS,24)
		   
                DO i=1,24,1
		     DO j=1,24,1
		        MASS(i,j)=MASS(i,j)+
	1                        Jxxi*w(g1)*w(g2)*w(g3)*LOCALMASS(i,j)
	            ENDDO
		  ENDDO   
	 

	      ENDDO
	   ENDDO  		
	ENDDO

!	DO i=1,24,1
!	   WRITE(6,100) (MASS(i,j),j=1,24)
!	ENDDO
!100    FORMAT (24(1x,f7.3))

	
	RETURN
	END
