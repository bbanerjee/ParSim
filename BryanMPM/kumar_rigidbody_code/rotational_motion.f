	SUBROUTINE rotational_motion(volume,density,rog,
     '       ornti,aveli,aacci,etorque,
     '       time_initial,time_final,dt,output_interval,rotation)
	
	
!***********************************************************************
!	Written by Kumar Mithraratne
!	© Auckland Bioengineering Institute
!	V1.1 July 2013
!***********************************************************************
	
	IMPLICIT NONE
	
	INCLUDE 'rigidbody1.cmn'
		
	REAL*8 volume(nbm),density(nbm),rog(nbm,3),
     '         ornti(nbm,3,3),aveli(nbm,3),aacci(nbm,3),
     '         etorque(nbm,3),
     '         time_initial,time_final,dt,output_interval,
     '         rotation(nbm,0:time_steps,18)	
	!INTEGER
	
	REAL*8 mass(nbm),
     '         ornt(0:1,9),theta(0:1,3),avel(0:1,3),aacc(0:1,3),
     '         time,out_time,iten(nbm,3,3),iiten(nbm,3,3),
     '         R(3,3),R0(3,3),Rx(3,3),Ry(3,3),Rz(3,3),
     '         C1(3,3),C2(3,3)
	INTEGER nb,nt,nv,i,j,nc,nc1,nc2,nj,comp_time_steps,ndt,k,
     '	        time_counter
     
	!initialise arrays
	DO nb=1,nbm,1
	   DO nt=0,time_steps,1
	      DO nv=1,15,1
	         rotation(nb,nt,nv)=0.0d0
	      ENDDO
	   ENDDO
	   DO i=1,3,1
	      DO j=1,3,1
	         iten(nb,i,j)=0.0d0
	      ENDDO
	   ENDDO
	ENDDO
	DO k=0,1,1
	   nc=0
	   DO nj=1,3,1
	      theta(k,nj)=0.0d0
	      avel(k,nj)=0.0d0
	      aacc(k,nj)=0.0d0
	      DO nv=1,3,1
	         nc=nc+1
		 ornt(k,nc)=0.0d0
	      ENDDO
	   ENDDO
	ENDDO
		
	!populate arrays with initial conditions
	DO nb=1,nbt,1
	   mass(nb)=volume(nb)*density(nb)
	   DO i=1,3,1
	      DO j=1,3,1
	         IF (i.EQ.j) THEN
	             iten(nb,i,j)=mass(nb)*(rog(nb,i))**2.0d0
		     iiten(nb,i,j)=1.0d0/iten(nb,i,j)
		 ENDIF    
	      ENDDO
	   ENDDO

           !update rotation array for time=0
	   nc1=0
	   DO nj=1,3,1
	      DO nv=1,3,1
	         nc1=nc1+1
	         rotation(nb,0,nc1)=ornti(nb,nj,nv)  !1 -9 orientaion of principal axes	
	      ENDDO
	      nc2=nj+9
	      rotation(nb,0,nc2)=0.0d0               !10 - 12 theta
	      nc2=nj+12
	      rotation(nb,0,nc2)=aveli(nb,nj)	     !13 - 15 angular velocity						
	      nc2=nj+15
	      rotation(nb,0,nc2)=aacci(nb,nj)        !15 - 18 angular accelaration
	   ENDDO	   
	ENDDO
 
        comp_time_steps=int((time_final-time_initial)/dt+dt/2.0d0)

	!DO nb=1,nbt,1
	DO nb=1,1,1
	   time=0.0d0   
	   out_time=0.0d0
	   time_counter=0            	   	   
	   nc=0
	   DO nj=1,3,1
	      DO nv=1,3,1
	         nc=nc+1
		 ornt(0,nc)=ornti(nb,nj,nv)
	      ENDDO	 
	      avel(0,nj)=aveli(nb,nj)
	      aacc(0,nj)=aacci(nb,nj)
	      theta(0,nj)=0.0d0	      
	   ENDDO   
	   
	   DO ndt=1,comp_time_steps,1

	      IF(abs(time-out_time).LT.dt/2.0d0) THEN  !update rotation array	
		 IF (time.NE.0.0d0) THEN
	             nc1=0
	             DO nj=1,3,1
	                DO nv=1,3,1
	                   nc1=nc1+1
	                   rotation(nb,time_counter,nc1)=
     '                     ornt(1,nc1)	
	                ENDDO
	                nc2=nj+9
	                rotation(nb,time_counter,nc2)=theta(1,nj)
	                nc2=nj+12
	                rotation(nb,time_counter,nc2)=avel(1,nj)
	                nc2=nj+15
	                rotation(nb,time_counter,nc2)=aacc(1,nj)			
	             ENDDO
		 ENDIF
		 time_counter=time_counter+1		 
	         out_time=out_time+output_interval
	      ENDIF
                
              avel(1,1)=avel(0,1)+(dt*avel(0,2)*avel(0,3))*
     '                  (iten(nb,2,2)-iten(nb,3,3))/(iten(nb,1,1))
              avel(1,2)=avel(0,2)+(dt*avel(0,3)*avel(0,1))*
     '                  (iten(nb,3,3)-iten(nb,1,1))/(iten(nb,2,2))
              avel(1,3)=avel(0,3)+(dt*avel(0,1)*avel(0,2))*
     '                  (iten(nb,1,1)-iten(nb,2,2))/(iten(nb,3,3))

	      DO nj=1,3,1
	         theta(1,nj)=dt*avel(1,nj)+theta(0,nj)
	      ENDDO
	     	      		
	      R0(1,1)=ornt(0,1)
	      R0(1,2)=ornt(0,4)
	      R0(1,3)=ornt(0,7)
	      R0(2,1)=ornt(0,2)
	      R0(2,2)=ornt(0,5)
	      R0(2,3)=ornt(0,8)
	      R0(3,1)=ornt(0,3)
	      R0(3,2)=ornt(0,6)
	      R0(3,3)=ornt(0,9)
	      
	      Rz(1,1)=cos(theta(1,3)*deg2rad)
	      Rz(1,2)=-sin(theta(1,3)*deg2rad)
	      Rz(1,3)=0.0d0
	      Rz(2,1)=sin(theta(1,3)*deg2rad)
	      Rz(2,2)=cos(theta(1,3)*deg2rad)
	      Rz(2,3)=0.0d0
	      Rz(3,1)=0.0d0
	      Rz(3,2)=0.0d0
	      Rz(3,3)=1.0d0
	      
	      Ry(1,1)=cos(theta(1,2)*deg2rad)
	      Ry(1,2)=0.0d0
	      Ry(1,3)=sin(theta(1,2)*deg2rad)
	      Ry(2,1)=0.0d0
	      Ry(2,2)=1.0d0
	      Ry(2,3)=0.0d0
	      Ry(3,1)=-sin(theta(1,2)*deg2rad)
	      Ry(3,2)=0.0d0
	      Ry(3,3)=cos(theta(1,2)*deg2rad)
	      
	      Rx(1,1)=1.0d0
	      Rx(1,2)=0.0d0
	      Rx(1,3)=0.0d0
	      Rx(2,1)=0.0d0
	      Rx(2,2)=cos(theta(1,1)*deg2rad)
	      Rx(2,3)=-sin(theta(1,1)*deg2rad)
	      Rx(3,1)=0.0d0
	      Rx(3,2)=sin(theta(1,1)*deg2rad)
	      Rx(3,3)=cos(theta(1,1)*deg2rad)
	      
	      CALL atimesb(Ry,Rx,C1,3)
	      CALL atimesb(Rz,C1,C2,3)
	      CALL atimesb(R0,C2,R,3)
	      
	      ornt(1,1)=R(1,1)
	      ornt(1,2)=R(2,1)
	      ornt(1,3)=R(3,1)
	      ornt(1,4)=R(1,2)
	      ornt(1,5)=R(2,2)
	      ornt(1,6)=R(3,2)
	      ornt(1,7)=R(1,3)
	      ornt(1,8)=R(2,3)
	      ornt(1,9)=R(3,3)
	      	      
	      time=time+dt
	      
	      nc=0
	      DO nj=1,3,1
	         avel(0,nj)=avel(1,nj)
		 theta(0,nj)=theta(1,nj)
		 DO nv=1,3,1
		    nc=nc+1
		    ornt(0,nc)=ornt(1,nc)
		 ENDDO
	      ENDDO
	      
	   ENDDO !ndt	
	   
	   !update rotation array for time=T	 
	   nc1=0
	   DO nj=1,3,1
	      DO nv=1,3,1
	         nc1=nc1+1
	         rotation(nb,out_time_steps_reqd,nc1)=ornt(1,nc1)  !1 -9 orientaion of principal axes	
	      ENDDO
	      nc2=nj+9
	      rotation(nb,out_time_steps_reqd,nc2)=theta(1,nj)     !10 - 12 theta
	      nc2=nj+12
	      rotation(nb,out_time_steps_reqd,nc2)=aveli(1,nj)	   !13 - 15 angular velocity						
	      nc2=nj+15
	      rotation(nb,out_time_steps_reqd,nc2)=aacci(1,nj)     !15 - 18 angular accelaration
	   ENDDO	   

	ENDDO !nb
        
	
	RETURN
	END


