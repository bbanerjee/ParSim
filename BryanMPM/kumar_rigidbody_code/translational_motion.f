	SUBROUTINE translational_motion(volume,density,
     '             posni,lveli,lacci,eforce,bforce,
     '             time_initial,time_final,dt,output_interval,
     '             translation)
	
		
!***********************************************************************
!	Written by Kumar Mithraratne
!	© Auckland Bioengineering Institute
!	V1.1 July 2013
!***********************************************************************
	
	IMPLICIT NONE
	
	INCLUDE 'rigidbody1.cmn'
	
	
	REAL*8 volume(nbm),density(nbm),
     '         posni(nbm,3),lveli(nbm,3),lacci(nbm,3),
     '         eforce(nbm,3),bforce(nbm,3),
     '         time_initial,time_final,dt,output_interval,
     '         translation(nbm,0:time_steps,9)
	!INTEGER
	
	REAL*8 mass(nbm),posn(0:1,3),lvel(0:1,3),lacc(0:1,3),
     '         time,out_time
	INTEGER nb,nt,nv,k,nj,comp_time_steps,ndt,
     '	        time_counter
     
             
	!check whether enough memory has been allocated for translation(-,-,-)
	out_time_steps_reqd=int((time_final-time_initial+0.000001d0)/
     '                          output_interval) !out_time_steps_reqd is set here
	IF (out_time_steps_reqd.GT.time_steps) THEN
           WRITE(6,'(" Increase time_steps in rigidbody1.cmn to",i6,
     '           " !!!")')out_time_steps_reqd
           STOP
	ENDIF     

	!initialise arrays
	DO nb=1,nbm,1
	   mass(nb)=0.0d0
	   DO nt=0,time_steps,1
	      DO nv=1,9,1
	         translation(nb,nt,nv)=0.0d0
	      ENDDO
	   ENDDO
	ENDDO
        DO k=0,1,1
	   DO nj=1,3,1
	      posn(k,nj)=0.0d0
	      lvel(k,nj)=0.0d0
	      lacc(k,nj)=0.0d0
	   ENDDO
	ENDDO   

	!update translation array for time=0
	DO nb=1,nbt,1
	   mass(nb)=volume(nb)*density(nb)
	   DO nv=1,3,1
	      nj=nv
	      translation(nb,0,nv)=posni(nb,nj)  !1 - 3 position
	   ENDDO
	   DO nv=4,6,1
	      nj=nv-3
	      translation(nb,0,nv)=lveli(nb,nj)  !4 - 6 linear velocity
	   ENDDO
	   DO nv=7,9,1
	      nj=nv-6
	      translation(nb,0,nv)=lacci(nb,nj)  !7 - 9 linear accelaration
	   ENDDO	
	ENDDO
		
	comp_time_steps=int((time_final-time_initial)/dt+dt/2.0d0)
	
	!DO nb=1,nbt,1
	DO nb=1,1,1
	   time=0.0d0   
	   out_time=0.0d0
	   time_counter=0 

	   DO nj=1,3,1
	      posn(0,nj)=posni(nb,nj)
	      lvel(0,nj)=lveli(nb,nj)	      
	      lacc(0,nj)=lacci(nb,nj)
	      posn(1,nj)=0.0d0
	      lvel(1,nj)=0.0d0	      
	      lacc(1,nj)=0.0d0
	   ENDDO   
	   
	   DO ndt=1,comp_time_steps,1

	      IF(abs(time-out_time).LT.dt/2.0d0) THEN !update translation array	
		 IF (time.NE.0.0d0) THEN
		    DO nj=1,3,1
		       nv=nj
		       translation(nb,time_counter,nv)=posn(1,nj)
		       nv=nj+3
		       translation(nb,time_counter,nv)=lvel(1,nj)
		    ENDDO
		 ENDIF
		 time_counter=time_counter+1		 
	         out_time=out_time+output_interval
	      ENDIF
                
	      DO nj=1,3,1
	         lvel(1,nj)=(dt/mass(nb))*
     '                      (eforce(nb,nj)+mass(nb)*bforce(nb,nj))+   
     '                      lvel(0,nj)	   	   
	         posn(1,nj)=dt*lvel(1,nj)+posn(0,nj)
	      ENDDO
	     	      		
	      time=time+dt
	      DO nj=1,3,1
	         lvel(0,nj)=lvel(1,nj)
		 posn(0,nj)=posn(1,nj)
	      ENDDO
	      
	   ENDDO !ndt	
	
	   !update translation array for time=T 
           DO nj=1,3,1
	      nv=nj
	      translation(nb,out_time_steps_reqd,nv)=posn(1,nj)
	      nv=nj+3
	      translation(nb,out_time_steps_reqd,nv)=lvel(1,nj)
	   ENDDO
	   	   
	ENDDO !nb
	
	
	RETURN
	END


