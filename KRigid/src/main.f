        PROGRAM main

!***********************************************************************
!        Written by Kumar Mithraratne
!        © Auckland Bioengineering Institute
!        V1.1 July 2013
!***********************************************************************
        
        IMPLICIT NONE
        
        INCLUDE 'rigidbody1.cmn'
        
        REAL*8 density(nbm),rog(nbm,3),volume(nbm)
        REAL*8 posni(nbm,3),lveli(nbm,3),lacci(nbm,3)
        REAL*8 ornti(nbm,3,3),aveli(nbm,3),aacci(nbm,3)
        REAL*8 eforce(nbm,3),etorque(nbm,3),bforce(nbm,3)
        REAL*8 time_initial,time_final,dt,output_interval
        REAL*8 translation(nbm,0:time_steps,9)
        REAL*8 rotation(nbm,0:time_steps,18)                
        !INTEGER 
        CHARACTER*100 rbname(nbm)
        
                
        CALL inertial_properties(rbname,volume,density,rog)

        CALL initial_conditions(posni,lveli,lacci,ornti,aveli,aacci)
        
        CALL external_forces(eforce,etorque)
        
        CALL body_forces(bforce)        
        
        CALL time_parameters(time_initial,time_final,dt,
     '       output_interval)
        
        CALL translational_motion(volume,density,
     '       posni,lveli,lacci,eforce,bforce,
     '       time_initial,time_final,dt,output_interval,translation)
        
        CALL rotational_motion(volume,density,rog,
     '       ornti,aveli,aacci,etorque,
     '       time_initial,time_final,dt,output_interval,rotation)
                
        CALL output(rbname,translation,rotation,output_interval)        
        
                
        STOP
        END
